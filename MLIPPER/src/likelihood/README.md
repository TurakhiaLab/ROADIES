# Likelihood Computation

This directory contains the GPU calculations that answer two questions:

1. How well does the current tree explain the alignment?
2. How well would a query fit on each candidate edge?

Both start by building reusable likelihood summaries ([CLVs](#core-concepts)). The code then combines those summaries either at the root for one whole-tree score or around candidate edges for placement scores.

## The Short Version

```text
observed tip characters
          |
          v
build CLVs with NodeOps
          |
          +--> combine at root --------> current-tree log likelihood
          |
          +--> combine around an edge
                 + query CLV ----------> placement log likelihood
```

The CUDA names follow the equivalent libpll operations where practical, but a reader does not need to know libpll before following the workflow below.

## Core Concepts

- A [CLV](../tree/README.md#what-is-a-clv) stores likelihoods for every site, rate category, and state at a node or directional edge.
- An [upward CLV](../tree/README.md#why-are-there-upward-and-outside-messages) summarizes descendants in postorder.
- An outside message summarizes the rest of the tree on the parent side of an edge; applying that edge's PMAT produces the downward message expressed at the child.
- A scaler prevents numerical underflow; pattern weights restore multiplicity after repeated alignment columns are compressed.
- A [PMAT](../pmatrix/README.md#mathematical-path) describes how character-state probabilities change across one branch and rate category.

For example, consider the edge connecting internal node `X` to its parent `R`:

```text
                         R
                        / \
                       X   Y
                      / \ / \
                     A  B C  D

Upward contribution at X                 Outside/downward contribution toward X

     A ──┐                                      C ──┐
         ├──> X                                 D ──┴──> Y ──> R
     B ──┘                                                   |
                                                             | outside(X)
                                                             v
                                                             X

     upward(X) summarizes {A, B}             outside/downward path carries {C, D}
```

The two messages describe opposite sides of the same edge:

```text
contribution from {A, B}                         contribution from {C, D}
          |                                                 |
          v                                                 v
      upward(X)                                    outside(X), at R side
          |                                                 |
          |                                      PMAT for edge R--X
          |                                                 |
          |                                                 v
          +-------------------------------------- downward(X), at X side
```

The upward pass visits children before parents, so it can combine `A` and `B` into `upward(X)`. The preorder pass first combines the parent-side and sibling contributions into `outside(X)`, then transports that message through the `R--X` transition matrix to produce `downward(X)`. The outside and downward messages contain evidence from the same taxa, but their state coordinates lie on opposite ends of the branch and their numerical values need not be equal.

A CLV is not itself a tree likelihood. It is an intermediate summary. A final log likelihood appears only after the relevant CLVs, PMATs, equilibrium frequencies, rate weights, scalers, and pattern weights are combined.

All functions consume a non-owning `DeviceTree`. Allocation and traversal schedule preparation happen in `tree/tree_generation_device.cu`.

## Function Guide

### `partials.cuh/.cu`

GPU counterparts of `pll_update_partials`:

- `InitializeTipPartialsKernel`: decode observed tip characters into the [starting CLVs](../tree/README.md#what-is-a-clv).
- `UpdatePartialsUpwardKernel`: execute the postorder [`UP` NodeOps](../tree/README.md#what-is-a-clv-operation) for one tree.
- `UpdatePartialsUpwardLevelKernel`: evaluate independent operations from one traversal level in a batched grid.
- `UpdatePartialsDownwardKernel` and `UpdatePartialsDownwardLevelKernel`: execute [`DOWN_LEFT`/`DOWN_RIGHT` NodeOps](../tree/README.md#what-is-a-clv-operation).
- `BuildTreeEdgeOutsideWarpSiteKernel`, `UpdateTreeUpwardWarpSiteKernel`, and `RefreshTreeChildDownWarpSiteKernel`: run DNA+G4 single-operation primitives during sequential branch optimization.

### `root_loglikelihood.cuh/.cu`

- `likelihood::root::compute_root_loglikelihood()`: combine one root CLV with equilibrium frequencies, rate weights, scalers, and pattern weights. This produces the [current-tree score](#the-short-version) and plays the role of libpll's `pll_compute_root_loglikelihood`.

### `placement_likelihood.cuh/.cu`

- `likelihood::placement::update_midpoint_partials()`: combine the [two sides of an edge](../tree/README.md#why-are-there-upward-and-outside-messages) with a query to build the candidate attachment shown in [the short version](#the-short-version).
- `likelihood::placement::compute_edge_loglikelihoods()`: score a batch of candidate edges and produce their placement log likelihoods, analogous to batched `pll_compute_edge_loglikelihood`.

## GPU Parallelism and Outputs

### CLV updates

The ordinary upward and downward kernels parallelize primarily over alignment sites. One thread owns one site at a time and writes the complete `rate_categories * states` CLV slice and its scaler for that site. The grid uses a grid-stride loop, so the same thread processes later sites separated by the total number of threads.

```text
thread/site input:  child or directional CLVs + branch PMATs + one NodeOp
thread/site output: target_clv[site][rate][state] + target_scaler[site]
```

The serial-schedule kernels give each site thread the complete traversal-ordered NodeOp list. This preserves the dependency between a child operation and a later parent operation. The levelized kernels launch a two-dimensional grid: `blockIdx.y` selects one independent NodeOp from the current traversal level, while `blockIdx.x` and `threadIdx.x` divide the sites. Operations within a level may run concurrently, but the host launches levels in dependency order.

```text
serial schedule:    one thread -> one site -> all NodeOps in order
levelized schedule: grid.y      -> one independent NodeOp
                    thread.x    -> one site for that NodeOp
```

The DNA+G4 single-operation fast path assigns 16 lanes to one site: four rate categories multiplied by four states. Each lane writes one `(rate, state)` CLV value. Within each group of four state lanes, warp shuffle instructions perform a maximum reduction; the resulting maximum determines the numerical-scaling shift for that rate category.

### Root likelihood reduction

Root scoring launches 256-thread blocks. Threads read sites with a grid-stride loop and accumulate weighted per-site log likelihoods in `double` registers. The result is reduced in two stages:

```text
sites
  -> each thread accumulates several sites
  -> warp reductions
  -> one sum per block in block_totals[block]
  -> one final 256-thread block reduces block_totals
  -> total[0] = whole-tree log likelihood
```

The first kernel outputs one `double` per block. The final kernel outputs exactly one `double`, which is copied to the CPU. Rate categories and states are combined inside the per-site calculation; they are not separate grid dimensions.

### Placement midpoint and edge likelihood

The midpoint kernel uses a two-dimensional grid. `grid.y` selects one candidate edge operation, while each x-thread owns one site. It reads the upward and outside CLVs plus the two half-edge PMATs and writes `edge_midpoint_clv[target][site][rate][state]` and the corresponding scaler.

Final edge scoring assigns one 256-thread block to each candidate operation. Threads in that block traverse different sites, accumulate their local weighted log likelihoods in `double`, and use a block reduction to produce one value:

```text
grid.y = candidate operation
block threads = sites belonging to that candidate
thread output = local sum across its site stride
block reduction output = placement_log_likelihood[candidate]
```

There is no cross-candidate reduction because every candidate must retain its own score for ranking.

## Normal Call Paths

Full-tree likelihood (the upper output path in the diagram):

```text
PrepareTreeClvOperations()
  -> UpdateTreeClvs*()
  -> likelihood::root::compute_root_loglikelihood()
```

Here `PrepareTreeClvOperations()` creates the reusable [`NodeOpInfo` schedule](../tree/README.md#what-is-a-clv-operation), and `UpdateTreeClvs*()` executes that schedule to refresh the CLVs.

Placement likelihood (the lower output path in the diagram):

```text
EvaluatePlacementCandidates()
  -> likelihood::placement::update_midpoint_partials()
  -> optimize pendant/proximal lengths
  -> likelihood::placement::compute_edge_loglikelihoods()
```

Do not launch partial kernels directly from workflow code. Prefer the host orchestration functions documented in the [tree function guide](../tree/README.md#function-guide) and [placement call flow](../placement/README.md#placement-call-flow); they maintain traversal and scaler invariants.

## Correctness Invariants

- `DeviceTree.sites`, `states`, and `rate_cats` define every CLV stride.
- The PMAT orientation must agree with `downward_pmat_indexing`.
- A topology or branch-length change invalidates specific CLVs and/or PMATs.
- Tip-preserving D&C paths must not overwrite the precomputed CLVs stored in [virtual boundary tips](../tree/README.md#dc-boundary-tips).
- Likelihood accumulation is double precision even when `fp_t` is float.
- Any new fast path must retain the generic path and be compared against it.
- Reduction order is part of numerical behavior: changing thread/block mapping can change low-order floating-point bits even when the mathematical sum is unchanged.
