# Placement and Branch Derivatives

This directory scores query attachments and implements the derivative-based branch-length updates shared by placement and existing-tree optimization.

## The placement question

For one query and one candidate edge, placement asks which attachment lengths give the highest likelihood:

```text
                         original candidate edge

parent-side endpoint P ------------- child-side endpoint D
                         length L

                                  becomes

parent-side endpoint P ---- proximal ---- J ---- distal ---- child-side endpoint D
                                          |
                                       pendant
                                          |
                                       query Q
```

`J` is the proposed attachment node. The proximal edge points toward the parent/root side, the distal edge points toward the original child/subtree side, and the pendant edge leads to the query.

During placement scoring, the original candidate-edge length `L` is divided into proximal and distal pieces. MLIPPER optimizes the pendant length and the split position; the remaining split length is derived from `distal = L - proximal`. It does not independently optimize all three lengths, because doing so would also change the length of the existing reference-tree edge.

The parent-side and child-side evidence are the [outside and upward CLV messages](../tree/README.md#why-are-there-upward-and-outside-messages). These messages let every candidate attachment be scored without rebuilding the whole reference tree.

## Newton–Raphson Branch Updates

For one branch length `x`, Newton–Raphson uses the slope and curvature of the log likelihood:

```text
current length:       x
first derivative:     dL/dx
second derivative:    d2L/dx2

Newton proposal:      x_new = x - (dL/dx) / (d2L/dx2)
```

The first derivative says which direction improves the likelihood.

The second derivative controls how large the proposed step is.

A derivative near zero means the current length is near a stationary point.

The raw formula is not safe enough by itself. MLIPPER follows an EPA-ng-style safeguarded update: it maintains legal lower and upper bounds, limits the maximum displacement, uses a bracket midpoint when the Newton proposal is unsafe, stops at convergence or a hard bound, and rejects non-finite numerical states. The likelihood is evaluated after a complete placement smoothing pass, and a pass that makes the score worse is rolled back.

## Placement Branch-Update Order

Each candidate edge starts with the existing edge split at its midpoint and a default pendant length. A smoothing pass then updates branches in this order:

```text
baseline candidate
  |
  +--> initialize proximal = L / 2
  |    derive distal = L - proximal
  |    initialize pendant
  |    build PMATs and compute baseline likelihood
  |
  +--> 1. Newton-optimize pendant
  |       |
  |       +--> rebuild query/pendant PMAT
  |
  +--> 2. Newton-optimize the proximal split coordinate
  |       |
  |       +--> rebuild proximal PMAT
  |       +--> derive distal = L - proximal
  |       +--> rebuild distal PMAT
  |
  +--> 3. score the complete candidate attachment
          |
          +--> improved or equal: keep lengths for the next pass
          +--> worse/non-finite:   restore the previous best lengths
```

Pendant is updated first, so the subsequent split update sees the new query-side branch. Proximal is updated second, and distal changes as its complement on the original edge. This is a sequential smoothing order, not a simultaneous three-variable Newton step. MLIPPER repeats the pass according to the placement tuning settings, then ranks candidate edges by their retained likelihoods.

## Files

- `placement.cuh/.cu`: placement candidate evaluation, ranking, result collection, and host-level full-tree/selected-edge branch optimizers.
- `derivative.cuh/.cu`: CUDA derivative kernels and single-edge Newton update primitives.

The files currently remain together because placement and tree-edge optimization reuse substantial derivative and PMAT machinery. Workflow code should call the host functions in the headers, not individual internal kernels.

## Placement Result

`RawPlacementResult` contains the best edge and optimized distal/pendant branch lengths. Placement-only mode can serialize ranked placements as jplace output, but output formatting belongs to [`io/`](../io/README.md), not to this scoring component.

`target_id` is a node ID representing the edge from that node to its current parent. Existing full-tree IDs are not renumbered by a traversal rebuild, but the edge represented by `(target_id, parent[target_id])` can change after a topology edit. Sector-local IDs also belong to a different namespace from full-tree IDs; see the [node-ID guide](../tree/README.md#node-ids-traversal-order-and-stable-labels).

## Function Guide

- `EvaluatePlacementCandidates()`: run the [placement branch-update order](#placement-branch-update-order) for each candidate [`NodeOpInfo` edge operation](../tree/README.md#what-is-a-clv-operation), optionally refine the local child side, and return the best/ranked placements.
- `RunAcceptedFullTreeSequentialBranchLengthOptimization()`: perform the [sequential update policy](#branch-optimization-semantics). Each accepted edge refreshes its PMAT and affected CLVs before the next edge.
- `RunSelectedTreeEdgeJacobiBranchLengthOptimization()`: perform the [Jacobi update policy](#branch-optimization-semantics) on a selected edge set.
- `OptimizeSingleTreeEdgeFromCurrentClvs()`: generic single-edge Newton update from the directional CLVs surrounding an existing edge.
- `OptimizeSingleTreeEdgeFromCurrentClvsWarpSite()`: specialized DNA+G4 version used when the device/model shape supports it.
- `LikelihoodDerivativePendantKernel` and `LikelihoodDerivativeProximalKernel`: calculate the derivatives used by the [Newton–Raphson updates](#newtonraphson-branch-updates) for the two independently updated placement coordinates.

## Placement Call Flow

```text
tree::EvaluatePlacementQueries()
  -> EvaluatePlacementCandidates()
       -> build midpoint partials
       -> optimize pendant/proximal lengths
       -> compute candidate likelihoods
       -> rank candidates
  -> optionally commit the selected query into the CPU tree
```

`PlacementScratchOverride` allows D&C/Local SPR workspaces to provide already allocated scratch buffers. Every pointer is non-owning, and each capacity field must cover the corresponding operation or node count.

## GPU Parallelism During Placement

Placement exposes parallelism in two directions: independent candidate edges and independent alignment sites within each candidate. The PMAT mapping is documented in [How GPU Work Is Assigned](../pmatrix/README.md#how-gpu-work-is-assigned); the likelihood mapping and output buffers are documented in [GPU Parallelism and Outputs](../likelihood/README.md#gpu-parallelism-and-outputs).

The high-level assignment is:

| Stage                               | GPU assignment                                                  | Output                                        |
| ----------------------------------- | --------------------------------------------------------------- | --------------------------------------------- |
| Initialize pendant/proximal lengths | One thread per node                                             | Initial branch-length arrays                  |
| Build pendant PMATs                 | One thread per`(candidate operation, rate category)`          | One scratch PMAT per candidate/rate pair      |
| Build proximal/distal PMATs         | One thread per`(node, rate category)`                         | Node-indexed scratch PMAT arrays              |
| Build midpoint CLVs                 | `grid.y` selects a candidate and each x-thread selects a site | Candidate-edge midpoint CLVs and scalers      |
| Pendant Newton update               | One block per candidate; threads divide that candidate's sites  | One optimized pendant length per target edge  |
| Proximal Newton update              | One block per candidate; threads divide that candidate's sites  | One optimized proximal length per target edge |
| Final edge score                    | One block per candidate; threads divide sites                   | One log likelihood per candidate              |
| Keep best state                     | One thread per candidate operation                              | Accepted/restored lengths and active flags    |

### Derivative reduction inside one candidate block

Pendant and proximal Newton kernels build one sumtable per candidate. Within that candidate's block, each thread accumulates first and second derivatives over a strided subset of sites. Both derivative values follow the same reduction path so they describe the same site set:

```text
per-site derivative contributions
  -> each thread accumulates local df and ddf
  -> shuffle reduction inside each warp
  -> one (df, ddf) pair per warp in shared memory
  -> first warp reduces the warp pairs
  -> thread 0 updates the candidate's shared Newton state
  -> repeat until convergence or iteration limit
  -> write new_branch_length[target_id]
```

Candidate blocks do not reduce with one another. They optimize independent hypothetical attachments, and the CPU ranks their final per-candidate likelihood outputs afterward.

### Existing-tree branch optimization

The selected-edge Jacobi kernel launches one block per tree node and inactive blocks exit early; the threads in an active block divide that edge's sites and reduce its derivatives to one proposed branch length. Sequential optimization handles one edge at a time so the next edge sees accepted PMAT and CLV updates.

For a large single-edge problem, the cooperative multi-block path divides sites across several 256-thread blocks. Each block first reduces to `partial_gradient[block]` and `partial_hessian[block]`; after a grid-wide barrier, grid thread 0 sums those block outputs and performs the Newton update. Another grid-wide barrier publishes the new branch value before the next iteration. The DNA+G4 fast path instead assigns 16 lanes to each site, one lane for every `(rate category, eigen component)` pair, before following the same block-then-grid reduction structure.

## Branch Optimization Semantics

- Sequential/Gauss-Seidel: later edges see earlier accepted updates. Used for full-tree final branch optimization.
- Jacobi: selected edges are optimized from one common state before the group is committed. Used for local or sector refinement.
- Branch bounds and tolerances come from `optimize/optimization_types.hpp`.
- A branch update is not complete until its PMAT and affected CLVs are current.

## Common Failure Modes

- Reusing a candidate edge solely because its `target_id` still exists after a topology edit; the node may now have a different parent, and sector-local IDs require explicit mapping.
- Using stale midpoint/downward CLVs after changing a branch or topology.
- Passing scratch buffers sized for fewer operations than `num_ops`.
- Comparing a local candidate likelihood as if it were a full-tree acceptance audit.
- Adding a new branch bound locally instead of using the shared optimization policy.
