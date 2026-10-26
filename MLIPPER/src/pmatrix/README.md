# Transition Matrices

This directory converts a substitution model and branch length into transition probability matrices (PMATs). Likelihood and derivative kernels consume these matrices but do not construct the model eigendecomposition themselves.

## What a PMAT Means

A PMAT answers: if a branch starts in state `i`, what is the probability that it ends in state `j` after evolving along this branch?

```text
parent-side state i
        |
        | branch length t
        | rate category r
        v
child-side state j

PMAT[r][i][j] = probability of i -> j across that branch
```

For DNA there are four states, so one rate category produces a `4 x 4` matrix. With four Gamma rate categories, the same branch has four matrices because a fast category and a slow category experience different effective evolutionary time.

A very short branch produces a matrix close to the identity matrix: retaining the same state is likely. As branch length increases, substitution probabilities spread across the other states. The [CLV operations](../tree/README.md#what-is-a-clv-operation) use these probabilities when moving likelihood evidence across a branch.

## Mathematical Path

```text
frequencies + GTR exchangeability rates
  -> Q                           instantaneous substitution process
  -> eigendecomposition          reusable representation of Q
  -> exp(Q * rate * length)      one PMAT for one branch/rate category
  -> CLV updates                 move evidence across branches
  -> likelihood
```

`Q` describes instantaneous substitution rates and belongs to the model. A PMAT is `Q` evolved for a particular branch length and rate multiplier, so model changes can invalidate every PMAT while one branch-length change invalidates only that branch's PMATs.

## How MLIPPER Builds the GTR Q Matrix

For DNA, index the states as `A`, `C`, `G`, and `T`. Let `pi_i` be the equilibrium frequency of state `i`, and let the six GTR exchangeability parameters satisfy `r_ij = r_ji`.

For two different states, MLIPPER constructs the instantaneous rate from state `i` to state `j` as:

$$
Q_{ij} = r_{ij}\pi_j, \qquad i \ne j
$$

The reverse direction uses the same exchangeability but the destination frequency of the other state:

$$
Q_{ji} = r_{ij}\pi_i
$$

Each diagonal entry is the negative sum of the other entries in its row:

$$
Q_{ii} = -\sum_{j \ne i} Q_{ij}
$$

Therefore every row sums to zero, which is required for a continuous-time Markov rate matrix. Negative diagonal entries describe the rate of leaving a state; non-negative off-diagonal entries describe rates of changing into other states.

The raw GTR parameters do not yet define the scale of a branch-length unit. MLIPPER calculates the equilibrium-weighted mean substitution rate:

$$
\mu = -\sum_i \pi_i Q_{ii}
$$

It then divides every entry by `mu`. After this normalization, branch length `t = 1` means one expected substitution per site under the model.

For example, before normalization the four-state matrix has this structure:

```text
        destination
          A          C          G          T
from A   Q_AA       r_AC*pi_C  r_AG*pi_G  r_AT*pi_T
from C   r_AC*pi_A  Q_CC       r_CG*pi_G  r_CT*pi_T
from G   r_AG*pi_A  r_CG*pi_C  Q_GG       r_GT*pi_T
from T   r_AT*pi_A  r_CT*pi_C  r_GT*pi_G  Q_TT

Q_ii = negative sum of the other entries in row i
```

The implementation of these equations is `build_gtr_q_matrix()` in `tree/tree_generation.cpp`.

## Why Eigenvectors Make PMAT Construction Cheaper

The transition matrix for branch length `t` and rate-category multiplier `r` is a matrix exponential:

$$
P(t,r) = \exp(Qrt)
$$

Computing a general matrix exponential separately for every branch and rate category would repeat expensive matrix work. Instead, MLIPPER decomposes `Q` once:

$$
Q = V\Lambda V^{-1}
$$

`V` contains the eigenvectors and `Lambda` is diagonal with eigenvalues `lambda_k`. Powers and exponentials of a diagonal matrix act on each diagonal value independently, giving:

$$
P(t,r) = V\,\operatorname{diag}\!\left(e^{\lambda_k r t}\right)V^{-1}
$$

```text
model changes
  -> recompute Q
  -> compute V, eigenvalues, and V^-1 once

for each branch length t and rate r
  -> compute exp(lambda_k * r * t) for each eigenvalue
  -> multiply V * diagonal exponential * V^-1
  -> obtain PMAT(t, r)
```

This works especially well because the rate matrix `Q` constructed by GTR is
reversible with respect to its equilibrium frequencies `pi`. That reversibility
is exactly the detailed-balance condition:

$$
\pi_i Q_{ij} = \pi_j Q_{ji}
$$

If `D = diag(pi)`, then the similar matrix `S = D^(1/2) Q D^(-1/2)` is symmetric. A real symmetric matrix has real eigenvalues and orthogonal eigenvectors, so MLIPPER can use LAPACK's stable symmetric eigensolver. It transforms the resulting eigenvectors back to obtain `V` and `V^-1` for the original `Q`.

The eigenvectors describe reusable substitution directions determined by the model. Branch length and Gamma rate do not change those directions; they only change how strongly each direction decays through `exp(lambda_k * r * t)`. This is why one eigendecomposition can be reused for every branch until the substitution rates or equilibrium frequencies change.

For short branches, MLIPPER computes `exp(x)` as `1 + expm1(x)` in the PMAT implementation. `expm1` preserves more precision when `x` is close to zero, which is common for small branch lengths.

## Function Guide

- `pmat.h/.cpp`
  - `EigResult`: owns the eigenvalues, eigenvectors, and inverse eigenvectors from the [reusable decomposition of Q](#why-eigenvectors-make-pmat-construction-cheaper).
  - `gtr_eigendecomp_cpu()`: uses GTR reversibility to convert `Q` into that reusable symmetric eigendecomposition.
  - `pmatrix_from_triple()`: constructs `exp(Q * rate * length)` on the CPU for the branch transition described in [What a PMAT Means](#what-a-pmat-means).
- `pmat_gpu.cuh/.cu`
  - `pmatrix_from_triple_gpu()`: is a GPU-only helper called inside CUDA kernels; it performs the same per-matrix construction for a small state space, currently up to 16 states.
  - `build_all_branch_pmats_gpu()`: is a CPU-callable GPU wrapper that applies the [global rebuild policy](#when-pmats-must-be-rebuilt) to every branch.
  - `build_single_branch_pmat_gpu()`: is a CPU-callable GPU wrapper that applies the [single-branch rebuild policy](#when-pmats-must-be-rebuilt) after one sequential branch update.

## How GPU Work Is Assigned

The GPU code does not use the same thread-to-data mapping for every rebuild. A full-tree or placement rebuild exposes many independent PMATs, whereas a single-branch rebuild has only a few matrices and therefore parallelizes their entries instead.

### All-branch rebuild: one thread per PMAT

`build_all_branch_pmats_gpu()` launches blocks of 256 threads. The flattened work list contains one item for every `(node, rate category)` pair:

```text
item_count = node_count * rate_categories

flat item 0  -> node 0, rate 0 -> one complete states x states PMAT
flat item 1  -> node 0, rate 1 -> one complete states x states PMAT
...
flat item R  -> node 1, rate 0 -> one complete states x states PMAT
```

Each thread converts its flat item index back into `node_id` and `rate_id`, then calls `pmatrix_from_triple_gpu()`. That thread computes all entries of its small matrix serially. For the maintained four-state DNA model, one thread computes 16 output values, while many branch/rate matrices are processed concurrently across the grid.

```text
CUDA grid
  block 0: threads 0..255   -> first 256 branch/rate PMATs
  block 1: threads 0..255   -> next 256 branch/rate PMATs
  ...

inside one thread
  for row in states
    for column in states
      sum over eigenvalue index
```

### Single-branch rebuild: threads cooperate on matrix entries

`build_single_branch_pmat_gpu()` launches one block of 128 threads because only one branch is changing. It first computes the `rate_categories * states` diagonal `expm1` values in shared memory. After synchronization, the block flattens all `(rate category, row, column)` output entries and distributes them across its threads:

```text
element_count = rate_categories * states * states

thread 0 -> output element 0, then 128, 256, ...
thread 1 -> output element 1, then 129, 257, ...
...
```

With four rate categories and four DNA states, there are `4 * 4 * 4 = 64` output entries. The first 64 threads each calculate one PMAT entry, and the remaining threads have no matrix entry to write. This mapping is still preferable to assigning one thread per PMAT here, because a single updated branch provides only four independent PMATs.

### Placement PMATs: one thread per candidate/rate PMAT

Placement launches blocks of 128 threads and uses the many candidate edges to expose parallel work. Pendant PMAT work is flattened by `(candidate operation, rate category)`, while proximal and distal PMAT work is flattened by `(node, rate category)`. Each thread calls `pmatrix_from_triple_gpu()` and computes one complete small PMAT. These matrices live in placement scratch buffers and do not overwrite the reference-tree PMAT array.

The important distinction is:

```text
many independent matrices available -> one thread computes one complete PMAT
only one branch being updated        -> one block divides the PMAT entries among threads
```

## Ownership and Layout

`EigResult` owns CPU `double` vectors. The corresponding `DeviceTree` fields (`d_lambdas`, `d_V`, `d_Vinv`) are non-owning GPU pointers. A PMAT contains `rate_categories * states * states` values per node slot.

The PMAT array uses the child node ID as the branch slot:

```text
parent P
   |
   | branch represented by PMAT slot C
   |
child C
```

In other words, slot `C` represents the branch from child `C` to its parent `P`. The root slot does not represent an ordinary parent edge, but it remains allocated so every node ID has the same fixed-stride address calculation.

For how a PMAT is consumed, see [What is a CLV operation?](../tree/README.md#what-is-a-clv-operation).

## When PMATs Must Be Rebuilt

- substitution rates or equilibrium frequencies changed;
- Gamma alpha changed the category rate multipliers;
- a branch length changed;
- subtree/sector extraction remapped full-tree branches into different local node IDs.

Use the all-branch builder after a global model change or complete topology repack. Use the single-branch builder when one existing-tree branch length changes during [sequential optimization](../placement/README.md#branch-optimization-semantics). Placement uses separate scratch PMATs for its [proximal, distal, and pendant branches](../placement/README.md#the-placement-question) because candidate scoring must not overwrite the reference tree's PMATs.

Rebuilding too little leaves stale transition probabilities and therefore produces a likelihood for the wrong model/tree state. Rebuilding every PMAT after every single-edge update is mathematically safe but discards the intended performance benefit of a local update.

Numerical cleanup in `pmat.cpp` clamps small negative values and normalizes rows. Changes here affect every workflow and require likelihood regression.
