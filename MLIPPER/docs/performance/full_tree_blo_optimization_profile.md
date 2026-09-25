# Full-Tree BLO Optimization Performance Profiling

This profile used gene638 to distinguish CUDA launch overhead from insufficient
GPU parallelism in the RAxML-style sequential full-tree branch-length
optimizer. It was recorded on 2026-07-31.

Status: historical performance evidence. The measurements explain past design
decisions; they are not current performance guarantees or setup instructions.

Terminology note: this report uses the historical name `mid_base` for the
parent-side edge message. Current code calls this the edge-outside message and
stores it in `d_edge_outside_clv`; there is no separate `mid_base` buffer.

## Benchmark configuration

| Property | Value |
|---|---:|
| Dataset | gene638 |
| Tips | 2,096 |
| Nodes | 4,191 |
| Original sites | 1,859 |
| Compressed patterns | 1,215 |
| States | 4 |
| Rate categories | 4 |
| Model | cold `GTR{1/1/1/1/1/1}+FU{0.25/0.25/0.25/0.25}+G4m{1.0}` |
| Precision | double |
| GPU | NVIDIA RTX A6000-class, approximately 84 SMs |
| Newton maximum | 30 iterations per branch |
| Multi-block configuration | 8 blocks per branch |

All profiles below cover one sequential full-tree sweep. Nsight Systems adds
instrumentation overhead, so kernel percentages should be used for attribution
and the unprofiled runs should be used for end-to-end timing.

## GPU kernel time by operation

### Version 1: single-block branch Newton

| Operation | Instances | GPU time | Share |
|---|---:|---:|---:|
| Single-branch sumtable + Newton | 4,190 | 1.0064 s | 66.4% |
| Full downward CLV operation | 8,632 | 0.3932 s | 25.9% |
| Upward CLV operation | 2,599 | 0.0658 s | 4.3% |
| Single-branch PMAT rebuild | 4,190 | 0.0422 s | 2.8% |
| Tip initialization | 6 | 0.0078 s | 0.5% |
| Root likelihood | 4 | 0.0002 s | <0.1% |
| **Total GPU kernel time** |  | **1.5158 s** | **100%** |

Unprofiled one-sweep BLO time was **1.5161 s**, and wall time was **1.94 s**.

### Version 2: eight-block branch Newton, full post-update downward refresh

| Operation | Instances | GPU time | Share |
|---|---:|---:|---:|
| Full downward CLV operation | 8,632 | 0.4013 s | 48.1% |
| Multi-block sumtable + Newton | 4,190 | 0.3162 s | 37.9% |
| Upward CLV operation | 2,599 | 0.0665 s | 8.0% |
| Single-branch PMAT rebuild | 4,190 | 0.0432 s | 5.2% |
| Tip initialization | 6 | 0.0078 s | 0.9% |
| Root likelihood | 4 | 0.0002 s | <0.1% |
| **Total GPU kernel time** |  | **0.8352 s** | **100%** |

Unprofiled one-sweep BLO time was **0.8041 s**, and wall time was **1.21 s**.

The multi-block Newton kernel reduced its own GPU time by **3.18x**, while the
complete BLO sweep improved by about **1.89x**. The full downward refresh then
became the largest operation.

### Version 3: eight-block branch Newton with incremental downward refresh

| Operation | Instances | GPU time | Share |
|---|---:|---:|---:|
| Multi-block sumtable + Newton | 4,190 | 0.3136 s | 45.1% |
| Full pre-optimization downward operation | 4,442 | 0.2297 s | 33.1% |
| Upward CLV operation | 2,599 | 0.0662 s | 9.5% |
| Single-branch PMAT rebuild | 4,190 | 0.0429 s | 6.2% |
| Incremental `mid_base -> child.down` refresh | 2,094 | 0.0345 s | 5.0% |
| Tip initialization | 6 | 0.0078 s | 1.1% |
| Root likelihood | 4 | 0.0002 s | <0.1% |
| **Total GPU kernel time** |  | **0.6948 s** | **100%** |

Unprofiled one-sweep BLO time was **0.6771 s**, and wall time was **1.11 s**.

The incremental update skips the post-update downward message for tip edges and
reapplies only the accepted edge PMAT to cached `mid_base` for internal edges.

### Version 4: warp-site branch Newton with incremental downward refresh

The dedicated warp-site function assigns 16 lanes to each DNA+G4 site, with
one lane per rate/eigen component. gene638 therefore exposes 76 blocks per
branch instead of the previous eight-block grid.

| Operation | Instances | GPU time | Share |
|---|---:|---:|---:|
| Full pre-optimization downward operation | 4,442 | 0.2300 s | 40.9% |
| Warp-site sumtable + Newton | 4,190 | 0.1803 s | 32.1% |
| Upward CLV operation | 2,599 | 0.0664 s | 11.8% |
| Single-branch PMAT rebuild | 4,190 | 0.0429 s | 7.6% |
| Incremental `mid_base -> child.down` refresh | 2,094 | 0.0348 s | 6.2% |
| Tip initialization | 6 | 0.0078 s | 1.4% |
| Root likelihood | 4 | 0.0002 s | <0.1% |
| **Total GPU kernel time** |  | **0.5624 s** | **100%** |

Unprofiled one-sweep BLO time was **0.5513 s**, and wall time was **1.00 s**.

The warp-site function reduces Newton GPU time from `0.3136 s` to `0.1803 s`
(**1.74x**) without changing the one-sweep or converged likelihood.

### Version 5: warp-site Newton and branch-local CLV updates

This version also replaces the small generic tree-operation launches inside
the sequential traversal. Before each branch optimization it constructs only
the required `mid_base`; after accepting the branch it refreshes only the
affected child-down or parent-up message. Each specialized kernel assigns the
16 DNA+G4 state/rate values of a site to cooperating lanes.

| Operation | Instances | GPU time | Share |
|---|---:|---:|---:|
| Warp-site sumtable + Newton | 4,190 | 0.1810 s | 53.3% |
| Initial/final full downward operations | 252 | 0.0514 s | 15.2% |
| Single-branch PMAT rebuild | 4,190 | 0.0418 s | 12.3% |
| Initial/final full upward operations | 504 | 0.0262 s | 7.7% |
| Branch-local `mid_base` construction | 4,190 | 0.0163 s | 4.8% |
| Branch-local parent-up refresh | 2,095 | 0.0084 s | 2.5% |
| Tip initialization | 6 | 0.0078 s | 2.3% |
| Branch-local child-down refresh | 2,094 | 0.0062 s | 1.8% |
| Root likelihood | 4 | 0.0002 s | <0.1% |
| **Total GPU kernel time** |  | **0.3394 s** | **100%** |

Unprofiled one-sweep BLO time was **0.3339 s**, and wall time was **0.77 s**.

The branch-local CLV kernels reduce the one-sweep BLO time from `0.5513 s` to
`0.3339 s` (**1.65x**) without changing its likelihood
(`-353284.028063`).

## CUDA launch measurements

| Version | Regular launches | Cooperative launches | Total launches | CUDA launch API time |
|---|---:|---:|---:|---:|
| Single block | 19,621 | 0 | 19,621 | 1.3927 s |
| Multi-block | 15,431 | 4,190 | 19,621 | 0.7432 s |
| Incremental multi-block | 13,335 | 4,190 | 17,525 | 0.6061 s |
| Incremental warp-site | 13,335 | 4,190 | 17,525 | 0.4801 s |
| Warp-site Newton and CLVs | 13,335 | 4,190 | 17,525 | 0.2723 s |

The CUDA API column is **not pure launch overhead** and must not be added to GPU
kernel time. CUDA launches are asynchronous; when the queue or GPU is busy, a
launch API call can absorb queue backpressure. This is why the API totals can
be comparable to the GPU execution time without producing an equally large
idle interval on the GPU.

To measure actual launch starvation, the interval from the first branch Newton
kernel to the last branch Newton kernel was examined directly:

| Version | Sweep GPU window | Kernel-busy time | Inter-kernel idle gaps | Idle share |
|---|---:|---:|---:|---:|
| Single block | 1,443.1 ms | 1,429.8 ms | 13.3 ms | 0.9% |
| Multi-block | 762.9 ms | 749.3 ms | 13.6 ms | 1.8% |
| Incremental multi-block | 620.9 ms | 609.1 ms | 11.8 ms | 1.9% |
| Incremental warp-site | 488.6 ms | 476.6 ms | 12.0 ms | 2.5% |
| Warp-site Newton and CLVs | 266.1 ms | 253.7 ms | 12.3 ms | 4.6% |

The GPU is therefore not spending a large fraction of the sweep waiting for
the host to submit the next kernel. Removing launch latency alone cannot close
the performance gap to RAxML-NG.

## Converged benchmark comparison

| Method | Final log-likelihood | BLO compute time | Full wall time |
|---|---:|---:|---:|
| Tuned global L-BFGS | -353258.586956 | 3.3417 s | 3.74 s |
| Sequential multi-block, incremental, 8 sweeps | -353258.583554 | 4.5846 s | 5.00 s |
| Sequential warp-site, incremental, 8 sweeps | -353258.583554 | 3.6004 s | 4.04 s |
| Sequential warp-site Newton and CLVs, 8 sweeps | -353258.583554 | 1.7918 s | 2.22 s |
| Sequential warp-site Newton and CLVs, RAxML-parity stop at 6 sweeps | -353258.584266 | 1.4187 s | 1.75 s |
| RAxML-NG branch optimization | -353258.584266 | included in wall | 0.918 s |

The L-BFGS compute time is `derivative_s + line_search_s` (`0.0665 + 3.2752`
seconds), not `line_search_s` alone.

## Conclusion

The sequential GPU implementation is **not primarily limited by launch
starvation**:

- The warp-site layout exposes 76 blocks per branch on gene638 instead of the
  original eight-block Newton or approximately five-block generic CLV grids.
- The GPU is executing kernels for 95.4% of the measured sequential sweep
  window; inter-kernel gaps account for only 12.3 ms.
- Specialized branch-local CLV work now accounts for only 9.1% of GPU kernel
  time (`mid_base`, upward refresh, and child-down refresh combined).
- Newton is the largest remaining operation at 53.3%, followed by the
  initial/final whole-tree CLV passes and the accepted-branch PMAT rebuilds.
- RAxML-NG's CPU implementation is well matched to these small sequential
  coordinate updates.

The specialized CLV path reduces the converged eight-sweep BLO stage from
`3.6004 s` to `1.7918 s` and is now faster than the `3.3417 s` global L-BFGS
baseline, while preserving its final likelihood. The next local targets are
the Newton kernel and eliminating or fusing the per-accepted-branch PMAT
rebuild; reducing launch count alone has limited upside at the measured 4.6%
idle share.

## Sweep convergence and RAxML-NG comparison

The three `after brlen` values in the RAxML-NG 2.0.2 debug log are three BLO
**calls**, not three individual full-tree sweeps. The calls continue from the
branch lengths accepted by the previous call. MLIPPER's per-sweep trace maps
to those values as follows:

| MLIPPER sweep | Log-likelihood | Improvement | RAxML-NG milestone |
|---:|---:|---:|---|
| 1 | -353284.028063470 | 369.786252 |  |
| 2 | -353260.095370842 | 23.932693 |  |
| 3 | -353258.873279094 | 1.222092 | Initial BLO: -353258.873273 |
| 4 | -353258.594410224 | 0.278869 |  |
| 5 | -353258.585367932 | 0.009042 | First model-loop BLO: -353258.585368 |
| 6 | -353258.584266213 | 0.001102 | Final RAxML-NG: -353258.584266 |
| 7 | -353258.583804107 | 0.000462 |  |
| 8 | -353258.583554002 | 0.000250 |  |

RAxML-NG parity therefore requires six MLIPPER sweeps on this input. Eight
sweeps are not required for parity and improve the likelihood by only another
`0.000712`. A configurable per-sweep stopping check was added through
`MLIPPER_BRANCH_SWEEP_LH_EPSILON`; `0.005` stops after sweep 6 and reduces BLO
time from `1.7918 s` to `1.4187 s`, with final likelihood
`-353258.584266` and full wall time `1.75 s`.

Reducing the maximum Newton iterations from 30 to 8 produces the same
six-sweep likelihood and essentially the same runtime (`1.4275 s` versus
`1.4369 s`). The branch kernel already exits when its derivative or step
converges, so the nominal iteration cap is not the main cost. An experimental
split of sumtable construction and Newton iteration was also slower on
gene638: its best tested one-sweep time was `0.3397 s`, compared with roughly
`0.3255 s` for the fused kernel. The additional launch and global sumtable
round trip outweighed the smaller Newton synchronization grid.

## Profile artifacts

- `/tmp/mlipper_gene638_raxmlseq_s1_profile.nsys-rep`
- `/tmp/mlipper_gene638_multiblock_s1_profile.nsys-rep`
- `/tmp/mlipper_gene638_incremental_s1_profile.nsys-rep`
- `/tmp/mlipper_gene638_warpsite_s1_profile.nsys-rep`
- `/tmp/mlipper_gene638_warpclv_s1_profile.nsys-rep`
- `/tmp/mlipper_gene638_raxmlseq_s1_stats.txt`
- `/tmp/mlipper_gene638_multiblock_s1_stats.txt`
- `/tmp/mlipper_gene638_incremental_s1_stats.txt`
- `/tmp/mlipper_gene638_warpsite_s1_stats.txt`
- `/tmp/mlipper_gene638_warpclv_s1_stats.txt`

The `/tmp` artifacts are machine-local and should be archived separately if
the raw traces need to be retained long term.
