# RNASim 10K cross-dataset experiment

Status: historical experiment record. Dataset paths and commands below are not
part of the current MLIPPER interface. Start with the root `README.md`; use this
document only to understand the recorded cross-dataset conclusions.

## Dataset selection

The search under `/data/swalia/dipper` identified RNASim 10K as the useful
independent dataset near 10K taxa:

- alignment: `/data/swalia/dipper/DIPPER_data/RNASim/dataset/10K/data.aligned.fa.gz`
- reference: `/data/swalia/dipper/DIPPER_data/RNASim/dataset/10K/data.treefile.nwk`
- 10,000 taxa, 8,835 sites, 82.29% gaps

This is a substantially different regime from the existing AliSim 10K case.
The `10000` ALISIM validation files contain 10,000 sites rather than 10,000
taxa. The AliSim `iterative/10000` alignment did not have an unambiguous paired
reference tree, so it was not used for an nRF experiment.

## Baselines

All RF values below were recomputed with RAxML-NG against the RNASim reference.
The current gap-fixed DIPPER binary was used to rebuild the DIPPER trees.

| Tree | RF | nRF |
|---|---:|---:|
| Old Sumit DIPPER conventional NJ | 4,916 | 0.245874 |
| Gap-fixed DIPPER conventional NJ | 3,280 | 0.164049 |
| Gap-fixed DIPPER strict placement (`-m 1 -K -1`) | 4,974 | 0.248775 |
| FastTree GTR+CAT final | 2,152 | 0.107632 |

The gap fix therefore repairs conventional NJ strongly on this highly gapped
alignment, but does not repair strict placement routing. The useful MLIPPER
starting point is the gap-fixed conventional NJ tree, not the placement tree.

FastTree was run with `-nt -gtr -pseudo -nosupport`. It required 53 ME-NNI
rounds, two radius-10 ME-SPR rounds, and seven ML-NNI rounds. Wall time was
36m33s. Its phase checkpoints show:

| FastTree phase | nRF |
|---|---:|
| Initial NJ | 0.159848 |
| After ME refinement | approximately 0.133 |
| Best/final ML-NNI | 0.107632 |

## Common likelihood model

RAxML-NG estimated a common GTR+Gamma4 model on the FastTree topology:

```text
GTR{0.996763/1.083480/1.126987/0.875124/1.243605/1.000000}
+FU{0.251128/0.244384/0.281239/0.223248}+G4m{0.881405}
```

The gap-fixed DIPPER NJ topology was branch-optimized under this fixed model.
It contained 469 near-zero branches in 423 weak components.

## MLIPPER D&C refinement

The blind refinement used only the current tree and alignment, with no truth or
FastTree guidance:

- weak-component scheduling at branch length `1e-5`
- streaming boundary context and D&C tip budget 300
- exact local NNI plus local SPR
- adaptive radii 2, 3, and 4
- top K 32 and up to three local rounds
- backbone moves enabled inside the extracted sector

All accepted components in these runs found moves at radius 2. The first ten
attempts processed four anchored components. A continuation of 30 attempts
processed ten more anchored components.

| Tree | RF | nRF | Common-model logL |
|---|---:|---:|---:|
| Gap-fixed DIPPER NJ start | 3,280 | 0.164049 | -7,495,060.205167 |
| MLIPPER after first four components | 3,202 | 0.160148 | -7,494,407.082962 |
| MLIPPER after ten additional components | 2,962 | 0.148144 | -7,493,055.355539 |
| FastTree final | 2,152 | 0.107632 | -7,488,265.765531 |

The 14 accepted components improve every recorded RF checkpoint monotonically.
Overall, MLIPPER recovers 318 RF and improves the audited full-tree likelihood
by 2,004.849628. It does not yet close the remaining 810-RF gap to FastTree.

## Streaming whole-tree ML-NNI coverage

To test whether the remaining gap was caused by weak-edge targeting, a second
workflow scheduled boundary-aware 300-tip sectors over the whole tree. Each
sector enumerated its internal NNI edges, retained top K 32 candidates for exact
ML validation, and accepted only positive-gain moves. Sector selection used the
center of the largest connected component of uncovered internal edges.

| Tree | RF | nRF | Common-model logL |
|---|---:|---:|---:|
| Gap-fixed DIPPER NJ start | 3,280 | 0.164049 | -7,495,060.205167 |
| Whole-tree ML-NNI checkpoint 1 | 2,434 | 0.121737 | -7,490,623.246343 |
| Whole-tree ML-NNI checkpoint 2 | 2,270 | 0.113534 | -7,489,742.731622 |
| FastTree final | 2,152 | 0.107632 | -7,488,265.765531 |

The two checkpoints recover 1,010 RF relative to DIPPER and leave a 118-RF
gap to FastTree. Both RF and likelihood improve, so the result is not an RF-only
effect. This is substantially stronger than short/weak-component scheduling
(RF 2,962), and identifies global search coverage as the dominant limitation of
that workflow.

The coverage manifest proves that sector edges entered candidate enumeration;
it does not prove exact rescoring of both NNI alternatives for every edge,
because exact validation is limited to top K 32 per sector. The checkpoints
therefore demonstrate effective global candidate coverage, not exhaustive ML
NNI convergence.

### Is every sector required?

Per-sector RF auditing shows that useful moves are distributed across the full
first sweep. The first 30 sectors recover only 442 of the final 846 RF units
(52%), the first 60 recover 706 (84%), and sectors 61--100 recover another 140.
Stopping the first sweep early would therefore discard substantial signal.

The second sweep has much lower marginal return. Its first 80 sectors recover
154 of 164 RF units, while the final 20 recover only 10. Nevertheless, a simple
rolling early-stop is not reliable: sectors 81--90 worsen RF by two, whereas
sectors 91--100 subsequently recover 12. Sector likelihood gain correlates with
RF improvement (Pearson r=0.71 in the first sweep and r=0.60 in the second), but
does not identify every beneficial sector.

The resulting policy should be:

- require one coverage-complete first sweep;
- prioritize later sectors by a truth-independent likelihood/support proxy;
- apply a time or candidate-validation budget to later sweeps instead of
  requiring full coverage every time;
- retain a final zero-move coverage audit only when strict ML-NNI convergence
  is required.

The complete 200-sector trajectory is recorded in
`tmp/rnasim10k/cross-dataset-v1/streaming-whole-tree-sector-trajectory.csv`.

### Standardized large-sector sweep

The validated fast preset uses fixed split ownership, adaptive anchors, 300
owned edges per sector, a 700-tip boundary-aware subtree, top K 2048, radius-2
exact NNI, and speculative conflict-free validation waves. On RNASim 10K, one
coverage-complete 41-sector sweep produced:

| Tree | RF | nRF | Common-model logL | Wall time |
|---|---:|---:|---:|---:|
| C++-only MLIPPER whole-tree NNI | 2,536 | 0.126838 | -7,490,878.902231 | 166.71 s |

The split scheduler and direct ownership path are now callable from
`MLIPPER_session_main --whole-tree-nni`. The binary performs its own CUDA
context recycle every three sectors while retaining the CPU tree, alignment,
model, split ownership, and coverage state. The complete sweep therefore needs
no Python scheduler, intermediate tree checkpoint, or process restart. Its RF
differs by four from the earlier Python-scheduled result (RF 2,532) because the
C++ scheduler uses a different deterministic sector order.

The production NNI path no longer contains the rejected ME/RELL/beam,
negative-gain bridge, approximate-ranking, or query-group SPR experiments.
Full-tree oracle, placement-audit, legacy-routing validation, and weak-component
experiment drivers were also removed after the standardized path achieved
byte-identical seven-sector output before and after cleanup.

```bash
./MLIPPER_session_main \
  --tree-alignment data.aligned.fa \
  --tree dipper.nwk \
  --best-model model.raxml.bestModel \
  --whole-tree-nni \
  --whole-tree-nni-validation-devices 0,1 \
  --gpu-id 0 \
  --write-tree refined.nwk
```

### Parallel exact-NNI validation

Exact validation can be split into node-disjoint NNI waves with
`MLIPPER_LOCAL_SPR_PARALLEL_VALIDATION=1`.  Validation lanes are assigned with
`MLIPPER_LOCAL_SPR_DEVICES` (for example, `0,1,2,3`). Every candidate in
a wave is evaluated against the same base tree.  Individually positive moves
are then combined and accepted only after a full Gamma4 likelihood audit of the
combined topology.  This mode currently requires
`MLIPPER_LOCAL_SPR_SKIP_SYMMETRIC_BASELINE=1`, beam width 1, no RELL gate, and no
negative-gain bridge moves.

Two RNASim 10K 300-tip sectors reproduced the sequential topology exactly:

| Anchor | Sequential validation / wall | Two-GPU validation / wall | Pairwise RF |
|---|---:|---:|---:|
| `L423468` | 1.83 s / 6.06 s | 1.62 s / 5.15 s | 0 |
| `L1006047` | 1.50 s / 5.61 s | 1.28 s / 4.79 s | 0 |

The rebuilt-versus-validated likelihood discrepancy was at most
`9.4e-10`.  On the first sector, a one-GPU wave took 2.33 s, confirming that
the second device, rather than wave batching alone, provides the validation
speedup.

### Large-core fixed-ownership NNI optimization

The corrected fixed-ownership sweep was first run with 100 owned edges, 300
selected tips, and six anchors per sector. It covered all 9,997 starting
internal edges in 125 sectors, took 442 seconds, and reached RF2534
(`nRF=0.126738`). Halo edges supplied boundary context but could not propose
NNIs.

Increasing the sector to 300 owned edges and 700 selected tips reduced repeated
boundary construction and likelihood audits. Since accepted NNIs can fragment
later starting cores, ownership and order remain fixed while 12 anchor tips are
recomputed from the current connected components before each sector. With a
persistent-session restart every three sectors, the complete 41-sector sweep
took 318.03 seconds and reached RF2520 (`nRF=0.126038`) after 1,047 accepted
NNIs. This is 28.0% faster and 14 RF better than the 100-edge schedule. The
validated settings are available as `--fast-nni-preset`.

Larger 500-edge/1,100-tip sectors were rejected: exact validation exhausted a
48 GB GPU on the second sector. Increasing the boundary site batch from 1,024
to 2,048 did not improve wall time, and 4,096 sites exhausted GPU memory.
Direct-quartet ranking reproduced the same three-sector topology but did not
change wall time, confirming that boundary construction and exact validation,
not candidate ranking, dominate this workflow.

Profiling the first 300-edge sector showed 566.6 ms in candidate ranking,
5,823.9 ms in per-candidate exact validation, and 189.5 ms in the final rebuild.
The standard fast path therefore groups ranking-positive, node-disjoint NNIs
into conflict-free waves and performs one complete Gamma4 smoothing and rebuilt
likelihood audit for each combined wave. A wave is committed only when its
combined likelihood gain is positive.

This speculative-wave path completed all 9,997 edges in 170.38 seconds, accepted
1,102 NNIs, and reached RF2532 (`nRF=0.126638`). Relative to per-candidate exact
validation it is 46.4% faster, with a 12-RF (`0.000600` nRF) quality cost and a
186.624064 lower final log-likelihood. Relative to the original 100-edge fixed
sweep it is 61.5% faster and two RF better. It is now enabled by
`--fast-nni-preset`; `--no-speculative-nni-waves` restores per-candidate exact
validation.

Persistent execution without periodic restart still exhausted a 48 GB GPU in
sector five, so the standard preset retains a restart every three sectors until
the cross-command GPU allocation retention is removed.

Direct-quartet ranking also has an optional multi-device backend enabled with
`MLIPPER_LOCAL_SPR_PARALLEL_RANKING=1`.  It replicates only the active sector's
directional CLV pools, partitions internal-edge ranges across all configured
devices, and merges their candidate lists before top-K selection.  On two PCIe
RTX A6000s, ranking score time improved from about 374 ms to 249 ms (33%), but
refreshing the secondary replica cost 239--359 ms per sector, so total ranking
time increased.  It is therefore intentionally opt-in; it may become useful
for larger sectors or NVLink-connected devices.  Candidate validation remains
the profitable default multi-GPU phase.

### Edge-owned scheduling under a ten-minute budget

The scheduler was extended with disjoint connected edge cores and six
farthest-point anchors per core. The 300-tip sectors retain directional
boundary context, while core ownership and halo coverage are reported
separately. Exact validation used top K 32 on CUDA devices 0 and 1. Separate
MLIPPER processes were required because the persistent stream reproducibly ran
out of GPU memory on the seventh, 631-node compressed sector.

| Stage | Sectors | Wall time | RF | nRF | Accepted NNI |
|---|---:|---:|---:|---:|---:|
| Complete edge-coverage sweep | 96 | 478.76 s | 2,472 | 0.123637 | 1,169 |
| Additional prioritized sectors | 20 | 89.98 s | 2,394 | 0.119736 | 86 |
| Final time-budget sectors | 6 | 26.62 s | 2,392 | 0.119636 | 18 |
| **Combined** | **122** | **595.36 s** | **2,392** | **0.119636** | **1,273** |

This satisfies the measured sub-ten-minute requirement and is within 0.006102
nRF of the older two-sweep checkpoint (`nRF=0.113534`). The final likelihood
audit was `-7,489,975.355830`, an improvement of 5,084.849337 over the starting
DIPPER tree. It remains 232.624208 log-likelihood units below the older
checkpoint, so this is a time-bounded approximation rather than full ML-NNI
convergence.

Static ownership did not reduce the complete sweep to the theoretical 32--41
sectors. Immediate NNI commits fragment the starting cores and leave small
surviving components for tail cleanup. The first sweep therefore still used 96
sectors and 3.56-fold aggregate edge coverage. A frozen proposal wave with
conflict-aware batch commit remains the next scheduler change if stricter
convergence must fit within the same budget.

The faster retained schedule repartitions after 60 sectors, then runs 50
sectors on the updated topology. A persistent MLIPPER process is proactively
restarted every four sectors: this reuses initialization within a block while
avoiding the reproducible CUDA allocator OOM seen in longer-lived sessions.

| Stage | Sectors | Wall time | RF / nRF |
|---|---:|---:|---:|
| First partition | 60 | 273.97 s | 2,580 / 0.129039 |
| Repartitioned continuation | 50 | 189.38 s | 2,354 / 0.117735 |
| **Combined** | **110** | **463.35 s** | **2,354 / 0.117735** |

This tree is identical (pairwise RF zero) to the corresponding non-persistent
60+50 result. Relative to the original 595.36-second schedule, runtime is 22.17%
lower and nRF improves by 0.001901. Restarting every five sectors was not
retained: a 15-sector test took 73.21 seconds versus an estimated 72.80 seconds
for restart-four and provides less memory headroom.

## Interpretation

This dataset provides a positive result that was not visible on the original
AliSim experiment. Boundary-aware MLIPPER D&C refinement can make substantial,
likelihood-valid topology improvements when it starts from the gap-fixed DIPPER
NJ tree. The remaining gap is primarily search coverage: the current workflow
has repaired only 14 weak components, while FastTree performs global ME search
and seven whole-tree ML-NNI rounds.

The experiment also rejects one candidate workflow: on an alignment with 82%
gaps, strict DIPPER placement remains much worse than conventional NJ even
after the gap fix. It should not be used as the starting topology here.

## Reproducibility artifacts

All generated trees, logs, checkpoint RF calculations, model files, and sweep
summaries are under:

`tmp/rnasim10k/cross-dataset-v1/`

The two MLIPPER sweep summaries are:

- `tmp/rnasim10k/cross-dataset-v1/mlipper-weak-sweep-r234/summary.json`
- `tmp/rnasim10k/cross-dataset-v1/mlipper-weak-sweep-r234-v2/summary.json`

Streaming whole-tree checkpoints and per-sector coverage manifests are under:

- `tmp/rnasim10k/cross-dataset-v1/streaming-whole-tree-nni-v3/`
- `tmp/rnasim10k/cross-dataset-v1/streaming-whole-tree-nni-v4/`
