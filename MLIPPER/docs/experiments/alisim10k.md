# ALISIM 10K / 2K Experiment Catalog

This document consolidates the major DIPPER, MLIPPER, and FastTree experiments
run on the ALISIM 10K simulation and its induced 2K benchmark. It separates
clean, blind results from exploratory or implementation-diagnostic results.

Status: historical experiment record. Paths, binaries, and commands below may
not match the current checkout. Use the root `README.md` for the supported CLI
and this document only for provenance and interpretation of recorded results.

## 1. Benchmark definitions

### Induced 2K benchmark

- Source: 10K simulated alignment and golden tree.
- Taxa used: 2,000.
- Alignment length: 7,199 sites.
- Placement construction: a 1K backbone plus 1K query tips.
- RF denominator for fully resolved 2K trees: 3,994.
- Truth: `tmp/alisim10k/pure-dipper-original-control-now.truth-induced.nwk`.
- Full MSA: `tmp/alisim10k/pure-dipper-original-control-now.full2000.fa`.

The clean DIPPER 2K baseline is a **placement tree**, not a conventional NJ
tree on all 2K taxa. It was generated with native DIPPER 0.1.5:

```bash
dipper -i m -o t -m 1 --add \
  -I full2000.fa \
  -t dipper_backbone1000.nwk \
  -O dipper-binary-add-dipper1k-2k.nwk
```

The historical 1K backbone cannot currently be reproduced exactly by running
the current DIPPER conventional-NJ mode on the stored 1K MSA: the two 1K trees
differ by RF436. Its historical input/version settings should therefore be
treated as unresolved provenance.

### Full 10K routing benchmark

- Starting backbone: 1K taxa.
- Queries: 9K taxa.
- Primary routing diagnostic: whether a method selects the golden-compatible
  attachment region for each query.
- These routing results must not be mixed with the induced 2K RF tables.

## 2. Result-status labels

| Label | Meaning |
|---|---|
| **Clean/blind** | No golden labels were used to propose, rank, or accept moves. |
| **Diagnostic** | Golden truth was used only after the tree was produced to measure RF or classify errors. |
| **Exploratory/oracle-assisted** | Golden labels affected the starting tree or move selection; RF cannot be reported as an unassisted method result. |
| **Implementation diagnostic** | Intended to test correctness, parity, or runtime rather than tree accuracy. |
| **Incomplete** | The run failed or was interrupted before a validated final tree was written. |

## 3. Authoritative clean 2K baselines

| Method | Construction | RF | nRF | Status |
|---|---|---:|---:|---|
| Native DIPPER control | Historical DIPPER 1K backbone + native placement of 1K queries | 851 | 0.213070 | Clean/blind; partly unresolved |
| Native DIPPER binary/smoothed | Binary representation of DIPPER construction | 882 | 0.220831 | Clean/blind |
| FastTree NJ only | De novo NJ on all 2K; no NNI, SPR, or ML-NNI | 810 | 0.202804 | Clean/blind |
| FastTree NJ + ME | NJ followed by ME-NNI/SPR; no ML-NNI | 874 | 0.218828 | Clean/blind |
| FastTree complete | De novo NJ + ME-NNI/SPR + ML-NNI | **724** | **0.181272** | Clean/blind; best verified 2K result |

Important interpretation:

- FastTree's ME phase is not RF-monotone: RF810 becomes RF874 before ML-NNI
  reaches RF724.
- This does not mean FastTree accepts negative moves under its active
  criterion. It means ME and RF are different criteria, and the ME phase
  changes the topology basin available to the later ML search.

Artifacts:

- `tmp/alisim10k/baseline-identify-rf.raxml.rfDistances`
- `tmp/alisim10k/fasttree-nj-only-2k.log`
- `tmp/alisim10k/fasttree-stage-rf.raxml.rfDistances`
- `tmp/alisim10k/full2000.fasttree-gtr.nwk`

## 4. Backbone substitution and native DIPPER placement

The same 1K query set was placed onto two different 1K backbones with the
native DIPPER CLI.

| 1K starting backbone | Backbone RF | Final 2K RF | Final nRF |
|---|---:|---:|---:|
| Historical DIPPER backbone | approximately 369-386, depending on representation | **851** | **0.213070** |
| FastTree de novo 1K backbone | **330** | 1287 | 0.322233 |
| FastTree 1K topology with fixed-topology Gamma4 branch optimization | 330 | 1254 | 0.313971 |

Conclusions:

- A lower induced-backbone RF does not imply better compatibility with
  DIPPER's sequential distance placement.
- The FastTree topology remains RF330 after query placement, so DIPPER does
  not destroy that backbone. The additional error is query-related.
- Relative to the historical DIPPER backbone run, the raw FastTree-backbone
  run increases query-only wrong splits from 166 to 323 and mixed-extension
  errors from 32 to 100.
- Gamma4 branch normalization recovers only 33 RF, so branch-model mismatch is
  not the primary explanation.

Artifacts:

- `tmp/alisim10k/dipper-binary-add-dipper1k-2k.nwk`
- `tmp/alisim10k/dipper-binary-add-fasttree1k-2k.nwk`
- `tmp/alisim10k/dipper-binary-add-fasttree1k-gamma4bl-2k.nwk`

## 5. Wrapper-path placement experiments

These paths are integration-wrapper experiments and are not native DIPPER CLI
results.

| Wrapper path | Original 1K input | FastTree 1K input | Notes |
|---|---:|---:|---|
| `MLIPPER_PURE_DIPPER_BASELINE=1` | RF2522 | RF1516 | D&C shadow-tree path |
| `MLIPPER_PURE_DIPPER_ITERATIVE=1` | not the native control | RF1302 | Iterative wrapper path |
| Iterative wrapper after Gamma4 backbone optimization | - | RF1276 | Different path from native `--add` |

These numbers explain wrapper behavior only. They must not replace RF851,
RF1287, or RF1254 in the native-DIPPER comparison.

## 6. FastTree refinement from placement trees

FastTree requires a binary `-intree`. Resolving DIPPER multifurcations adds
arbitrary zero-information splits, so both the original and binary-start RF
are reported.

| Starting construction | Original RF | Binary FastTree input RF | Refined RF | nRF | FastTree de novo reference |
|---|---:|---:|---:|---:|---:|
| Historical DIPPER backbone + DIPPER queries | 851 | 880 | **762** | 0.190786 | 724 |
| FastTree 1K backbone + DIPPER queries | 1287 | 1310 | **780** | 0.195293 | 724 |

The refinement used FastTree's default full-tree trajectory:

- up to 44 rounds of ME-NNI;
- two ME-SPR rounds with maximum SPR length 10;
- up to 22 ML-NNI rounds, converging after three on these runs;
- GTR+CAT20 branch/model optimization.

Clean DIPPER binary-start split transition, RF880 to RF762:

| Transition | Count |
|---|---:|
| Original wrong splits removed | 358 / 440 |
| Original wrong splits retained | 82 / 440 |
| Missing truth splits recovered | 148 / 440 |
| Missing truth splits still missing | 292 / 440 |
| Originally correct splits lost | 91 |
| New wrong splits introduced | 297 |

FastTree can rewrite most DIPPER errors, but it converges to a different basin
from de novo FastTree. The residual RF is not simply a set of immutable DIPPER
splits.

Artifacts:

- `tmp/alisim10k/fasttree-refine-dipper1k-dipper2k.default.nwk`
- `tmp/alisim10k/fasttree-refine-ft1k-dipper2k.default.nwk`
- `tmp/alisim10k/dipper-binary-vs-fasttree-with-refinement-split-audit.tsv`

## 7. Split-by-split DIPPER versus FastTree audit

| Split class | Count |
|---|---:|
| Correct in both | 1475 |
| Correct only in DIPPER | 71 |
| Correct only in FastTree | 158 |
| Missing from both | 293 |
| Wrong in both | 82 |
| Wrong only in DIPPER | 318 |
| Wrong only in FastTree | 278 |

FastTree's 127-RF advantage consists of 87 net additional correct splits and
40 fewer wrong splits.

The 158 truth splits found by FastTree but missed by DIPPER are distributed
across the tree:

| Origin | Count |
|---|---:|
| Query-only | 51 |
| Backbone-only side | 61 |
| Backbone-pendant mixed | 20 |
| Mixed; backbone projection already present | 13 |
| Backbone projection itself missing | 13 |

The difference is therefore not attributable to one isolated backbone or
placement failure mode.

Artifacts:

- `tmp/alisim10k/dipper-vs-fasttree-split-audit.tsv`
- `tmp/alisim10k/dipper-vs-fasttree-split-audit.summary.txt`

## 8. Post-placement query-only and local SPR pilots

This earlier experiment family started at RF902 (`nRF=0.225839`) and tested
limited query-only/local repair.

| Experiment | Final RF | nRF | Change from RF902 |
|---|---:|---:|---:|
| Single target T102 | 902 | 0.225839 | 0 |
| Single target T1155 | 900 | 0.225338 | -2 |
| Single target T1352 | 904 | 0.226340 | +2 |
| All clusters, unrestricted acceptance | 915 | 0.229094 | +13 |
| All clusters, gain threshold 5 | 902 | 0.225839 | 0 |
| All clusters, gain threshold 10 | 898 | 0.224837 | -4 |
| All-cluster SPR, gain threshold 10 | **896** | **0.224337** | **-6** |
| Repeated gain-10 rounds | 898 | 0.224837 | -4; no further change |

Conclusion: exact boundary context makes local moves valid, but query-only
repair produces only small RF changes and can regress when applied broadly.

Artifacts: `tmp/alisim10k/postplacement*-truth-rf.raxml.rfDistances`.

## 9. Backbone-only SPR pilots

A 1K backbone was refined with prune size limited to at most four tips.

| Pass | Backbone RF | nRF |
|---|---:|---:|
| Radius 1 | 381 | 0.191073 |
| Radius 2 | 373 | 0.187061 |
| Radius 3 | **365** | **0.183049** |
| Radius 4 | 365 | 0.183049 |

The improvement saturates at radius 3. Replacing or modifying the backbone
before placement did not translate into a comparable final-tree improvement;
one iterative pre-SPR construction ended at RF886.

Artifacts:

- `tmp/alisim10k/backbone-prespr-max4-r*-truth-rf.raxml.rfDistances`
- `tmp/alisim10k/pure-dipper-iterative-prespr1000-truth-rf.raxml.rfDistances`

## 10. MLIPPER whole-tree short-edge NNI/SPR sequence

### Provenance warning

The RF806 starting tree is **not clean**. It descends from an RF818 tree whose
four quartet moves were selected from rows labeled `truth=improving`, followed
by two MLIPPER NNI moves. The following series demonstrates search behavior
and implementation capability, not blind end-to-end accuracy.

The whole-tree mode did not use individual DIPPER placement anchors. It swept
the complete 2K tree, restricted central edges to branch length at most
`1e-6`, retained top-K 32 candidates, and performed exact likelihood
validation.

| Step | Move | Start RF | Final RF | nRF | Accepted behavior |
|---:|---|---:|---:|---:|---|
| 1 | Exact NNI, radius 1 | 810 | 806 | 0.201803 | One exploratory improvement stage |
| 2 | SPR radius 2 | 806 | 800 | 0.200300 | Two candidates accepted |
| 3 | Exact NNI | 800 | 798 | 0.199800 | One candidate accepted |
| 4 | SPR radius 2 | 798 | 796 | 0.199299 | One candidate accepted |
| 5 | Exact NNI | 796 | 796 | 0.199299 | Converged for that NNI pass |
| 6 | SPR radius 3 | 796 | 790 | 0.197797 | Two candidates accepted |
| 7 | Exact NNI | 790 | 790 | 0.197797 | RF-neutral topology change |
| 8 | SPR radius 3 | 790 | **788** | **0.197296** | One candidate accepted |

Later bridge, tabu, larger top-K, and repeated radius-2/3 attempts did not
improve beyond RF788 in this exploratory lineage.

Artifacts:

- `tmp/alisim10k/internal1e6-step*.log`
- `tmp/alisim10k/internal1e6-step*-rf.raxml.rfDistances`
- `tmp/alisim10k/rf788-shortedge1e6-*.log`

## 10.1 Guided overlapping-sector trajectory from the clean DIPPER tree

This clean experiment used FastTree only as a topology guide. It did not use
the truth tree for proposal generation or likelihood acceptance.

1. Boundary-consensus sectors failed: the only clearly positive-likelihood
   proposal changed RF880 to RF882.
2. Actual-taxon sectors with a fixed constrained boundary succeeded. A blind
   small-sector sweep changed RF880 to RF874.
3. Rebuilding the same sector under 4, 8, 12, and 16 outgroup contexts and
   choosing the pairwise-RF medoid produced stable moves. Three accepted
   sectors changed RF874 to RF866, RF840, and RF832.
4. Replacing bounded sectors whose boundary split was shared with the
   FastTree guide reduced the best RF to 796. The trajectory exhausted all
   shared-boundary sectors of at most 800 tips; the next available shared
   boundary contained 1,688 tips.
5. Gamma4 smoothing followed by whole-tree short-edge MLIPPER refinement
   reached **RF792, nRF=0.198297**. Exact NNI and SPR radii 1, 2, and 5 then
   converged; radius 4 supplied the final two-RF improvement.
6. The remaining guide disagreement had no shared boundary between 312 and
   1,688 tips. An adaptive 1,688-tip escalation replaced that bounded sector
   with the guide-induced topology. After identical Gamma4 smoothing, logL
   improved from -128297.550470 to -127773.616093 and the result reached
   **RF722, nRF=0.180771**.

The clean FastTree-style reference is RF762, nRF=0.190786, so the adaptive
trajectory exceeds the requested accuracy target by 40 RF. This is an upper
scale control rather than the final scalable design: the last sector contains
1,688 of 2,000 tips. It proves that the residual is cross-sector topology and
not an immutable MLIPPER likelihood barrier, while leaving sector-size
reduction as the next implementation problem.

Artifacts:

- `tmp/alisim10k/actual-sector-blind-sweep/`
- `tmp/alisim10k/sector-stability-T448-b160/`
- `tmp/alisim10k/sector-stability-T1831-b256-after-T448/`
- `tmp/alisim10k/sector-stability-T1257-b160-after-rf840/`
- `tmp/alisim10k/guided-shared-sector-trajectory/`
- `tmp/alisim10k/guided-rf796-mlipper-refine/`
- `tmp/alisim10k/guided-large-sector-1688/`

## 11. Weak-edge correlation

The short-edge association was reproduced on clean trees and is not dependent
on the oracle-assisted RF806 lineage.

| Tree | Error rate for `b <= 1e-6` | Error rate for longer edges | Wrong splits that are short |
|---|---:|---:|---:|
| Clean DIPPER RF851 | 78.1% | 9.3% | 62.3% |
| Clean binary/smoothed RF882 | 72.7% | 3.2% | 89.3% |

On the exploratory RF806 tree, short length plus low SH-like support selected
445 edges in 314 connected components; 324 of the 445 selected edges were
wrong. Weak status is therefore a strong error-localization signal.

It is not a direction signal:

- 46 of the 71 correct splits unique to DIPPER are also at most `1e-6`;
- changing every short edge would destroy correct weak resolutions as well as
  repair incorrect ones.

## 12. Weak-component exhaustive topology inference

Connected weak regions were represented by boundary ports while every outside
subtree contributed an exact directional boundary CLV.

| Ports | Unrooted binary candidates |
|---:|---:|
| 4 | 3 |
| 5 | 15 |
| 6 | 105 |

Selected Gamma4 exhaustive results on the exploratory RF806 tree:

| Component | Ports | Gamma4 winner RF | Best candidate RF | Winner margin |
|---|---:|---:|---:|---:|
| rank85 | 4 | 806 | 804 | `7.04e-7` |
| rank257 | 4 | 804 | 804 | `4.56e-7` |
| rank310 | 4 | 806 | 804 | `1.33e-6` |
| rank28 | 5 | 804 | 802 | `6.12e-8` |
| rank83 | 5 | 806 | 804 | `3.64e-7` |
| rank12 | 6 | 802 | 800 | `2.43e-7` |
| rank15 | 6 | 804 | 800 | `2.91e-10` |

The reference-compatible topology is usually present in the candidate set,
but the full-data ML winner is not reliably the reference resolution.

## 13. RELL and SH-like stability tests

| Component | Candidates | Winner RELL support | Runner-up | Decision |
|---|---:|---:|---:|---|
| rank257 | 3 | 0.5919 | 0.4053 | Defer |
| rank15 | 105 | at most 0.1500 | 0.0849 for full-data runner-up | Defer |

The likelihood margins per site are at numerical-tie scale. Requiring high
RELL support turns the method into a reliable ambiguity detector but produces
almost no topology correction.

Across all 230 four-port weak components:

- nine full-data winners had RELL frequency at least 0.9;
- seven also excluded every alternative at the 5% level;
- only two proposed topology changes;
- both changes were RF-neutral;
- the strict workflow remained at RF806.

## 14. ML-constrained ME/NJ bootstrap

The tested rule was: ML rejects incompatible candidates; ME/NJ bootstrap
provides direction only within the ML-equivalent set.

Rank257 pilot, support for the RF806 start / RF804 ML winner / third topology:

| Port representation | Distance | Supports | Decision |
|---|---|---|---|
| One nearest representative | p-distance | 0.7923 / 0.0262 / 0.1815 | Defer |
| One nearest representative | JC | 0.6947 / 0.0300 / 0.2754 | Defer |
| Up to eight representatives | p-distance | 0.8141 / 0.0001 / 0.1857 | Defer |
| Up to eight representatives | JC | 0.7444 / 0.0002 / 0.2555 | Defer |

Batch result over 230 four-port components:

- agreement and support filtering selected 27 components;
- 16 proposed topology changes;
- unconstrained ME changes would predict RF816 from RF806;
- 14 of 16 changes passed the ML non-rejection rule;
- their net RF effect was zero, predicting RF806;
- choosing the full-data Gamma4 winner instead would worsen RF by four.

Conclusion: ME and ML often disagree inside statistically ambiguous regions,
and neither reliably selects the reference binary resolution.

## 15. Placement reranking and dual-region experiments

### Full 10K routing result

| Method | Golden-compatible attachment count | Rate |
|---|---:|---:|
| Original DIPPER | 7,753 / 9,000 | 86.14% |
| MLIPPER reassign winner | 7,390 / 9,000 | 82.11% |

MLIPPER reassign reduced correct attachment selection by 363 queries, or 4.03
percentage points. Successful execution of the reranker did not imply better
routing.

### Tested reranking variants

- current-backbone and RAxML-estimated model parameters;
- larger top-K candidate sets;
- MASH-supported candidates;
- likelihood-margin-gated reassignment;
- DIPPER winner plus MLIPPER winner as a dual-region candidate set;
- exact edge likelihood and branch-length diagnostics for regression queries.

General outcome:

- model-parameter replacement did not reverse the routing regression;
- forcing MLIPPER winners was unsafe;
- margin gating reduced the number of regressions but provided little total RF
  improvement;
- dual-region comparison is safer than replacing the DIPPER winner, but the
  tested local scorer still lacked a reliable direction in weak regions.

Legacy full-10K RF files under `tmp/alisim10k/reassign*-rfDistances` were
produced by multiple incompatible or partially corrected paths. They should
not be combined into a single accuracy table without reconstructing each
command and taxon set.

## 16. Branch-length and model experiments

| Experiment | Outcome |
|---|---|
| RAxML-estimated evolutionary parameters for MLIPPER reassignment | Did not make MLIPPER routing exceed DIPPER. |
| Fixed-topology Gamma4 optimization of FastTree 1K backbone before DIPPER | RF1287 to RF1254; insufficient. |
| Global branch optimization before local NNI/SPR | Improved numerical consistency but did not supply topology direction. |
| Lowering branch-length floor / double precision | Necessary for very short branches, but not sufficient to identify the correct binary split. |
| CAT versus Gamma4 comparisons | Model mismatch matters for score parity, but does not explain the main placement/topology gap. |

Branch-length optimization must always be followed by a fresh tree pack and
likelihood rescore before accepting a topology.

## 17. Correctness bugs found during experimentation

### Fresh repack/rescore after smoothing

Earlier code did not always rebuild, repack CLVs, and rescore both the smoothed
baseline and candidate before acceptance. Some apparent gains disappeared
after this was fixed. All accepted topology transactions must be validated on
the exact rebuilt state.

### Approximate-positive filtering

A weak-component candidate with approximate gain near `-0.98` had exact gain
near `+0.01`. Approximate-negative candidates cannot all be discarded in weak
regions; proposal ranking and exact acceptance must remain separate.

### Radius-5 zero-length split bug

The first clean radius-5 run accepted three likelihood-positive moves with a
rebuilt likelihood gain of `+42.8889`, then failed to serialize because
splitting a zero-length regraft edge produced a distal length of `-1e-6`.
The topology path manually computed `distal = total - proximal` instead of
using the existing zero-safe `normalize_split_branch_lengths()` helper. The
path was corrected before the validated 10K rerun described below.

### Post-placement sector round override

The post-placement D&C driver silently forced `local_spr_rounds` to one unless
whole-tree mode was active. `MLIPPER_POSTPLACEMENT_MULTIROUND_SECTOR_SPR=1`
now permits the requested number of rounds while retaining the legacy
single-round default. A radius-5 sector rerun reached round 2/3 and stopped
when the second round accepted no candidate, confirming both the override and
the convergence stop.

## 18. Bug-fixed 10K sector refinement

The bug-fixed DIPPER start has RF5219 (`nRF=0.261028`) in its original
polytomous representation. MLIPPER serializes a binary tree, so optimizer
deltas are measured against the corresponding binary baseline, RF5274
(`nRF=0.263779`); comparing a resolved output directly with RF5219 mixes
topology changes with arbitrary resolutions of the input polytomies.

| Stage | RF | nRF | Exact result |
|---|---:|---:|---|
| Binary bug-fixed DIPPER baseline | 5274 | 0.263779 | Starting point |
| Top-100 short-edge radius-2 SPR | 5270 | 0.263579 | 11 positive transactions, `+2.6771` logL |
| Broad exact NNI on the same sectors | 5270 | 0.263579 | 4 positive transactions, `+1.1585` logL; RF-neutral |
| Radius-5 sector SPR, top 100 | **5252** | **0.262679** | 26 positive transactions, `+15.7076` logL |

Radius 5 therefore repairs 22 RF relative to the binary baseline, but remains
786 RF worse than the established DIPPER+FastTree result (RF4466). Positive
Gamma4 gain is not an RF-direction oracle: several accepted sector
transactions increased RF, and the largest gain (`+8.8274`) improved RF by
only four. The full-sector rebuilt-likelihood rollback remains mandatory;
74/100 radius-5 sectors failed that transaction gate.

The supported optimizer configuration is consequently adaptive rather than a
uniform radius-5 sweep: use radius 2 as the cheap proposal stage, escalate
selected unresolved sectors to radius 5, retain streaming boundary CLVs and
symmetric branch smoothing, and commit only a non-negative rebuilt sector
transaction. Exact NNI is not a useful standalone stage for this 10K start.
Multi-round sector SPR is available for revisiting profitable sectors, but is
not enabled globally because the measured revisit added only about `7e-6`
logL before the next round converged.

### 10K actual-taxon guided-sector replication

The 2K actual-taxon sector workflow was scaled to all 10,000 taxa. Truth was
used only after each output was written to compute RF. The first trajectory
used the binary de novo FastTree topology as a guide; a second trajectory used
the binary DIPPER+FastTree-refined topology as guidance. Every replacement
preserved a guide-shared boundary and was limited to at most 1,750 real tips.

| Stage | RF | nRF | Fixed-model Gamma4 logL |
|---|---:|---:|---:|
| Binary bug-fixed DIPPER baseline | 5274 | 0.263779 | -370378.003523 |
| De novo-guide bounded sectors, round 13 | 4594 | 0.229769 | -370330.227313 |
| Refined-guide bounded sectors, round 5 | **4486** | **0.224367** | **-370330.227575** |
| DIPPER+FastTree-style reference | 4476 | 0.223867 | not used for acceptance |

The final guided-sector tree is within 10 RF (`0.000500` nRF) of the
FastTree-style reference and improves the baseline likelihood by 47.775948.
Its likelihood differs from the de novo-guide checkpoint by only 0.000262
(`1.6e-8` per site), so the second trajectory remains in the Gamma4
equivalence set. This reproduces the 2K accuracy trend without the 84%-of-tree
adaptive escalation: the largest 10K sector contained 1,735 taxa.

This result is a FastTree-guided replication, not an independent replacement
for FastTree. It demonstrates that bounded D&C sector transactions can carry
the useful FastTree topology trajectory while retaining the DIPPER/MLIPPER
tree outside each active sector.

Artifacts:

- `tmp/alisim10k/guided10k-binary-trajectory/`
- `tmp/alisim10k/guided10k-binary-eval/`
- `tmp/alisim10k/guided10k-refined-guide-trajectory/`
- `tmp/alisim10k/guided10k-refined-guide-eval/`

### 10K MLIPPER-only blind D&C sector search

To separate MLIPPER refinement from FastTree guidance, a new topology-only
D&C route selects components solely from short edges in the current DIPPER
tree. MLIPPER enumerates radius-2 SPR candidates inside a 300-tip sector and
validates them with streaming boundary messages and rebuilt Gamma4
likelihood. Radius 3 and 4 are attempted only when the smaller radius accepts
no move. FastTree and the truth tree are not loaded until the sweep ends.

The first five blind attempts accepted three radius-2 moves. They reduced RF
from 5274 to 5264 (`nRF 0.263779` to `0.263279`). After the same fixed-model
RAxML branch optimization, log likelihood changed from `-370378.003523` to
`-370378.003517`, an effectively tied gain of about `0.000006`. This is a
genuine MLIPPER-only improvement, but it is still far from the guided result
and does not replace FastTree.

This run also exposed and corrected a symmetric-smoothing acceptance bug. A
candidate could beat a separately smoothed baseline copy while still scoring
below the currently accepted tree. Acceptance now requires it to beat both,
and a rejected baseline copy no longer overwrites the current per-site
likelihood state.

A follow-up tested sector-local tip placement before SPR. In the second-ranked
weak component, unrestricted regrafting of its three terminal tips accepted
two moves (`+0.020558` sector log likelihood), but they were RF-neutral.
Subsequent radius-2 SPR reached RF5272 (`nRF 0.263679`). A matched SPR-only
control reached the same RF with the same accepted SPR move, so single-tip
placement did not change the search basin. Coordinated subtree/group placement
is required for the next sectorial-search prototype.

Artifacts:

- `tmp/alisim10k/mlipper-only-blind-sweep10k/summary.json`
- `tmp/alisim10k/mlipper-only-blind-sweep10k/final.nwk`
- `tmp/alisim10k/mlipper-only-blind-sweep10k/final-rf.raxml.rfDistances`
- `tmp/alisim10k/mlipper-only-blind-sweep10k/final-eval.raxml.bestTree`

### Edge-owned whole-tree ML-NNI replication

The faster RNASim workflow was replicated without FastTree or truth guidance:
connected 300-tip edge-owned sectors, six farthest-point anchors, top K 32,
two-GPU exact validation, persistent-session restart every four sectors, and a
60-sector sweep followed by a 50-sector repartitioned continuation. The old
RF5274 bug-fixed start had been removed during temporary-file cleanup, so the
retained blind-refinement tree (RF5264) was branch-optimized with a newly fitted
GTR+Gamma4 model and used as the explicit starting point. Its prior 10-RF gain
is not counted as part of this experiment.

| Stage | Sectors | Wall time | RF | nRF | Accepted NNI | Audited logL |
|---|---:|---:|---:|---:|---:|---:|
| Re-fitted start | - | 115.65 s preparation | 5,264 | 0.263279 | - | -370378.009333 |
| First partition | 60 | 381.52 s | 5,258 | 0.262979 | 7 | -370360.322751 |
| Repartitioned continuation | 50 | 318.28 s | 5,258 | 0.262979 | 2 | -370355.307264 |
| **Refinement total** | **110** | **699.80 s** | **5,258** | **0.262979** | **9** | **-370355.307264** |

The topology refinement took 11m39.8s; including model preparation took
13m35.4s. Likelihood improved by 22.702069 under the re-fitted model, but RF
improved by only six. The second sweep's meaningful gain (`+5.015484`) was
RF-neutral, and its other accepted move gained only `2.2e-6`. In contrast to
RNASim, broad exact-NNI coverage therefore does not escape the AliSim DIPPER
search basin. Additional NNI sectors are not justified; a placement or
larger-subtree SPR proposal mechanism is required here.

### Batched weak-sector tip reinsertion

A 2K pilot removed 72 tips from a blindly selected short-edge sector, then
reinserted them in batches of five. Placement and topology refinement used
separate GPU sessions; the SPR phase retained a 300-tip D&C sector with
directional boundary messages and allowed only the latest query tips to move.

| Binary-tree method | RF | nRF |
|---|---:|---:|
| Original DIPPER tree | 851 | 0.213070 |
| Five-tip reinsertion only | 849 | 0.212569 |
| Reinsertion + query-only adaptive SPR, uncontrolled re-binarization | 870 | 0.217827 |
| Reinsertion + query-only adaptive SPR, zero-only collapse | **847** | **0.212068** |

The original SPR run was confounded by D&C output restoring zero-length parser
resolutions. Setting both writers to collapse only exactly zero-length output
edges makes every zero-move batch split-identical before and after SPR. Under
this corrected rule, five accepted query moves improve RF849 to RF847: six
splits are replaced, with one correct split removed and two correct splits
added. The gain is real but small.

An apparent reinsertion-only RF682 result was invalid for binary comparison:
the legacy small-tip writer collapsed 350 internal edges below `1e-6`, whereas
the D&C writer emitted all binary resolutions. The main-like path now supports
`MLIPPER_MAINLIKE_COLLAPSE_EPSILON=0`, which was used for the table above.

Artifacts:

- `tmp/alisim10k/smalltip-sector-reinsert2k/reinsert-mainlike-batch5-binary-only.nwk`
- `tmp/alisim10k/smalltip-sector-reinsert2k/reinsert-query-only-spr-binary-v5.nwk`
- `tmp/alisim10k/smalltip-sector-reinsert2k/batched-query-only-spr-binary-v5/`

## 19. Runtime and GPU optimization experiments

The detailed gene638 branch-length optimization profile is maintained in
[`../performance/full_tree_blo_optimization_profile.md`](../performance/full_tree_blo_optimization_profile.md).
Key results:

| Implementation | One-sweep BLO time | Full wall time |
|---|---:|---:|
| Original single-block Newton | 1.5161 s | 1.94 s |
| Eight-block Newton | 0.8041 s | 1.21 s |
| Incremental multi-block | 0.6771 s | 1.11 s |
| Warp-site incremental | 0.5513 s | 1.00 s |
| Warp-site plus branch-local CLVs | **0.3339 s** | **0.77 s** |

The final version reduced one-sweep BLO time by about 4.5x. GPU launch
starvation was not the main bottleneck; the measured kernel window was 95.4%
busy.

## 20. Current conclusions

1. Short internal branches are a strong error-localization signal, not a
   reliable topology-direction signal.
2. Exact boundary CLVs solve the missing-outside-context problem without
   storing the entire tree's CLVs.
3. Exhaustive local enumeration solves candidate coverage for 4-6 ports, but
   ML, RELL, and ME bootstrap often cannot identify one stable binary answer.
4. Query-only and small-radius local repair improve RF only modestly.
5. FastTree's de novo construction trajectory remains the best verified 2K
   result at RF724.
6. Full FastTree refinement can reduce a clean binary DIPPER start from RF880
   to RF762, showing that most DIPPER errors are movable, but the starting
   topology still leads to a different final basin.
7. A larger SPR radius is worth testing, but radius alone does not address
   top-K proposal pruning, multi-move barriers, or ambiguous acceptance.
8. Sectorial search is the most relevant next architectural experiment:
   rebuild a weak connected sector jointly, retain exact boundary CLVs, use a
   cheap criterion for proposals, and validate the complete replacement under
   rebuilt Gamma4 likelihood.

## 21. Recommended reporting table

Use this short table in presentations unless a specific experimental family is
being discussed:

| Method | RF | nRF | Evidence class |
|---|---:|---:|---|
| FastTree de novo 2K | **724** | **0.181272** | Clean/blind |
| FastTree full refinement from binary clean DIPPER | 762 | 0.190786 | Clean/blind |
| FastTree full refinement from FastTree-1K + DIPPER | 780 | 0.195293 | Clean/blind |
| MLIPPER short-edge NNI/SPR best | 788 | 0.197296 | Exploratory/oracle-assisted start |
| Native DIPPER 1K + 1K placement | 851 | 0.213070 | Clean/blind |
| Native FastTree-1K + DIPPER placement | 1287 | 0.322233 | Clean/blind |

Do not present RF788 as the best clean MLIPPER result, and do not compare the
wrapper RF2522/RF1516 values directly with native DIPPER RF851/RF1287.
