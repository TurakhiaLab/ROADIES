# D&C SPR Build Plan

## Goal

Add a whole-tree, divide-and-conquer SPR refinement stage after D&C NNI. The
new stage must search beyond one-NNI topology changes, reuse directional sector
likelihood boundaries, and support experiments testing whether regions changed
by NNI are enriched for errors that require SPR to repair.

This is distinct from the existing small-tip `runLocalSPR()` workflow. That
workflow builds repair anchors from recently inserted query tips; D&C SPR must
derive its prune ownership and search regions from the current whole tree.

## Proposed workflow

```text
D&C NNI to its configured stopping condition
  -> record every accepted NNI region
  -> rebuild the current whole-tree split index
  -> partition owned prune edges into connected sectors
  -> search non-NNI-equivalent SPR moves at increasing radii
  -> replay accepted local moves on the authoritative full tree
  -> rebuild stale split/sector state after each accepted move
  -> run final model and branch-length optimization
```

The correctness-first implementation accepts at most one revalidated SPR move
before rebuilding the whole-tree split index and sector partition. Later work
may accept multiple edge-disjoint moves from one sector after equivalence is
demonstrated.

## Radius semantics

For the usual SPR distance measured between the original attachment region and
the regraft edge on the edge graph:

```text
radius 0: reconnect to the original edge; no topology change
radius 1: move to an adjacent edge; normally equivalent to one NNI
radius 2+: may contain changes that cannot be made by one NNI
```

The implementation's current radius counting may include the prune endpoint,
suppressed parent, or another path endpoint. Therefore, a numeric radius must
not be assumed to mean a mathematical SPR distance until an enumeration test
establishes the mapping.

The authoritative definition of an SPR-only candidate is:

> The candidate topology cannot be obtained from the current topology by one
> NNI move.

Candidate filtering should compare canonical split sets, not rely only on the
configured radius. The initial experimental schedule should therefore begin at
the smallest empirically verified non-NNI radius, followed by a schedule such
as `{2, 4, 8, 16}` after accounting for any radius offset.

## Public API

Keep insertion repair and whole-tree refinement separate:

```cpp
struct DivideAndConquerSPROptions {
    int core_edges = 80;
    int initial_radius = 2;
    int maximum_radius = 8;
    int radius_step = 1;
    int max_sweeps_per_radius = 1;
    int topk_per_sector = 128;
    int branch_smoothing_sweeps = 1;
    bool exclude_nni_equivalent_moves = true;
};

struct DivideAndConquerSPRResult {
    int radius_levels = 0;
    int sweeps = 0;
    size_t sectors = 0;
    int accepted_moves = 0;
    bool converged = false;
};

DivideAndConquerSPRResult runDivideAndConquerSPR(
    const DivideAndConquerSPROptions& options);
```

`runLocalSPR()` remains the recently inserted-tip repair path.
`runDivideAndConquerSPR()` becomes the whole-tree sector path.

## Reusing the D&C NNI architecture

Reuse the following mechanisms from `runDivideAndConquerNNI()`:

- canonical split fingerprints for identities that survive node renumbering;
- connected whole-tree edge partitioning;
- bounded subtree extraction with directional boundary messages;
- local-to-global move replay;
- sector cache invalidation;
- OOM recovery by splitting only the failing sector and retrying it.

Generalize subtree preparation from an NNI-specific operation to a topology
sector operation:

```cpp
prepare_topology_sector_with_directional_boundaries(
    full_tree,
    owned_prune_edges,
    halo_radius,
    ...);
```

The extracted sector has two regions:

- **Owned core:** edges allowed to initiate a prune operation.
- **Halo:** additional topology needed for regraft search, likelihood context,
  and branch smoothing.

Directional virtual tips represent omitted full-tree components. They provide
likelihood input but must never be selected as a prune or regraft endpoint.

The extraction halo must cover the SPR radius plus sufficient neighborhood for
edge suppression, regraft context, and branch smoothing.

## Decoupling move type from execution domain

The current refinement code infers directional-subtree execution from
`move_type == NNI`. That does not support SPR inside a D&C sector. Introduce an
explicit domain:

```cpp
enum class RefinementDomain {
    ResidentFullTree,
    DirectionalSector,
};
```

The refinement engine should then receive both independent choices:

```text
move type: SPR or NNI
domain: resident full tree or directional sector
```

SPR sector execution also needs an explicit `owned_spr_prune_root_ids` set.
NNI central-edge ownership and SPR prune-edge ownership should remain separate
because they have different meanings.

## Sector ownership

Assign each legal SPR candidate to exactly one sector using its prune edge:

```text
prune edge: must belong to the sector's owned core
regraft edge: may belong to the owned core or halo
directional boundary node: never a legal endpoint
```

This prevents duplicate work while allowing moves across a sector's core
boundary. Candidate records should carry canonical full-tree identities:

```cpp
struct OwnedSPRCandidate {
    SplitFingerprint prune_split;
    SplitFingerprint regraft_split;
    int radius = 0;
    double approximate_gain = 0.0;
};
```

Node IDs are only temporary resolved coordinates. Before exact validation, map
the fingerprints against the current tree. If either split no longer exists,
discard the stale candidate and regenerate its sector.

## Candidate scoring and acceptance

Use two stages:

```text
enumerate legal owned-prune/regraft pairs
  -> approximate GPU placement-style scoring
  -> retain sector top-K positive candidates
  -> materialize candidate topology
  -> symmetrically smooth baseline and candidate branches
  -> exact likelihood comparison
  -> full-tree replay and audit
```

Approximate placement scoring optimizes the attachment coordinates: pendant
length and the proximal/distal split of the target edge. Exact validation must
smooth a broader, symmetric edge set for both the baseline and candidate:

- the prune edge and edges joined after suppression;
- both pieces of the split regraft edge;
- the new pendant edge;
- the prune-to-regraft path;
- a one-edge neighborhood around those structures.

Accept a move only if it improves both the symmetrically optimized baseline and
the currently accepted tree likelihood within the configured numerical
tolerance. After replay, rebuild and audit the full-tree likelihood before the
move becomes durable.

## Accepted NNI trace and hotspot hypothesis

The experimental hypothesis is:

> Regions in which NNI accepted a move are enriched for remaining topology
> errors that require a non-NNI SPR move.

Record each accepted NNI using topology-stable information:

```cpp
struct AcceptedNNITrace {
    int sweep = 0;
    int sector = 0;
    SplitFingerprint old_central_split;
    SplitFingerprint new_central_split;
    std::array<SplitFingerprint, 4> quartet_clades;
    double approximate_gain = 0.0;
    double exact_gain = 0.0;
};
```

The old split disappears after the move, so the four A/B/C/D clade
fingerprints define the persistent changed region. For every potential SPR
prune edge, calculate its edge-graph distance to the nearest accepted NNI
region.

NNI-hot regions should be a search priority, not the only search space. Some SPR
improvements require crossing an intermediate likelihood valley and may occur
far from every accepted NNI.

## Radius and stopping policy

Support two policies:

1. **Single sweep per radius:** predictable cost; advance after one whole-tree
   sweep.
2. **Converge per radius:** repeat a radius until no move is accepted, then
   increase it.

No accepted move at one radius does not imply that larger radii are exhausted.
Optional early stopping may require several consecutive empty radii. Any
accepted move invalidates the current partition; the first implementation must
repartition before continuing.

## Implementation phases

### Phase 0: radius oracle

- Enumerate candidates on small known trees.
- Canonicalize every resulting topology by its split set.
- Determine which configured radius first produces nontrivial moves.
- Mark candidates obtainable by exactly one NNI.
- Establish the numeric-radius to mathematical-distance mapping.

### Phase 1: correctness-first D&C SPR

- Add the public options, result, and session entry point.
- Support one fixed, verified non-NNI radius.
- Partition whole-tree prune edges using the NNI sector machinery.
- Extract directional sectors with owned core and halo.
- Accept at most one globally revalidated move before repartitioning.
- Reject virtual-boundary endpoints and stale fingerprints.

### Phase 2: radius schedule and instrumentation

- Add increasing radius schedules and both stopping policies.
- Record candidate, acceptance, likelihood, timing, and memory metrics.
- Record distance from every SPR candidate to an accepted NNI region.

### Phase 3: performance

- Reuse persistent sector workspaces across compatible sectors.
- Split and retry a sector after CUDA allocation failure.
- Accept multiple provably edge-disjoint moves before rebuilding.
- Rebuild only affected partitions after correctness equivalence is shown.

## Correctness tests

For trees of roughly 8--20 tips, compare D&C SPR with an exhaustive full-tree
oracle. Verify that:

- every legal owned prune/regraft pair is enumerated exactly once;
- sector-boundary moves are not omitted;
- virtual boundary nodes are never endpoints;
- NNI-equivalent candidates are identified independently of radius numbering;
- local replay and direct full-tree application produce identical split sets;
- local directional and resident full-tree likelihoods agree within tolerance;
- accepting a move never decreases audited full-tree likelihood;
- accepting one move invalidates and regenerates stale candidates;
- OOM sector splitting produces the same result as the unsplit search.

Also retain an NNI-equivalent SPR control. After D&C NNI convergence, this
control should accept no meaningful improvement. A positive result may expose
NNI coverage gaps, premature stopping, asymmetric branch optimization, or a
directional-boundary discrepancy.

## Main experiment matrix

Use identical starting trees, alignments, models, NNI settings, final optimizer,
precision, and hardware.

| Group | Refinement after D&C NNI |
| --- | --- |
| A | None: NNI baseline |
| B | NNI-equivalent SPR control only |
| C | Short non-NNI SPR |
| D | Increasing non-NNI SPR radii |
| E | NNI-hot regions only |
| F | Matched NNI-cold regions with the same search budget |
| G | Random regions with the same search budget |
| H | NNI-hot first, then global D&C SPR |
| I | Full global D&C SPR |

Run runtime measurements at least three times per dataset. Topology outputs and
likelihood audits should be deterministic or report observed variation.

## Metrics

Record per radius, sweep, and sector:

```text
radius and verified move class
sector count and owned prune-edge count
enumerated candidates
NNI-equivalent candidates filtered
approximately positive candidates
exactly validated candidates
accepted moves
distance to nearest accepted NNI region
log likelihood before and after
sector extraction, scoring, replay, and rebuild time
peak GPU memory
OOM split/retry count
```

With a reference tree, additionally record false-positive and missing splits,
RF distance, whether each SPR move fixes a reference split, and the minimum
correcting SPR radius. Without a reference tree, describe the outcome as an
SPR-improvable region rather than a proven topology error.

## Hotspot analysis

Primary descriptive quantities are:

```text
P(SPR-improvable or reference-error edge | distance to accepted NNI = d)
recall of all accepted/correcting SPR moves captured within distance <= k
likelihood or RF improvement per fixed search budget
```

Compare NNI-hot edges with matched cold controls rather than only random edges.
Match or adjust for branch length, support, subtree size, alignment entropy,
missingness, and dataset/tree identity. Report odds ratios with bootstrap
confidence intervals; a mixed-effects model can account for correlated edges
within one tree.

When practical, use training sites to generate NNI/SPR candidates and held-out
sites to validate gains. This separates true predictive enrichment from both
optimizers responding to the same likelihood noise.

## Success criteria

The implementation is ready for performance work only after:

- the radius oracle distinguishes one-NNI and non-NNI SPR moves;
- exhaustive and D&C candidate coverage agree on small trees;
- local-to-full replay preserves the expected canonical split set;
- every accepted move passes a non-decreasing full-tree likelihood audit;
- NNI-equivalent control behavior is understood;
- hotspot, matched-control, and global searches can be run under equal budgets;
- accepted-move traces and all experiment metrics are reproducible.
