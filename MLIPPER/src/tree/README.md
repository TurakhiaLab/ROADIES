# Tree State and Likelihood Messages

This directory manages MLIPPER's tree topology and the likelihood information attached to it. Start with the pictures and examples below. The API reference at the end is mainly for people changing the implementation.

## The short version

MLIPPER keeps the same tree in three forms:

1. `TreeBuildResult` is the editable CPU tree.
2. `HostPacking` is the same tree rearranged into flat CPU arrays.
3. `DeviceTree` is a view of the arrays stored on the GPU.

The CPU tree is the source of truth for topology. The packed arrays and GPU buffers are derived copies. After changing a parent, child, or branch length, the affected derived data must be refreshed before likelihood is evaluated.

## What is a CLV?

A conditional likelihood vector (CLV) is a compact summary of sequence data. Instead of repeatedly visiting every tip below a node, MLIPPER stores the answer to this question:

> If this node were A, C, G, or T, how likely is all the sequence data on this
> side of the tree?

The vector contains an answer for every compressed alignment site and Gamma rate category:

```text
values in one node CLV = sites × rate categories × states
```

For DNA, `states` is normally four. A tip starts from its observed DNA symbol. An internal node combines the summaries from its two children. A PMAT (transition-probability matrix) accounts for substitutions along each branch.

## Why are there upward and outside messages?

Calling them “up” and “down” can be confusing because the biological tree is unrooted. These names only describe the temporary traversal orientation used by the computation.

```text
                         R              computational root
                        / \
                       X   Y
                      / \ / \
                     A  B C  D

Inside message for X:

                     A ──┐
                         ├──> upward CLV(X)
                     B ──┘

Outside message for X:

                     C ──┐
                         ├──> CLV(Y) ──> root side ──> outside(X)
                     D ──┘

The edge above X separates inside {A, B} from outside {C, D}.
```

For the highlighted edge above node `X`:

- `X`'s upward CLV summarizes everything in the `A/B` subtree.
- The outside message for `X` summarizes everything on the other side: `C/D`.
- Together, the two messages contain the information needed to score that edge.

This is why placement and local topology search can test an edge without recomputing the entire tree.

| Buffer                  | Plain-language meaning                                             |
| ----------------------- | ------------------------------------------------------------------ |
| `d_clv_up`            | Sequence evidence inside each node's descendant subtree.           |
| `d_clv_down`          | Workspace used to pass evidence from the root side toward a child. |
| `d_edge_outside_clv`  | Sequence evidence on the opposite side of an edge.                 |
| `d_edge_midpoint_clv` | A prepared edge-side message used by placement and branch scoring. |
| scaler buffers          | Exponents that keep tiny likelihood values from becoming zero.     |

The root only selects message direction. It is not a claim about the biological root of the phylogeny.

## What is a CLV operation?

A `NodeOpInfo` is one small instruction, not a likelihood value. It tells a GPU kernel which three nodes participate in one update:

```text
parent P with children L and R
```

An upward operation means:

```text
read CLV(L), CLV(R), and their PMATs
                     ↓
               write CLV(P)
```

A downward-left operation means:

```text
read the message above P and the message from sibling R
                     ↓
       write the outside message for child L
```

`DOWN_RIGHT` is the mirror image: it uses sibling `L` to build the outside message for child `R`.

| `dir_tag`            | Target written by the operation | Information combined        |
| ---------------------- | ------------------------------- | --------------------------- |
| `CLV_DIR_UP`         | Parent                          | Left child + right child    |
| `CLV_DIR_DOWN_LEFT`  | Left child's edge               | Parent side + right sibling |
| `CLV_DIR_DOWN_RIGHT` | Right child's edge              | Parent side + left sibling  |

`OP_TIP_TIP`, `OP_TIP_INNER`, and related values select an optimized kernel variant based on whether the inputs are observed tips or computed internal CLVs. They do not represent different mathematical objectives.

## When is a full rebuild necessary?

A CLV is valid only for the topology, branch lengths, PMATs, model parameters, pattern weights, and tip data used to create it.

Use a full rebuild when loading a tree or globally changing the model. During a local SPR evaluation, most of the tree has not changed, so only messages whose paths cross the edited region need to be refreshed.

```text
Full rebuild                         Local SPR refresh

       R                                    R
       #                                    *
       P                                    P
      # #                                  * \
     X   Y                                X   Y
    # # # #                              / * / \
   A  B C  D                            A  B C  D

  # = every message rebuilt             * = changed path rebuilt
                                         other messages are reused
```

The important rule is:

> Reusing an unaffected CLV is an optimization. Reusing a CLV whose represented
> side of the tree changed is a correctness bug.

## Typical workflows

### Loading a tree

```text
parse Newick
  → build TreeBuildResult
  → build HostPacking and PMATs
  → allocate/copy DeviceTree
  → run upward and downward operations
  → likelihood state is ready
```

### Committing a placement or topology change

```text
edit TreeBuildResult links and branch lengths
  → rebuild traversal order
  → update HostPacking
  → refresh changed PMATs
  → prepare/upload affected operations
  → refresh affected CLVs and scalers
  → evaluate likelihood
```

All representations must describe the same tree before the final likelihood step.

### Evaluating local SPR candidates

1. Choose a subtree to prune.
2. Preserve CLVs that summarize unaffected parts of the full tree.
3. Recompute the path changed by suppressing the old attachment.
4. Score candidate regraft edges with their inside and outside messages.
5. Fully validate selected moves before committing them.

`UpdateTreeClvsAfterPrune()` performs step 3. It does not certify that a move is good; it only makes the likelihood state consistent with the temporary pruned topology.

## D&C boundary tips

A D&C sector intentionally omits most of the full tree. At every cut edge, MLIPPER creates a virtual boundary tip whose CLV summarizes the omitted side. It is best understood as a cached message entering the sector, not as an observed species or sequence.

For example, suppose the D&C sector keeps the subtree around `X` but omits the `Y/C/D` side:

```text
full tree

                         R
                        / \
                       X   Y
                      / \ / \
                     A  B C  D

cut the edge between R and the omitted Y side
                                  |
                                  v

sector tree                                      R
                                                / \
                                               X   V
                                              / \
                                             A   B
                                                    ^
                                                    |
                                                    | install V as a virtual tip
    C ──┐                                           |
        ├──> Y contribution ──> boundary CLV   ──>  V
    D ──┘
omitted full-tree component

V receives Y's directional contribution; it is not a new observed sequence.
```

`V` does not contain an observed DNA row. Its installed CLV is the already-computed directional likelihood contribution from the omitted component. To a likelihood kernel inside the sector, `V` behaves like a tip input, allowing the bounded sector to retain the statistical context of the full tree without copying every omitted node.

```text
ordinary tip:           observed characters -> initialize tip CLV
virtual boundary tip:   omitted subtree      -> reuse directional boundary CLV
```

When refreshing a sector workspace, these boundary-tip CLVs must be preserved. Replacing them with ordinary tip characters changes the likelihood problem.

## Tree representations in more detail

### `TreeBuildResult`

The mutable CPU topology contains `TreeNode` objects, a root ID, traversal orders, and tip-name lookup. `TreeNode::id` is the current array index. `stable_node_label` is retained across topology repacking and is used when an earlier result must find the same logical node again.

### Node IDs, traversal order, and stable labels

`TreeNode::id` is a storage index, and the invariant is `tree.nodes[node_id].id == node_id`. It is not a taxon number, traversal position, or permanent biological identifier.

The same logical node may be referred to through several different number spaces:

```text
                              one logical tip A
                                     |
             +-----------------------+-----------------------+
             |                       |                       |
             v                       v                       v
       node_id = 0          stable_node_label = 2      tip_index = 0
   index tree.nodes[]       survives workflow edits    index tip arrays/MSA
             |
             | copied into a D&C sector
             v
      sector_local_id = 5
   index sector.nodes[] and sector GPU buffers

These numbers may differ even though they refer to the same full-tree tip.
```

When MLIPPER first parses a Newick tree, it visits nodes in postorder and assigns consecutive IDs in that order. Children therefore receive smaller IDs than their parent in the initial tree, and the initial computational root is normally the last node:

```text
postorder visit: A, B, X, C, D, Y, R
assigned ID:     0, 1, 2, 3, 4, 5, 6

             R:6
            /   \
          X:2   Y:5
          / \   / \
        A:0 B:1 C:3 D:4
```

That child-before-parent pattern is convenient for the first upward CLV traversal, but code must not depend on it remaining true forever. `preorder` and `postorder` arrays define execution order; numeric node IDs only select array slots.

Committing a query appends two new nodes without renumbering existing nodes:

```text
old target edge: parent ---- target

new IDs:
  internal_id = tree.nodes.size()
  tip_id      = internal_id + 1

parent ---- new internal ---- target
                    |
                 new tip
```

`rebuild_traversals()` recalculates the preorder and postorder arrays after the edit, but it does not change node IDs. If insertion occurs above the current root, the new internal node becomes `root_id`; therefore code must read `root_id` rather than assume the root is `N - 1`.

For example, inserting tip `E` on the branch above `X` appends IDs 7 and 8, while existing node `X` keeps ID 2:

```text
before insertion                    after insertion

        R:6                                R:6
       /   \                              /   \
     X:2   Y:5                         J:7   Y:5
     / \   / \                         / \   / \
   A:0 B:1 C:3 D:4                   X:2 E:8 C:3 D:4
                                     / \
                                   A:0 B:1

new postorder: A:0, B:1, X:2, E:8, J:7, C:3, D:4, Y:5, R:6

The postorder sequence is correct, but its IDs are no longer increasing.
```

`stable_node_label` is a separate logical label. At initial session setup, existing nodes receive labels in preorder. Newly committed internal and tip nodes receive fresh monotonically increasing labels. This lets an earlier result find the same logical node after a workflow rebuilds or remaps storage. Array indexing must still use `id`, not `stable_node_label`.

D&C subtrees introduce another namespace. Selected full-tree nodes are copied into a dense local `subtree.nodes` array and receive local IDs. The plan's `local_to_global` and `global_to_local` mappings translate between the two namespaces. Virtual boundary tips are appended to the local subtree, receive local IDs, and use `stable_node_label = -1` because they represent incoming likelihood messages rather than persistent full-tree nodes.

```text
full-tree node ID  --mapping--> sector-local node ID
       842                         17

Never use 17 to index the full tree or 842 to index the sector workspace.
```

Finally, `tip_node_ids[tip_index]` maps a compact alignment/tip slot to its tree node ID. `tip_index` and `node_id` are not interchangeable because internal nodes occupy slots in `tree.nodes` and tips are not guaranteed to form one permanent contiguous ID range.

### `HostPacking`

Flat host arrays derived from `TreeBuildResult`, including parent/child IDs, branch lengths, tip characters, PMATs, traversal order, and pattern weights. It is not a second independently editable topology.

### `DeviceTree` and `OwnedDeviceTree`

`DeviceTree` is a non-owning descriptor passed to CUDA kernels. It contains raw pointers but does not control their lifetime. `OwnedDeviceTree` is move-only and owns those allocations.

Do not copy an `OwnedDeviceTree`, free memory through a plain `DeviceTree`, or retain a raw pointer across workspace reallocation.

### GPU buffer layout

The four main CLV pools use the same node-major layout:

```text
d_clv_up
d_clv_down
d_edge_midpoint_clv
d_edge_outside_clv

logical shape: [node][site][rate category][state]

flat offset = (((node * sites) + site) * rate_categories + rate) * states + state
values per node = sites * rate_categories * states
```

Keeping the layout identical lets a `NodeOpInfo` select an upward or downward pool without changing site/rate/state indexing. These pointers are separate logical pools; a kernel must use the pointer corresponding to the message meaning described in [Why are there upward and outside messages?](#why-are-there-upward-and-outside-messages).

Scaler layout depends on `per_rate_scaling`:

```text
per_rate_scaling = true:  [node][site][rate category]
per_rate_scaling = false: [node][site]
```

A scaler stores a power-of-two shift, not a likelihood. Kernels add inherited shifts when combining messages and apply the recorded correction when converting a scaled CLV into a log likelihood.

Other important layouts are:

```text
d_tipchars       [tip index][site]
d_tip_node_ids   [tip index] -> node ID
d_query_chars    [query index][site]
d_query_clv      [query index][site][rate category][state]
d_pmat           [child node ID][rate category][parent state][child state]
d_pattern_weights_u[site]
```

`N`, `tips`, and `placement_queries` describe live entries. `capacity_N`, `capacity_tips`, and `query_capacity` describe allocated storage. A kernel may process only live entries even when the backing allocation is larger. Reallocation can invalidate every raw pointer in an earlier `DeviceTree` copy, which is why the owning workspace must outlive all kernel launches that use its view.

## Function guide

Construction and synchronization:

- `build_tree_from_newick_with_pll()` parses a binary Newick tree and chooses a computational traversal root. It creates the [`TreeBuildResult`](#treebuildresult) CPU source of truth.
- `pack_host_arrays_from_tree_and_msa()` converts that editable tree into the flat [`HostPacking`](#hostpacking) representation.
- `fill_pmats_in_host_packing()` refreshes the branch transition matrices used when [child CLVs are combined](#what-is-a-clv).
- `allocate_device_tree_on_current_gpu()` creates the owned GPU allocation described under [`DeviceTree` and `OwnedDeviceTree`](#devicetree-and-owneddevicetree).
- `reload_device_tree_live_data*()` synchronizes an existing GPU allocation after the CPU topology or branch data changes; see the [topology-change workflow](#committing-a-placement-or-topology-change).
- `load_subtree_workspace()` creates or reuses the bounded GPU workspace used by [local SPR](#evaluating-local-spr-candidates) and [D&amp;C sectors](#dc-boundary-tips).

CLV updates:

- `PrepareTreeClvOperations()` builds the reusable, traversal-ordered list of [`NodeOpInfo` operations](#what-is-a-clv-operation) described above. Each entry says which parent and children to read and whether to perform an `UP`, `DOWN_LEFT`, or `DOWN_RIGHT` update.
- `UpdateTreeClvs()` performs the [full rebuild](#when-is-a-full-rebuild-necessary) used when every upward and outside message must be regenerated.
- `UpdateTreeClvsPreservingTipClvs()` performs a rebuild while preserving [D&amp;C virtual boundary messages](#dc-boundary-tips).
- `UpdateTreeClvsUpwardOnly()` runs only the `UP` operations explained in [the NodeOp direction table](#what-is-a-clv-operation). This is enough when only root likelihood is needed.
- `UpdateTreeClvsPrepared()` executes a previously prepared list of [`NodeOpInfo` operations](#what-is-a-clv-operation).
- `UpdateTreeClvsAfterPrune()` performs the [local refresh](#evaluating-local-spr-candidates) required after temporarily detaching a subtree.
- `BuildSingleTreeEdgeOutsideWarpSitePrepared()`, `UpdateSingleTreeClvUpwardWarpSitePrepared()`, and `RefreshSingleTreeChildDownPrepared()` perform one `outside`, `UP`, or `DOWN` message update at a time. These are the single-edge form of the [directional operations](#why-are-there-upward-and-outside-messages).

## How Tree Work Is Assigned on the GPU

Tree orchestration decides which NodeOps are independent; likelihood kernels decide how threads process their sites. For a full ordered traversal, one site thread executes every NodeOp in order and writes each operation's target CLV slice. For a levelized traversal, the host launches one dependency level at a time, `grid.y` selects an independent NodeOp in that level, and x-threads divide its sites.

```text
one NodeOp + one site
  reads:  source CLV slices, source scalers, and branch PMATs
  writes: target CLV[site][rate][state] and target scaler[site]
```

Tip initialization flattens `(tip, site)` pairs and assigns one pair to each thread; its output is the decoded rate/state CLV slice and a zero scaler for that tip/site. Query construction launches one selected query at a time, assigns one thread to each site, and writes that query's complete rate/state CLV slice. Copy kernels used by subtree workspaces assign threads across flattened CLV/scaler elements and write the corresponding destination workspace; they do not perform likelihood reductions.

Direct NNI context construction uses `grid.y` for the requested target-context operation and one x-thread per site. That thread copies or transforms the complete rate/state message for the site and writes the target upward CLV, edge-outside CLV, and scalers in the destination sector workspace. State loops and scaling maxima are local to that thread; this kernel has no cross-thread reduction.

The specialized DNA+G4 single-operation path uses 16 lanes per site, corresponding to four rates by four states. It uses four-lane maximum reductions to choose each rate category's scaling shift, then writes the 16 CLV components. See [GPU Parallelism and Outputs](../likelihood/README.md#gpu-parallelism-and-outputs) for the complete likelihood and reduction map.

Placement and D&C:

- `build_placement_query()` turns alignment rows into the query records used by the [placement workflow](#committing-a-placement-or-topology-change).
- `UploadPlacementOps()` uploads candidate-edge [`NodeOpInfo` instructions](#what-is-a-clv-operation) to the GPU.
- `EvaluatePlacementQueries()` uses the [inside and outside edge messages](#why-are-there-upward-and-outside-messages) to score queries and optionally commit them.
- `build_subtree_plan_from_anchors()` and `build_nni_sector_plan()` select a bounded part of the tree for the [D&amp;C workflow](#dc-boundary-tips).
- `build_subtree_from_plan()` constructs the sector's local [`TreeBuildResult`](#treebuildresult).
- `prepare_*_with_directional_boundaries()` attaches the [virtual boundary messages](#dc-boundary-tips) that represent the omitted full-tree context.

## Files

- `tree.hpp`: tree/state structures and orchestration declarations.
- `tree_generation.cpp`: Newick parsing and CPU tree construction.
- `tree_generation_device.cu`: GPU allocation, CLV updates, and placement commits.
- `tree_topology_utils.hpp`: traversal rebuilding and branch-optimization order.
- `divide_and_conquer.hpp/.cpp`: sector selection, extraction, and boundaries.

When debugging a likelihood mismatch after a topology edit, check the synchronization workflow first. Most failures come from one stale representation or directional message rather than from the tree link edit.
