# MLIPPER

MLIPPER is a GPU-accelerated phylogenetics program for constructing or
extending gene trees under a GTR+Gamma likelihood model.

## Workflows

| Goal                                       | Mode                     | Output      |
| ------------------------------------------ | ------------------------ | ----------- |
| Add query taxa to a backbone tree          | Small-tip commit         | Newick tree |
| Score placements without changing the tree | Placement only           | jplace      |
| Construct and refine a full tree           | Divide and conquer (D&C) | Newick tree |

MLIPPER requires an NVIDIA GPU. It builds gene trees; species-tree inference
and normalized Robinson-Foulds calculations are downstream tasks.

## Docker

Docker is the quickest way to obtain the validated CUDA and library
environment. The host still needs an NVIDIA driver, Docker Engine, and NVIDIA
Container Toolkit.

| Image                                      | Contents                                     | Use |
| ------------------------------------------ | -------------------------------------------- | --- |
| `wenchiehlo/mlipper:runtime` (also `latest`) | MLIPPER and required runtime libraries only | Running MLIPPER; recommended for most users |
| `wenchiehlo/mlipper:dev`                    | MLIPPER source, CUDA build tools, compilers, EPA-ng, RAxML-ng, and ASTRAL-Pro3 | Rebuilding MLIPPER and running development or downstream analysis |

The `latest` tag points to the same image as `runtime`. Use `dev` only when
you need to compile inside the container or use its additional analysis tools;
it is substantially larger than the runtime image.

Pull the published images:

```bash
docker pull wenchiehlo/mlipper:runtime
docker pull wenchiehlo/mlipper:dev
```

For the recommended runtime image, the shorter command is equivalent:

```bash
docker pull wenchiehlo/mlipper
```

Docker also pulls a missing image automatically when it is first run. Check
GPU visibility with:

```bash
docker run --rm --gpus 'device=0' \
  wenchiehlo/mlipper:runtime --help
```

To rebuild an image from the current checkout, run one of these commands from
the repository root:

```bash
docker build -f docker/Dockerfile.runtime \
  -t wenchiehlo/mlipper:runtime .

docker build -f docker/Dockerfile \
  -t wenchiehlo/mlipper:dev .
```

Both Dockerfiles download their CUDA base image and dependencies during the
build. See [docker/README.md](docker/README.md) for advanced options.

## Build from Source

The maintained source build uses Ubuntu 22.04, CUDA 12, C++17, Boost, TBB,
BLAS/LAPACK, zlib, and libpll-2 `0.4.0`. DIPPER is included in
`third_party/dipper`.

```bash
make -j4 MLIPPER
./MLIPPER --help
```

## Small-Tip Placement

The following example uses the included ROADIES iteration 4, gene 1561
fixture. Run it from the repository root:

```bash
mkdir -p data/small_short

./MLIPPER \
  --tree-alignment data/small_tip/small_short/reference.fa \
  --query-alignment data/small_tip/small_short/query.fa \
  --tree data/small_tip/small_short/backbone.nwk \
  --best-model data/small_tip/small_short/model.bestModel \
  --commit-to-tree data/small_short/committed_tree.nwk \
  --gpu-id 0
```

To run the same gene with the runtime image:

```bash
mkdir -p data/small_short

docker run --rm --gpus 'device=0' \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace/job" -w /workspace/job \
  wenchiehlo/mlipper:runtime \
  --tree-alignment /workspace/job/data/small_tip/small_short/reference.fa \
  --query-alignment /workspace/job/data/small_tip/small_short/query.fa \
  --tree /workspace/job/data/small_tip/small_short/backbone.nwk \
  --best-model /workspace/job/data/small_tip/small_short/model.bestModel \
  --commit-to-tree /workspace/job/data/small_short/committed_tree.nwk
```

## Placement Only

```bash
mkdir -p data/small_short

./MLIPPER \
  --tree-alignment data/small_tip/small_short/reference.fa \
  --query-alignment data/small_tip/small_short/query.fa \
  --tree data/small_tip/small_short/backbone.nwk \
  --best-model data/small_tip/small_short/model.bestModel \
  --jplace-out data/small_short/placements.jplace \
  --filter-acc-lwr 0.99 \
  --gpu-id 0
```

## Divide and Conquer

The 10K alignment is compressed in Git. Extract it and run the included
DIPPER fixture from the repository root:

```bash
mkdir -p data/dipper_10k
gzip -dc data/divide_and_conquer/10k/alignment.fa.gz \
  > data/dipper_10k/alignment.fa

./MLIPPER \
  --divide-and-conquer \
  --tree-alignment data/dipper_10k/alignment.fa \
  --best-model data/divide_and_conquer/10k/model.bestModel \
  --empirical-freqs \
  --dipper-starting-tree-mode divide-and-conquer \
  --write-tree data/dipper_10k/mlipper_tree.nwk \
  --gpu-id 0
```

`--dipper-starting-tree-mode` accepts `nj-placement` or
`divide-and-conquer`.

## Included Datasets

The repository includes ROADIES small-tip inputs and the original DIPPER
AliSim 10K and 20K aligned datasets under [`data`](data). See the data README
for the tip/site sizes and provenance.

The commands above extract the 10K alignment into the Git-ignored
`data/dipper_10k` output directory. To use the 20K dataset instead:

```bash
mkdir -p data/dipper_20k
gzip -dc data/divide_and_conquer/20k/alignment.fa.gz \
  > data/dipper_20k/alignment.fa
```

Then replace `10k` with `20k` in the D&C command.

## Inputs and GPU Selection

| Input                 | Small-tip           | D&C            |
| --------------------- | ------------------- | -------------- |
| `--tree-alignment`  | Reference alignment | Full alignment |
| `--query-alignment` | Query alignment     | Not used       |
| `--tree`            | Backbone tree       | Not used       |
| `--best-model`      | Required            | Required       |

Small-tip reference taxa must match the backbone-tree taxa. Reference and
query alignments must have the same aligned length and disjoint taxon names.
The maintained model scope is nucleotide GTR with 1–8 equal-weight Gamma
categories and no invariant-site component.

Use `--gpu-id N` to select a visible CUDA device or `--gpu-auto` to let MLIPPER
select one. Docker GPU ordinals follow the devices exposed to the container.

## Documentation

- [Starter Guide](STARTER_GUIDE.md) — architecture and codebase tour
- [ROADIES Integration](README_ROADIES.md) — per-gene input and invocation guide
- [Source Guide](src/README.md) — implementation map
- [Model Scope](docs/reference/supported_model_scope.md) — supported models
- [Docker Guide](docker/README.md) — image details and GPU visibility
