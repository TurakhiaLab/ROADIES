# MLIPPER for ROADIES

This document describes how ROADIES passes one gene at a time to MLIPPER.

## Current ROADIES integration

ROADIES calls MLIPPER only for GPU-backed placement, selected with
`--mode placement --gpu N` where `N > 0`. Placement runs without `--gpu` (or
with `--gpu 0`) use RAxML-NG instead. The GPU path invokes the MLIPPER
executable built inside the ROADIES checkout; it does not pull or start the
MLIPPER Docker image automatically.

The ROADIES placement command has the following form:

```bash
/path/to/MLIPPER \
  --tree-alignment GENE/ref.fa \
  --query-alignment GENE/query.fa \
  --tree GENE/backbone.nwk \
  --best-model GENE/gene.raxml.bestModel \
  --commit-to-tree GENE/mlipper_gene_tree.nwk
```

The exact executable path is managed by the ROADIES checkout. MLIPPER exits
with status `0` after writing a successful result; ROADIES should also verify
that the output tree is non-empty.

## Per-gene input and output

Each gene requires:

- a reference alignment containing the taxa in the backbone tree;
- a query alignment containing the taxa to insert;
- a backbone tree in Newick format;
- a RAxML-NG `*.raxml.bestModel` file.

The reference and query alignments must have the same alignment length and
must not share sequence names. The output is one committed gene tree in
Newick format.

The current `--best-model` workflow supports four-state DNA GTR models with
discrete Gamma rate heterogeneity. FreeRate, mixture, amino-acid, and non-GTR
models are not supported by this path. The public workflow requires
`pinv = 0.0`.

## Run with Docker instead

Docker is an alternative way to run one gene without installing MLIPPER and
its native dependencies on the host. It is not currently wired into ROADIES
itself.

Pull the published image:

```bash
docker pull wenchiehlo/mlipper:runtime
```

From the directory containing the four input files, run:

```bash
docker run --rm --gpus 'device=0' \
  --user "$(id -u):$(id -g)" \
  -v "$PWD:/workspace/job" \
  -w /workspace/job \
  wenchiehlo/mlipper:runtime \
  --tree-alignment /workspace/job/ref.fa \
  --query-alignment /workspace/job/query.fa \
  --tree /workspace/job/backbone.nwk \
  --best-model /workspace/job/gene.raxml.bestModel \
  --commit-to-tree /workspace/job/mlipper_gene_tree.nwk
```

The bind mount makes the host files available under `/workspace/job` and
writes the result back to the same host directory. Change `device=0` to the
GPU that should run the job.

To use Docker automatically inside ROADIES, the ROADIES placement workflow
would need an explicit container integration, such as a Docker call or a
Snakemake container directive. The current workflow does not provide that
integration.

## Build the runtime image

To rebuild the image from this repository instead of pulling it:

```bash
docker build \
  -f docker/Dockerfile.runtime \
  -t wenchiehlo/mlipper:runtime \
  .
```

The image uses the repository's pinned CUDA base image and libpll-2 `0.4.0`
release, builds MLIPPER in double precision, and packages the executable and
runtime libraries. Running the resulting image still requires an NVIDIA GPU,
a compatible host driver, and the NVIDIA Container Toolkit.

## GPU selection

MLIPPER uses one GPU per invocation. For a local executable, select a device
with `--gpu-id N`, restrict visibility with `CUDA_VISIBLE_DEVICES`, or use
`--gpu-auto` for shared admission across visible GPUs. For Docker, control
which devices enter the container with Docker's `--gpus` option.

See [`src/gpu/README.md`](src/gpu/README.md) for the detailed runtime controls.
