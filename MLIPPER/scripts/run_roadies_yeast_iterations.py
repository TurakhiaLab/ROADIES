#!/usr/bin/env python3
"""Run MLIPPER over one or more ROADIES yeast iterations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import queue
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
DEFAULT_DATA_ROOT = Path("/data/w4lo/MLIPPER_data/yeast")
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "output" / "roadies_yeast_iterations"


@dataclass(frozen=True)
class GeneJob:
    name: str
    ref_msa: Path
    query_msa: Path
    backbone_tree: Path
    best_model: Path


def gene_number(name: str) -> int:
    try:
        return int(name.removeprefix("gene_"))
    except ValueError as exc:
        raise RuntimeError(f"Unexpected gene directory name: {name}") from exc


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_jobs(data_root: Path, iteration: int) -> list[GeneJob]:
    placement_root = data_root / f"iter_{iteration}_placement" / "genes"
    backbone_root = data_root / f"iter_{iteration}_backbone" / "genes"
    if not placement_root.is_dir() or not backbone_root.is_dir():
        raise RuntimeError(f"Iteration {iteration} input directories are missing")

    jobs: list[GeneJob] = []
    gene_dirs = [path for path in placement_root.iterdir() if path.is_dir()]
    for gene_dir in sorted(gene_dirs, key=lambda path: gene_number(path.name)):
        gene = gene_dir.name
        job = GeneJob(
            name=gene,
            ref_msa=gene_dir / "iter0_output_msa_from_ref.fa",
            query_msa=gene_dir / "iter0_output_msa_from_query.fa",
            backbone_tree=backbone_root / f"{gene}_filtered.fa.aln.raxml.bestTree",
            best_model=gene_dir / "iter1_tree_output" / "gene_tree.raxml.bestModel",
        )
        missing = [path for path in vars(job).values() if isinstance(path, Path) and not path.is_file()]
        if missing:
            missing_text = ", ".join(str(path) for path in missing)
            raise RuntimeError(f"{gene} is incomplete; missing: {missing_text}")
        jobs.append(job)
    if not jobs:
        raise RuntimeError(f"Iteration {iteration} contains no runnable genes")
    return jobs


def write_manifest(path: Path, jobs: list[GeneJob]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "ref_msa", "query_msa", "backbone_tree", "best_model"])
        for job in jobs:
            writer.writerow(
                [job.name, job.ref_msa, job.query_msa, job.backbone_tree, job.best_model]
            )


def mlipper_command(
    binary: Path,
    job: GeneJob,
    output_tree: Path,
    gpu: int,
    args: argparse.Namespace,
) -> list[str]:
    return [
        str(binary),
        "--no-model-optimization",
        "--no-global-branch-optimization",
        "--tree-alignment",
        str(job.ref_msa),
        "--query-alignment",
        str(job.query_msa),
        "--tree",
        str(job.backbone_tree),
        "--best-model",
        str(job.best_model),
        "--empirical-freqs",
        "--commit-to-tree",
        str(output_tree),
        "--batch-insert-size",
        str(args.batch_size),
        "--local-spr-radius",
        str(args.local_spr_radius),
        "--local-spr-cluster-threshold",
        str(args.local_spr_cluster_threshold),
        "--local-spr-rounds",
        str(args.local_spr_rounds),
        "--gpu-id",
        str(gpu),
    ]


def run_iteration(
    iteration: int,
    jobs: list[GeneJob],
    binary: Path,
    output_root: Path,
    args: argparse.Namespace,
) -> dict[str, object]:
    iteration_dir = output_root / f"iter_{iteration}"
    genes_dir = iteration_dir / "genes"
    genes_dir.mkdir(parents=True, exist_ok=True)
    write_manifest(iteration_dir / "manifest.tsv", jobs)

    pending: queue.Queue[GeneJob] = queue.Queue()
    skipped: list[str] = []
    for job in jobs:
        tree = genes_dir / job.name / "mlipper_gene_tree.nwk"
        if args.resume and tree.is_file() and tree.stat().st_size > 0:
            skipped.append(job.name)
        else:
            pending.put(job)

    failures: list[tuple[str, int, int]] = []
    failure_lock = threading.Lock()
    print_lock = threading.Lock()
    completed = 0
    completed_lock = threading.Lock()
    stop_requested = threading.Event()
    active_processes: dict[int, subprocess.Popen[str]] = {}
    active_processes_lock = threading.Lock()

    def worker(gpu: int) -> None:
        nonlocal completed
        while not stop_requested.is_set():
            try:
                job = pending.get_nowait()
            except queue.Empty:
                return
            gene_dir = genes_dir / job.name
            gene_dir.mkdir(parents=True, exist_ok=True)
            output_tree = gene_dir / "mlipper_gene_tree.nwk"
            log_path = gene_dir / "mlipper.log"
            output_tree.unlink(missing_ok=True)
            started = time.monotonic()
            with log_path.open("w", encoding="utf-8") as log:
                process = subprocess.Popen(
                    mlipper_command(binary, job, output_tree, gpu, args),
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
                with active_processes_lock:
                    active_processes[gpu] = process
                returncode = process.wait()
                with active_processes_lock:
                    active_processes.pop(gpu, None)
            elapsed = time.monotonic() - started
            if returncode != 0 or not output_tree.is_file() or output_tree.stat().st_size == 0:
                output_tree.unlink(missing_ok=True)
                with failure_lock:
                    failures.append((job.name, gpu, returncode))
                status = "FAIL"
            else:
                with completed_lock:
                    completed += 1
                    done = completed
                status = f"DONE {done + len(skipped)}/{len(jobs)}"
            with print_lock:
                print(
                    f"[iter {iteration}][gpu {gpu}] {status} {job.name} ({elapsed:.1f}s)",
                    flush=True,
                )
            pending.task_done()

    threads = [threading.Thread(target=worker, args=(gpu,), daemon=False) for gpu in args.gpus]
    for thread in threads:
        thread.start()
    try:
        for thread in threads:
            while thread.is_alive():
                thread.join(timeout=0.5)
    except KeyboardInterrupt:
        stop_requested.set()
        with active_processes_lock:
            processes = list(active_processes.values())
        for process in processes:
            process.terminate()
        for process in processes:
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill()
        for thread in threads:
            thread.join()
        raise

    failures.sort(key=lambda row: gene_number(row[0]))
    with (iteration_dir / "failed.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "gpu", "exit_code"])
        writer.writerows(failures)

    successful_jobs = [
        job
        for job in jobs
        if (genes_dir / job.name / "mlipper_gene_tree.nwk").is_file()
        and (genes_dir / job.name / "mlipper_gene_tree.nwk").stat().st_size > 0
    ]
    aggregate = iteration_dir / "mlipper_all_gene_trees.tre"
    with aggregate.open("wb") as output:
        for job in successful_jobs:
            tree = genes_dir / job.name / "mlipper_gene_tree.nwk"
            content = tree.read_bytes()
            output.write(content)
            if content and not content.endswith(b"\n"):
                output.write(b"\n")

    comparison_status = "not_checked"
    matches = 0
    mismatches: list[str] = []
    missing_references: list[str] = []
    if args.reference_root is not None:
        comparison_status = "checked"
        reference_genes = args.reference_root / f"iter_{iteration}" / "genes"
        with (iteration_dir / "byte_comparison.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["gene", "status", "new_sha256", "reference_sha256"])
            for job in successful_jobs:
                new_tree = genes_dir / job.name / "mlipper_gene_tree.nwk"
                reference_tree = reference_genes / job.name / "mlipper_gene_tree.nwk"
                if not reference_tree.is_file():
                    missing_references.append(job.name)
                    writer.writerow([job.name, "missing_reference", sha256(new_tree), ""])
                elif new_tree.read_bytes() == reference_tree.read_bytes():
                    matches += 1
                    digest = sha256(new_tree)
                    writer.writerow([job.name, "match", digest, digest])
                else:
                    mismatches.append(job.name)
                    writer.writerow(
                        [job.name, "mismatch", sha256(new_tree), sha256(reference_tree)]
                    )

    summary: dict[str, object] = {
        "iteration": iteration,
        "expected_genes": len(jobs),
        "successful_genes": len(successful_jobs),
        "failed_genes": len(failures),
        "resumed_genes": len(skipped),
        "aggregate_sha256": sha256(aggregate),
        "byte_comparison": comparison_status,
        "byte_matches": matches,
        "byte_mismatches": mismatches,
        "missing_references": missing_references,
    }
    (iteration_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run complete ROADIES yeast iterations with one worker per GPU."
    )
    parser.add_argument("--iterations", type=int, nargs="+", default=[2, 3, 4, 5])
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--binary", type=Path, default=REPO_ROOT / "MLIPPER")
    parser.add_argument("--gpus", type=int, nargs="+", default=[0, 1])
    parser.add_argument("--reference-root", type=Path)
    parser.add_argument("--batch-size", type=int, default=5)
    parser.add_argument("--local-spr-radius", type=int, default=4)
    parser.add_argument("--local-spr-cluster-threshold", type=int, default=3)
    parser.add_argument("--local-spr-rounds", type=int, default=1)
    parser.add_argument("--resume", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument(
        "--dry-run", action="store_true", help="Validate and report inputs without running MLIPPER"
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.data_root = args.data_root.resolve()
    args.output_root = args.output_root.resolve()
    args.binary = args.binary.resolve()
    if args.reference_root is not None:
        args.reference_root = args.reference_root.resolve()
    if len(set(args.gpus)) != len(args.gpus) or any(gpu < 0 for gpu in args.gpus):
        raise SystemExit("--gpus must contain unique, non-negative device IDs")
    if len(set(args.iterations)) != len(args.iterations):
        raise SystemExit("--iterations must not contain duplicates")
    if any(iteration < 1 for iteration in args.iterations):
        raise SystemExit("--iterations values must be positive")
    if not args.gpus:
        raise SystemExit("At least one GPU is required")
    if args.reference_root is not None and not args.reference_root.is_dir():
        raise SystemExit(f"Reference root does not exist: {args.reference_root}")
    if args.batch_size < 1 or args.local_spr_radius < 1:
        raise SystemExit("Batch size and Local SPR radius must be positive")
    if args.local_spr_cluster_threshold < 0 or args.local_spr_rounds < 1:
        raise SystemExit("Cluster threshold must be non-negative and rounds must be positive")

    jobs_by_iteration = {
        iteration: discover_jobs(args.data_root, iteration) for iteration in args.iterations
    }
    for iteration, jobs in jobs_by_iteration.items():
        print(f"Iteration {iteration}: {len(jobs)} genes")
    if args.dry_run:
        return 0
    if not args.binary.is_file() or not os.access(args.binary, os.X_OK):
        raise SystemExit(f"MLIPPER binary is missing or not executable: {args.binary}")

    args.output_root.mkdir(parents=True, exist_ok=True)
    run_metadata = {
        "binary": str(args.binary),
        "binary_sha256": sha256(args.binary),
        "git_head": subprocess.run(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip(),
        "iterations": args.iterations,
        "gpus": args.gpus,
        "batch_size": args.batch_size,
        "local_spr_radius": args.local_spr_radius,
        "local_spr_cluster_threshold": args.local_spr_cluster_threshold,
        "local_spr_rounds": args.local_spr_rounds,
        "reference_root": str(args.reference_root) if args.reference_root else None,
    }
    (args.output_root / "run_metadata.json").write_text(
        json.dumps(run_metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    summaries = []
    for iteration in args.iterations:
        summaries.append(
            run_iteration(
                iteration,
                jobs_by_iteration[iteration],
                args.binary,
                args.output_root,
                args,
            )
        )
    (args.output_root / "summary.json").write_text(
        json.dumps(summaries, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 1 if any(summary["failed_genes"] for summary in summaries) else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("Interrupted; completed non-empty gene trees can be resumed.", file=sys.stderr)
        raise SystemExit(130)
