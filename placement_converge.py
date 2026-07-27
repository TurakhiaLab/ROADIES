import csv
import os
import subprocess
import time
import yaml
import argparse
from pathlib import Path
import shutil

ROADIES_ROOT = Path(__file__).resolve().parent


def parse_args():
    parser = argparse.ArgumentParser(
        description="Iterative ROADIES backbone + placement pipeline"
    )
    parser.add_argument(
        "--backbone",
        required=True,
        help="Path to directory of backbone/reference genome files",
    )
    parser.add_argument(
        "--query",
        required=True,
        help="Path to directory of query genome files",
    )
    parser.add_argument(
        "--config",
        default=str(ROADIES_ROOT / "config" / "config.yaml"),
        help="Path to ROADIES config.yaml (default: config/config.yaml next to this script)",
    )
    parser.add_argument(
        "--out-base-dir",
        default=str(ROADIES_ROOT / "roadies_iterations"),
        help="Base output directory for all iterations (default: roadies_iterations/ next to this script)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=48,
        help="Number of cores (default: 48)",
    )
    parser.add_argument(
        "--support-threshold",
        type=float,
        default=0.95,
        help="Local posterior probability threshold for convergence (default: 0.95)",
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=9,
        help="Maximum number of backbone+placement iterations (default: 9)",
    )
    parser.add_argument(
        "--gpu",
        type=int,
        default=0,
        help="Number of GPU devices to use (default: 0, CPU only)",
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume from the last completed iteration (reads time_stamps.csv)",
    )
    parser.add_argument(
        "--deep",
        action="store_true",
        help="Enable deep phylogeny mode (applies to both the backbone and placement stages)",
    )
    return parser.parse_args()


def run_roadies(roadies_script, mode, config_file, cores, gpu, no_clean=False, deep=False):
    cmd = [
        "python3", roadies_script,
        "--mode", mode,
        "--config", config_file,
        "--cores", str(cores),
        "--noconverge",
        "--gpu", str(gpu),
    ]
    if no_clean:
        cmd.append("--no-clean")
    if deep:
        cmd.append("--deep")
    print(f"Running ROADIES: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def update_gene_count(config_file, iteration, base_gene_count=1000):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    if iteration <= 2:
        gene_count = base_gene_count
    else:
        gene_count = base_gene_count * (2 ** (iteration - 2))

    config["GENE_COUNT"] = gene_count

    with open(config_file, "w") as f:
        yaml.safe_dump(config, f)

    print(f"[ITER {iteration}] GENE_COUNT set to {gene_count}")


def compute_percent_high_support(freq_file, threshold):
    count = 0
    with open(freq_file, "r") as file:
        csv_reader = csv.reader(file, delimiter="\t")
        rows = list(csv_reader)
        total_rows = len(rows)
        for i, row in enumerate(rows):
            if (i + 1) % 3 == 1:
                value = float(row[3])
                if value >= threshold:
                    count += 1
    percent_high_support = (count / (total_rows / 3)) * 100
    return percent_high_support


def update_config_yaml(config_file, out_dir=None, species=None, ref_dir=None):
    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    if out_dir:
        config["OUT_DIR"] = out_dir
        Path(out_dir).mkdir(parents=True, exist_ok=True)
    if species:
        config["GENOMES"] = species
    if ref_dir:
        config["REF_DIR"] = ref_dir

    with open(config_file, "w") as f:
        yaml.safe_dump(config, f)


def combine_iter(out_dir, run, cores, out_base_dir, roadies_dir):
    run_dir = Path(out_dir) / run
    astral_out = run_dir / f"{run}.nwk"
    astral_stats_out = run_dir / f"{run}_stats.nwk"

    if not astral_out.exists():
        master_gt = Path(out_dir) / "master_gt.nwk"
        master_map = Path(out_dir) / "master_map.txt"
        append_marker = run_dir / f"{run}.append_done"

        if not append_marker.exists():
            with open(run_dir / "genetrees/gene_tree_merged.nwk") as infile, open(master_gt, "a") as outfile:
                outfile.write(infile.read())
            with open(run_dir / "genes/mapping_combined.txt") as infile, open(master_map, "a") as outfile:
                outfile.write(infile.read())
            append_marker.touch()

        subprocess.run(["astral-pro3", "-t", str(cores), "-i", str(master_gt),
                        "-o", str(astral_out), "-a", str(master_map)], check=True, cwd=str(ROADIES_ROOT))
        subprocess.run(["astral-pro3", "-t", str(cores), "-u", "3", "-i", str(master_gt),
                        "-o", str(astral_stats_out), "-a", str(master_map)], check=True, cwd=str(ROADIES_ROOT))

    shutil.copy(astral_out, Path(out_base_dir) / roadies_dir / "roadies.nwk")
    shutil.copy(astral_stats_out, Path(out_base_dir) / roadies_dir / "roadies_stats.nwk")


def sampling_output_is_from_query(query_dir):
    """Return True if sampling_output.txt exists and its first sample belongs to the query directory."""
    path = ROADIES_ROOT / "sampling_output.txt"
    if not path.exists():
        return False
    query_path = Path(query_dir)
    section = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                section = line[1:].strip()
            elif section == "SAMPLES" and line:
                return bool(list(query_path.glob(f"{line}.fa*")))
    return False


def find_resume_point(out_base_dir):
    """Determine the next iteration to run and restore high_support_list from time_stamps.csv.

    A timestamp row is written only after combine_iter completes, so iterations
    with backbone-done-but-placement-not-done will have no entry and are correctly
    returned as the start iteration (Snakemake skips already-done backbone steps).
    """
    timestamps_file = Path(out_base_dir) / "time_stamps.csv"
    if not timestamps_file.exists():
        return 1, []

    high_support_list = []
    last_iter = 0
    with open(timestamps_file, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) >= 2:
                iter_num = int(parts[0])
                pct = float(parts[1])
                last_iter = max(last_iter, iter_num)
                while len(high_support_list) < iter_num:
                    high_support_list.append(0.0)
                high_support_list[iter_num - 1] = pct

    return last_iter + 1, high_support_list


def main():
    args = parse_args()

    # Resolve every user-supplied path against the caller's cwd now, before any
    # subprocess inherits/changes cwd elsewhere, so this script works from any directory.
    args.config = os.path.abspath(args.config)
    args.out_base_dir = os.path.abspath(args.out_base_dir)
    args.backbone = os.path.abspath(args.backbone)
    args.query = os.path.abspath(args.query)

    roadies_script = str(ROADIES_ROOT / "run_roadies.py")
    config_file = args.config
    out_base_dir = args.out_base_dir
    backbone_species = args.backbone
    query_species = args.query
    support_thr = args.support_threshold
    max_iterations = args.max_iterations
    cores = args.cores
    gpu = args.gpu
    deep = args.deep
    roadies_dir = "roadies_final"

    Path(out_base_dir, roadies_dir).mkdir(parents=True, exist_ok=True)

    if args.resume:
        iteration, high_support_list = find_resume_point(out_base_dir)
        print(f"Resuming from iteration {iteration} (found {iteration - 1} completed iteration(s))")
    else:
        iteration = 1
        high_support_list = []

    time_stamps = [time.time()]

    while iteration <= max_iterations:
        print(f"\n=== ITERATION {iteration} ===")

        update_gene_count(config_file, iteration, base_gene_count=1000)

        # --- Backbone ---
        backbone_out_dir = f"{out_base_dir}/iter_{iteration}_backbone"
        backbone_done = Path(backbone_out_dir, "roadies.nwk").exists()
        # Check for partial backbone before update_config_yaml creates the dir
        backbone_partial = (
            not backbone_done
            and Path(backbone_out_dir).exists()
            and any(Path(backbone_out_dir).iterdir())
        )
        if backbone_done:
            print(f"[ITER {iteration}] Backbone already complete, skipping.")
        else:
            (ROADIES_ROOT / "sampling_output.txt").unlink(missing_ok=True)
            update_config_yaml(config_file, out_dir=backbone_out_dir, species=backbone_species, ref_dir=None)
            run_roadies(roadies_script, mode="accurate", config_file=config_file,
                        cores=cores, gpu=gpu, no_clean=backbone_partial, deep=deep)

        # --- Placement ---
        placement_out_dir = f"{out_base_dir}/iter_{iteration}_placement"
        update_config_yaml(config_file, out_dir=placement_out_dir, species=query_species, ref_dir=backbone_out_dir)
        placement_done = Path(placement_out_dir, "genetrees", "gene_tree_merged.nwk").exists()
        if placement_done:
            print(f"[ITER {iteration}] Placement already complete, skipping to combine.")
        elif sampling_output_is_from_query(query_species):
            print(f"[ITER {iteration}] Resuming mid-placement (query sampling_output.txt found).")
            run_roadies(roadies_script, mode="placement", config_file=config_file,
                        cores=cores, gpu=gpu, no_clean=True, deep=deep)
        else:
            (ROADIES_ROOT / "sampling_output.txt").unlink(missing_ok=True)
            run_roadies(roadies_script, mode="placement", config_file=config_file,
                        cores=cores, gpu=gpu, no_clean=False, deep=deep)

        combine_iter(out_base_dir, f"iter_{iteration}_placement", cores, out_base_dir, roadies_dir)

        freq_file = str(ROADIES_ROOT / "freqQuad.csv")
        percent_high_support = compute_percent_high_support(freq_file, support_thr)
        print(f"Iteration {iteration}: Percent high support = {percent_high_support:.2f}%")

        curr_time = time.time()
        curr_time_l = time.asctime(time.localtime(curr_time))
        elapsed_time = curr_time - time_stamps[0]
        time_stamps.append(curr_time)
        high_support_list.append(percent_high_support)

        with open(Path(out_base_dir) / "time_stamps.csv", "a") as t_out:
            t_out.write(f"{iteration},{percent_high_support:.2f},{curr_time_l},{elapsed_time:.2f}\n")

        if ((iteration == 1 and percent_high_support == 100) or
                (iteration >= 2 and
                 (abs(percent_high_support - high_support_list[iteration - 2]) < 1
                  or percent_high_support == 100
                  or iteration == max_iterations))):
            print("Convergence criteria met. Stopping iterations.")
            break

        iteration += 1


if __name__ == "__main__":
    main()
