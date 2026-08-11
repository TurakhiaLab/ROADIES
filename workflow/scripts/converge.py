# converge.py is a script that iteratively runs ROADIES, wherein after each run, the resultant gene trees are concatenated into a master file, and input into ASTRAL-PRO.
# The program stops after configured number of ITERATIONS

# REQUIREMENTS: Activated conda environment with snakemake and ete3

import os, sys, glob
import argparse
import random
import subprocess
import signal
import shutil
from ete3 import Tree
from reroot import rerootTree
import yaml
from pathlib import Path
import time
import math
import csv


# function that finds the average distance between an array of trees and itself
def comp_tree(t1, t2):
    d = t1.compare(t2)
    return d["norm_rf"]


# Function to update the configuration file
def update_config(config_path, base_gene_count):
    with open(config_path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    # Update GENE_COUNT based on the iteration number
    config["GENE_COUNT"] = base_gene_count * 2

    # Save the updated configuration
    with open(config_path, "w") as file:
        yaml.dump(config, file)


# Function to read the initial GENE_COUNT from the config file
def read_initial_gene_count(config_path):
    with open(config_path) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    return config["GENE_COUNT"]


def format_run(iteration):
    """iteration_00, iteration_01, ..., iteration_10, iteration_11, ... -
    shared by converge_run/combine_iter/resume-scanning so they can never
    disagree on what an iteration's directory/file names look like."""
    return "iteration_" + str(iteration).zfill(2)


def find_resume_point(out_dir):
    """Scan out_dir for iterations that actually finished (i.e. have a final
    {run}.nwk ASTRAL output - the last file combine_iter writes, and the
    first thing converge_run reads back afterwards, so its presence proves
    that whole iteration completed) and reconstruct where to pick back up.

    Stops at the first gap: iterations always run strictly in order, so a
    missing {run}.nwk means that iteration never finished, regardless of
    what higher-numbered directories might exist from an earlier attempt.
    """
    percent_by_iteration = {}
    ts_path = os.path.join(out_dir, "time_stamps.csv")
    if os.path.exists(ts_path):
        with open(ts_path) as f:
            for line in f:
                parts = line.strip().split(",")
                if len(parts) < 3:
                    continue  # "Start time: ..." header lines
                try:
                    percent_by_iteration[int(parts[0])] = float(parts[2])
                except ValueError:
                    continue

    high_support_list = []
    iteration = 0
    while (
        os.path.exists(os.path.join(out_dir, format_run(iteration) + ".nwk"))
        and iteration in percent_by_iteration
    ):
        high_support_list.append(percent_by_iteration[iteration])
        iteration += 1
    return iteration, high_support_list


# Rules that spawn one Snakemake job per sampled locus (up to GENE_COUNT
# jobs), keyed by mode. Under --cluster these get grouped so they don't each
# turn into a separate sbatch submission.
PER_LOCUS_RULES = {
    "placement": ["pasta"],
    "accurate": ["pasta", "filtermsa", "raxmlng"],
    "balanced": ["pasta", "filtermsa", "fasttree"],
}
CLUSTER_GROUP_SIZE = 250  # loci per sbatch submission

# Multi-node SLURM execution args, opt-in via --cluster. Mirrors the resource
# profile used on this cluster for the fish-only 32k backbone run: 32
# concurrent sbatch jobs, each 1 node / 16 tasks / 4 cpus-per-task (64 cpus).
def cluster_snakemake_args(mode):
    args = ["--jobs", "32", "--latency-wait", "120", "--keep-going"]
    per_locus_rules = PER_LOCUS_RULES.get(mode, [])
    if per_locus_rules:
        args += ["--groups"] + [f"{rule}=group0" for rule in per_locus_rules]
        args += ["--group-components", f"group0={CLUSTER_GROUP_SIZE}"]
    args += [
        "--cluster",
        (
            "sbatch "
            "--job-name=ROADIES_run "
            "--partition=long "
            "--account=standard "
            "--nodes=1 "
            "--ntasks-per-node=16 "
            "--cpus-per-task=4 "
            "--time=8-0 "
            "--mem-per-cpu=4G "
            "--output=%x_%j.out "
            "--error=%x_%j.err"
        ),
    ]
    return args


# function to run snakemake with settings and add to run folder
def run_snakemake(
    cores, mode, out_dir, run, roadies_dir, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN, gpu, cluster
):

    # Set threads per instance dynamically
    num_threads = cores // fixed_parallel_instances

    cmd = ["snakemake"]
    if cluster:
        cmd += cluster_snakemake_args(mode)
    else:
        cmd += ["--cores", str(cores)]
    cmd += [
        "--config",
        "mode=" + str(mode),
        "config_path=" + str(config_path),
        "num_threads=" + str(num_threads),
        "deep_mode=" + str(deep_mode),
        "MIN_ALIGN=" + str(MIN_ALIGN),
        "gpu=" + str(gpu),
        "--use-conda",
        "--rerun-incomplete",
        "--conda-frontend", "conda"
    ]
    for i in range(len(cmd)):
        if i == len(cmd) - 1:
            print(cmd[i])
        else:
            print(cmd[i], end=" ")
    subprocess.run(cmd)
    # get the run output in folder
    os.system(
        "./workflow/scripts/get_run.sh {0} {1} {2}".format(out_dir, run, roadies_dir)
    )


# function to combine gene trees and mapping files from all iterations
def combine_iter(out_dir, iteration, cores, roadies_dir):
    run = format_run(iteration)

    # Rebuild the cumulative master files from scratch out of every
    # completed iteration's own saved output, rather than incrementally
    # appending onto them. An append is not safe to retry: if a previous
    # attempt at this same iteration crashed after appending but before
    # {run}.nwk was written, retrying with `cat >>` would double-count that
    # iteration's gene trees. Rebuilding from 0..iteration every time is
    # idempotent by construction - it doesn't matter how many times or in
    # what state this gets re-run, the result only depends on which
    # iteration directories exist on disk. Written to temp paths and
    # rename()'d into place so a crash mid-write never leaves a partial
    # master file that a later run might read as complete.
    master_gt_tmp = out_dir + "/master_gt.nwk.tmp"
    master_map_tmp = out_dir + "/master_map.txt.tmp"
    with open(master_gt_tmp, "w") as gt_out, open(master_map_tmp, "w") as map_out:
        for i in range(iteration + 1):
            run_i = format_run(i)
            with open(f"{out_dir}/{run_i}/gene_tree_merged.nwk") as f:
                gt_out.write(f.read())
            with open(f"{out_dir}/{run_i}/mapping.txt") as f:
                map_out.write(f.read())
    os.replace(master_gt_tmp, out_dir + "/master_gt.nwk")
    os.replace(master_map_tmp, out_dir + "/master_map.txt")

    # Same reasoning for the ASTRAL outputs: write to temp paths, only
    # promote to the final {run}.nwk name (the file resume uses as proof
    # this iteration is done) once both calls have actually succeeded.
    nwk_tmp = f"{out_dir}/{run}.nwk.tmp"
    stats_tmp = f"{out_dir}/{run}_stats.nwk.tmp"
    ret1 = os.system(
        "astral-pro3 -t {0} -i {1}/master_gt.nwk -o {2} -a {1}/master_map.txt".format(
            cores, out_dir, nwk_tmp
        )
    )
    ret2 = os.system(
        "astral-pro3 -t {0} -u 3 -i {1}/master_gt.nwk -o {2} -a {1}/master_map.txt".format(
            cores, out_dir, stats_tmp
        )
    )
    if ret1 != 0 or ret2 != 0:
        raise RuntimeError(
            f"astral-pro3 failed for {run} (exit status {ret1}, {ret2}) - "
            f"{run}.nwk was not written, so this iteration will be retried on resume."
        )
    os.replace(nwk_tmp, f"{out_dir}/{run}.nwk")
    os.replace(stats_tmp, f"{out_dir}/{run}_stats.nwk")

    os.system("cp {0}/{1}.nwk {2}/roadies.nwk".format(out_dir, run, roadies_dir))
    os.system("cp {0}/{1}_stats.nwk {2}/roadies_stats.nwk".format(out_dir, run, roadies_dir))
    # open both master files and get gene trees and mapping
    gt = open(out_dir + "/master_gt.nwk", "r")
    gene_trees = gt.readlines()
    gt.close()
    return gene_trees


# function for convergence run
def converge_run(
    iteration,
    cores,
    mode,
    out_dir,
    ref_exist,
    ref,
    roadies_dir,
    support_thr,
    config_path,
    fixed_parallel_instances,
    deep_mode,
    MIN_ALIGN,
    ref_path,
    gpu,
    grow,
    cluster
):
    # Per-iteration scratch space - always wiped fresh regardless of --clean,
    # since each iteration's Snakemake run needs a clean DAG/working
    # directory. This is unrelated to whether the ALL_OUT_DIR convergence
    # history (out_dir) gets preserved across a resume.
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))
    os.system("rm {0}".format('sampling_output.txt'))
    run = format_run(iteration)
    # run snakemake with specificed gene number and length
    if iteration >= 2:
        base_gene_count = read_initial_gene_count(
            config_path
        )  # Read initial GENE_COUNT value
        update_config(config_path, base_gene_count)
    run_snakemake(
        cores, mode, out_dir, run, roadies_dir, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN, gpu, cluster
    )
    # merging gene trees and mapping files
    gene_trees = combine_iter(out_dir, iteration, cores, roadies_dir)
    t = Tree(out_dir + "/" + run + ".nwk")
    # add species tree to tree list
    if ref_exist:
        reroottree = t
        rerootTree(ref, reroottree)
        # print(t)
        reroottree.write(outfile=out_dir + "/" + run + ".rerooted.nwk")
    # extract percentage of gene trees with support value more than from freqQuad.csv
    local_pp_values = []
    count = 0
    with open("freqQuad.csv", "r") as file:
        csv_reader = csv.reader(file, delimiter="\t")
        rows = list(csv_reader)  # Read all rows into a list
        total_rows = len(rows)  # Calculate the total number of rows
        for i, row in enumerate(rows):
            if (i + 1) % 3 == 1:
                value = float(row[3])
                if value >= support_thr:
                    count += 1
                local_pp_values.append(value)

    percent_high_support = (count / (total_rows / 3)) * 100

    # preserve this iteration's quartet-support data before the next
    # iteration's ASTRAL-Pro3 run overwrites the shared freqQuad.csv
    shutil.copy("freqQuad.csv", out_dir + "/" + run + "/freqQuad.csv")

    return percent_high_support, len(gene_trees), t


# main function
if __name__ == "__main__":
    # taking in arguments, have default values for most; information in README.md
    parser = argparse.ArgumentParser(
        prog="Converge",
        description="Script to continuously run snakemake with a small number of genes combining the gene trees after each run",
    )

    parser.add_argument("--cores", type=int, default=32, help="number of cores")
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Config file containing global variables",
    )
    parser.add_argument(
        "--mode",
        default="accurate",
        help="select modes of operations (fast, accurate, balanced)",
    )
    parser.add_argument(
        "--deep",
        action="store_true",
        help="specify if ROADIES will run in deep mode - to capture deeper phylogenetic timescales",
    )
    parser.add_argument(
        "--gpu",
        default="0",
        help="specify number of GPU cores",
    )
    parser.add_argument(
        "--grow",
        action="store_true",
        help="specify if you want to update your tree or grow your tree in placement mode",
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="submit Snakemake rule jobs to SLURM via sbatch for multi-node execution",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="delete the ALL_OUT_DIR convergence directory before running, for a "
             "genuine fresh start (default is to resume from the last iteration "
             "that fully completed - master_gt.nwk/master_map.txt get rebuilt "
             "from just the completed iterations' own saved output each time, so "
             "a crash mid-iteration can't leave duplicated/corrupted data behind)",
    )
    # assigning argument values to variables
    args = vars(parser.parse_args())
    config_path = args["config"]
    CORES = args["cores"]
    MODE = args["mode"]
    deep_mode = args["deep"]
    gpu = args["gpu"]
    grow = args["grow"]
    cluster = args["cluster"]
    clean = args["clean"]
    # read config.yaml for variables
    config = yaml.safe_load(Path(config_path).read_text())
    ref_exist = False
    ref = None
    if config["REFERENCE"] != None:
        ref_exist = True
        ref = Tree(config["REFERENCE"])
    genomes = config["GENOMES"]
    out_dir = config["ALL_OUT_DIR"]
    NUM_GENOMES = len(os.listdir(genomes))
    NUM_GENES = config["GENE_COUNT"]
    LENGTH = config["LENGTH"]
    MIN_ALIGN = max(4, math.ceil(0.1 * NUM_GENOMES))
    support_thr = config["SUPPORT_THRESHOLD"]
    roadies_dir = config["OUT_DIR"]
    fixed_parallel_instances = config["NUM_INSTANCES"]
    ref_path = config["REF_DIR"]
    if clean:
        os.system("rm -r {0}".format(out_dir))
        os.system("mkdir -p " + out_dir)
        iteration = 0
        high_support_list = []
    else:
        os.makedirs(out_dir, exist_ok=True)
        iteration, high_support_list = find_resume_point(out_dir)
        if iteration > 0:
            print(
                f"Resuming from iteration {iteration} "
                f"({len(high_support_list)} already-completed iteration(s) "
                f"found in {out_dir})"
            )
    sys.setrecursionlimit(2000)
    os.system("snakemake --unlock")
    if ref_exist:
        ref_dists = []
    start_time = time.time()
    start_time_l = time.asctime(time.localtime(time.time()))
    with open(out_dir + "/time_stamps.csv", "a") as t_out:
        t_out.write("Start time: " + str(start_time_l) + "\n")
    while True:
        percent_high_support, num_gt, outputtree = converge_run(
            iteration,
            CORES,
            MODE,
            out_dir,
            ref_exist,
            ref,
            roadies_dir,
            support_thr,
            config_path,
            fixed_parallel_instances,
            deep_mode,
            MIN_ALIGN,
            ref_path,
            gpu,
            grow,
            cluster
        )
        curr_time = time.time()
        curr_time_l = time.asctime(time.localtime(time.time()))
        high_support_list.append(percent_high_support)
        elapsed_time = curr_time - start_time
        with open(out_dir + "/time_stamps.csv", "a") as t_out:
            t_out.write(
                str(iteration)
                + ","
                + str(num_gt)
                + ","
                + str(percent_high_support)
                + ","
                + str(curr_time_l)
                + ","
                + str(elapsed_time)
                + "\n"
            )
        # if reference exists get distance between ref and roadies tree
        if ref_exist:
            ref_dist = comp_tree(ref, outputtree)
            ref_dists.append(ref_dist)
            with open(out_dir + "/ref_dist.csv", "a") as ref_out:
                ref_out.write(
                    str(iteration) + "," + str(num_gt) + "," + str(ref_dist) + "\n"
                )

        iteration += 1
        if ((iteration == 1) and (percent_high_support == 100)) or (
            (iteration >= 2)
            and (
                (abs(percent_high_support - high_support_list[iteration - 2]) < 1)
                or (percent_high_support == 100)
                or (iteration == 9)
            )
        ):
            break
