# noconverge.py is a script that runs ROADIES once and input gene trees into ASTRAL-PRO.

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


# function that finds the average distance between an array of trees and itself
def comp_tree(t1, t2):
    d = t1.compare(t2)
    return d["norm_rf"]


# Rules that spawn one Snakemake job per sampled locus (up to GENE_COUNT
# jobs), keyed by mode. Under --cluster, rules sharing a group_name are
# bundled into one sbatch submission of group_size jobs instead of one
# submission per job.
#
# Keep group_size small relative to --jobs: --group-components re-chunks
# ALL jobs assigned to a group id on every new scheduling wave, including
# ones already merged by a previous wave, so a group can silently balloon
# past its declared size (and its sbatch resource request) as more jobs are
# discovered. lastz_batch/pasta here comfortably exceed --jobs on their own
# and are left ungrouped for that reason - each job just requests its own
# declared threads/resources.
GROUPED_RULES = {
    "accurate": [
        ("pasta", "group0", 250),
        ("filtermsa", "group0", 250),
        ("raxmlng", "group0", 250),
    ],
    "balanced": [
        ("pasta", "group0", 250),
        ("filtermsa", "group0", 250),
        ("fasttree", "group0", 250),
    ],
}

# Multi-node SLURM execution args, opt-in via --cluster. The submit command is
# templated with {threads}/{resources.mem_mb} so each rule gets sized for what
# it actually needs, instead of every job requesting a flat hardcoded amount.
# --jobs/--cores/--resources below are the *ceiling* Snakemake uses to
# bin-pack grouped jobs - actual per-job requests still come from each rule's
# own threads/resources.
def cluster_snakemake_args(mode, partition, account, time_limit):
    args = [
        "--jobs", "200",
        "--cores", "64",
        "--resources", "mem_mb=256000",
        "--latency-wait", "120",
        "--keep-going",
        # greedy instead of Snakemake's default ilp scheduler: placement's
        # pasta stage has ~100k+ mutually independent jobs ready at once,
        # and ilp's per-round optimization over that many candidates became
        # the wall-clock bottleneck. greedy has nothing to optimize here
        # anyway, since every pasta job wants the same threads/mem_mb.
        "--scheduler", "greedy",
    ]
    grouped_rules = GROUPED_RULES.get(mode, [])
    if grouped_rules:
        args += ["--groups"] + [f"{rule}={group}" for rule, group, _ in grouped_rules]
        seen_groups = {}
        for _, group, size in grouped_rules:
            seen_groups[group] = size
        args += ["--group-components"] + [
            f"{group}={size}" for group, size in seen_groups.items()
        ]
    args += [
        "--executor",
        "cluster-generic",
        "--cluster-generic-submit-cmd",
        (
            "sbatch "
            "--job-name=ROADIES_run "
            f"--partition={partition} "
            f"--account={account} "
            "--nodes=1 "
            "--ntasks=1 "
            "--cpus-per-task={threads} "
            "--mem={resources.mem_mb}M "
            f"--time={time_limit} "
            "--output=%x_%j.out "
            "--error=%x_%j.err"
        ),
    ]
    return args


# function to run snakemake with settings and add to run folder
def run_snakemake(cores, mode, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN, gpu, cluster,
                   cluster_partition, cluster_account, cluster_time):

    # Set threads per instance dynamically
    num_threads = cores // fixed_parallel_instances

    cmd = ["snakemake"]
    if cluster:
        cmd += cluster_snakemake_args(mode, cluster_partition, cluster_account, cluster_time)
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
        "cluster=" + str(cluster),
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
    os.system("./workflow/scripts/get_run_noconverge.sh")


# function for convergence run
def converge_run(
    cores,
    mode,
    ref_exist,
    ref,
    trees,
    roadies_dir,
    config_path,
    fixed_parallel_instances,
    deep_mode,
    MIN_ALIGN,
    ref_path,
    gpu,
    grow,
    cluster,
    cluster_partition,
    cluster_account,
    cluster_time
):
    # run snakemake with specificed gene number and length
    run_snakemake(cores, mode, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN, gpu, cluster,
                   cluster_partition, cluster_account, cluster_time)
    # ASTRAL-Pro3 segfaults on 0 input gene trees rather than erroring cleanly,
    # which otherwise cascades into a confusing downstream ete3 NewickError.
    # Most often means every sampled locus failed MIN_ALIGN - too few
    # genomes/query sequences, or IDENTITY/COVERAGE too strict for this data.
    merged_gt_path = roadies_dir + "/genetrees/gene_tree_merged.nwk"
    if not os.path.exists(merged_gt_path) or os.path.getsize(merged_gt_path) == 0:
        raise RuntimeError(
            f"No gene trees were produced ({merged_gt_path} is missing or empty) - "
            "every sampled locus likely failed the MIN_ALIGN species-count filter. "
            "Try more genomes, a higher GENE_COUNT, or relaxing IDENTITY/COVERAGE in config.yaml."
        )
    if (mode == 'placement'):
        os.system(
            "cat {0}/genes/mapping.txt {1}/genes/mapping.txt >> {0}/genes/mapping_combined.txt".format(
                roadies_dir, ref_path
            )
        )
        if (grow):
            os.system(
                "astral-pro3 -t {1} --constraint {2}/roadies.nwk -i {0}/genetrees/gene_tree_merged.nwk -o {0}/roadies.nwk -a {0}/genes/mapping_combined.txt".format(
                    roadies_dir, cores, ref_path
                )
            )
            os.system(
                "astral-pro3 -t {1} -u 3 --constraint {2}/roadies.nwk -i {0}/genetrees/gene_tree_merged.nwk -o {0}/roadies_stats.nwk -a {0}/genes/mapping_combined.txt".format(
                    roadies_dir, cores, ref_path
                )
            )
        else:
            os.system(
                "astral-pro3 -t {1} -i {0}/genetrees/gene_tree_merged.nwk -o {0}/roadies.nwk -a {0}/genes/mapping_combined.txt".format(
                    roadies_dir, cores
                )
            )
            os.system(
                "astral-pro3 -t {1} -u 3 -i {0}/genetrees/gene_tree_merged.nwk -o {0}/roadies_stats.nwk -a {0}/genes/mapping_combined.txt".format(
                    roadies_dir, cores
                )
            )

    else:
        os.system(
            "astral-pro3 -t {1} -i {0}/genetrees/gene_tree_merged.nwk -o {0}/roadies.nwk -a {0}/genes/mapping.txt".format(
                roadies_dir, cores
            )
        )
        os.system(
            "astral-pro3 -t {1} -u 3 -i {0}/genetrees/gene_tree_merged.nwk -o {0}/roadies_stats.nwk -a {0}/genes/mapping.txt".format(
                roadies_dir, cores
            )
        )

    # keep this run's quartet-support data with its own output instead of
    # leaving it in the ROADIES install directory, where the next run would overwrite it
    if os.path.exists("freqQuad.csv"):
        shutil.copy("freqQuad.csv", roadies_dir + "/freqQuad.csv")

    gt = open(roadies_dir + "/genetrees/gene_tree_merged.nwk", "r")
    gene_trees = gt.readlines()
    gt.close()
    t = Tree(roadies_dir + "/roadies.nwk")
    # add species tree to tree list
    trees.append(t)
    if ref_exist:
        rerootTree(ref, t)
        # print(t)
        t.write(outfile=roadies_dir + "/roadies_rerooted.nwk")
    # create bootstrapping trees
    return len(gene_trees)


# main function
if __name__ == "__main__":
    # taking in arguments, have default values for most; information in README.md
    parser = argparse.ArgumentParser(
        prog="Noconverge",
        description="Script to run snakemake for one iteration",
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
        "--clean",
        action="store_true",
        help="delete the output directory before running, for a genuine fresh start "
             "(default is to leave existing output alone and let Snakemake's "
             "--rerun-incomplete resume it - safer against accidental double-launches)",
    )
    parser.add_argument(
        "--cluster",
        action="store_true",
        help="submit Snakemake rule jobs to SLURM via sbatch for multi-node execution",
    )
    # assigning argument values to variables
    args = vars(parser.parse_args())
    config_path = args["config"]
    CORES = args["cores"]
    MODE = args["mode"]
    deep_mode = args["deep"]
    gpu = args["gpu"]
    grow = args["grow"]
    clean = args["clean"]
    cluster = args["cluster"]
    # read config.yaml for variables
    config = yaml.safe_load(Path(config_path).read_text())
    ref_exist = False
    ref = None
    if config["REFERENCE"] != None:
        ref_exist = True
        ref = Tree(config["REFERENCE"])
    genomes = config["GENOMES"]
    NUM_GENOMES = len(os.listdir(genomes))
    NUM_GENES = config["GENE_COUNT"]
    LENGTH = config["LENGTH"]
    MIN_ALIGN = max(4, math.ceil(0.1 * NUM_GENOMES))
    roadies_dir = config["OUT_DIR"]
    fixed_parallel_instances = config["NUM_INSTANCES"]
    ref_path = config.get("REF_DIR")
    cluster_partition = config.get("CLUSTER_PARTITION", "long")
    cluster_account = config.get("CLUSTER_ACCOUNT", "standard")
    cluster_time = config.get("CLUSTER_TIME", "8-0")
    if clean:
        os.system("rm -r {0}".format(roadies_dir))
        os.system("mkdir {0}".format(roadies_dir))
        os.system("rm {0}".format('sampling_output.txt'))
    else:
        os.makedirs(roadies_dir, exist_ok=True)
    sys.setrecursionlimit(2000)
    os.system("snakemake --unlock --config config_path={0}".format(config_path))
    # initialize lists for runs and distances
    time_stamps = []
    if ref_exist:
        ref_dists = []
    # list of roadies trees after each iteration
    trees = []
    start_time = time.time()
    start_time_l = time.asctime(time.localtime(time.time()))
    time_stamps.append(start_time)
    with open(roadies_dir + "/time_stamps.csv", "a") as t_out:
        t_out.write("Start time: " + str(start_time_l) + "\n")
    num_gt = converge_run(
        CORES,
        MODE,
        ref_exist,
        ref,
        trees,
        roadies_dir,
        config_path,
        fixed_parallel_instances,
        deep_mode,
        MIN_ALIGN,
        ref_path,
        gpu,
        grow,
        cluster,
        cluster_partition,
        cluster_account,
        cluster_time
    )
    curr_time = time.time()
    curr_time_l = time.asctime(time.localtime(time.time()))
    to_previous = curr_time - time_stamps[len(time_stamps) - 1]
    time_stamps.append(curr_time)
    elapsed_time = curr_time - start_time
    with open(roadies_dir + "/time_stamps.csv", "w") as t_out:
        t_out.write(
            str(num_gt) + "," + str(curr_time_l) + "," + str(elapsed_time) + "\n"
        )
    # if reference exists get distance between ref and roadies tree
    if ref_exist:
        ref_dist = comp_tree(ref, trees[0])
        ref_dists.append(ref_dist)
        with open(roadies_dir + "/ref_dist.csv", "a") as ref_out:
            ref_out.write(str(num_gt) + "," + str(ref_dist) + "\n")
    print("Species tree created")
