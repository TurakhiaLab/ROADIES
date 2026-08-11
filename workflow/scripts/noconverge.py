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
# jobs), keyed by mode. Under --cluster these get grouped so they don't each
# turn into a separate sbatch submission. Each entry is (rule, group_name,
# group_size) - rules sharing a group_name get bundled into the same sbatch
# submissions.
#
# --group-components snowballs when a group has far more members than fit in
# one Snakemake scheduling wave (i.e. more than roughly --jobs): the DAG's
# incremental update loop re-runs _update_group_components() on every new
# wave of discovered jobs, and each call re-chunks *all* jobs assigned to
# that group id so far - including ones already merged into a group by a
# PREVIOUS call. A group-of-4 formed in wave 1 becomes one "component" going
# into wave 2's chunking, gets merged with another already-formed group-of-4,
# and keeps compounding across waves. Confirmed live: lastz_batch was
# configured group-of-4 (16 cpu/62.5G requested) but actually executed with
# 16 members (needed ~64cpu/256G), causing the sbatch allocation to OOM-kill
# 6 of 16 concurrent lastz_40 processes. placement's lastz_batch (352 total
# jobs) and pasta (up to 32000) both vastly exceed --jobs 48, so both are
# left ungrouped here - each job requests exactly its own declared
# threads/resources with no bin-packing/merge risk. accurate/balanced modes
# haven't hit this in practice yet but share the identical mechanism (same
# many-disconnected-jobs-under-a-jobs-cap shape) - treat their grouping as
# equally suspect before relying on it for a real run.
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
# it actually needs, instead of every job (lastz included) requesting a flat
# 64 cpus / 256G. --cores/--resources below are the *ceiling* Snakemake uses
# to bin-pack grouped jobs (see pasta's per-locus resources in placement.smk):
# 64 threads / 256000 mem_mb reproduces the old fixed per-group-job request
# (4 concurrent pasta instances @ 16 threads / 64G each), just derived
# explicitly instead of accidentally via a hardcoded string applied to every
# rule. lastz_batch's per-job resources (pair_align_placement_batch.smk) are
# tiny in comparison, so a group of 4 of them fits well inside that same
# ceiling without ever needing to split into multiple layers.
def cluster_snakemake_args(mode):
    args = [
        # Raised 48 -> 200 alongside right-sizing pasta's resources (8cpu/2G
        # for real placements, 1cpu/500M for the touch/cp fallback - was a
        # flat 16cpu/64G for everything). 48 concurrent jobs at the old
        # footprint could already saturate a lot of cluster capacity; at the
        # new footprint it's tiny, so the --jobs cap (not node availability)
        # was going to be the limiting factor. Checked cluster-wide headroom
        # via `sinfo -p long` before picking this: ~3062 idle cpus across
        # the partition at the time (shared with other users, so this is
        # advisory not exclusive - actual concurrency still depends on
        # fairshare/what else is running).
        "--jobs", "200",
        "--cores", "64",
        "--resources", "mem_mb=256000",
        "--latency-wait", "120",
        "--keep-going",
        # Snakemake defaults to --scheduler ilp, which solves an integer
        # linear program over every currently-ready job each round to
        # optimize resource packing. Fine when few jobs are ready at once,
        # but placement's pasta stage has ~100k+ mutually-independent jobs
        # all ready simultaneously (nothing depends on anything else), so
        # ILP has to optimize over a huge candidate pool every round.
        # Confirmed via sstat on a live 128k-locus run: driver CPU time was
        # ~48% of wall-clock with rounds taking 1-2 minutes to select only
        # 15-36 jobs, even though individual pasta jobs (many are instant
        # touch/cp fallbacks) run in seconds. greedy has nothing meaningful
        # to optimize here anyway - every pasta job wants the same
        # threads/mem_mb - so its speed is pure upside for this workload.
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
            "--partition=long "
            "--account=standard "
            "--nodes=1 "
            "--ntasks=1 "
            "--cpus-per-task={threads} "
            "--mem={resources.mem_mb}M "
            "--time=8-0 "
            "--output=%x_%j.out "
            "--error=%x_%j.err"
        ),
    ]
    return args


# function to run snakemake with settings and add to run folder
def run_snakemake(cores, mode, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN, gpu, cluster):

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
    cluster
):
    # run snakemake with specificed gene number and length
    run_snakemake(cores, mode, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN, gpu, cluster)
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
    ref_path = config["REF_DIR"]
    if clean:
        os.system("rm -r {0}".format(roadies_dir))
        os.system("mkdir {0}".format(roadies_dir))
        os.system("rm {0}".format('sampling_output.txt'))
    else:
        os.makedirs(roadies_dir, exist_ok=True)
    sys.setrecursionlimit(2000)
    os.system("snakemake --unlock")
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
        cluster
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
