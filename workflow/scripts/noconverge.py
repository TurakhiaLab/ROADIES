# noconverge.py is a script that runs ROADIES once and input gene trees into ASTRAL-PRO.

# REQUIREMENTS: Activated conda environment with snakemake and ete3
# USAGE: `python workflow/scripts/noconverge.py -c {# of cores} --out_dir {converge output directory} --config {config file}`

import os, sys, glob
import argparse
import random
import subprocess
import signal
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


# function to run snakemake with settings and add to run folder
def run_snakemake(cores, mode, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN):

    # Set threads per instance dynamically
    num_threads = cores // fixed_parallel_instances

    cmd = [
        "snakemake",
        "--cores",
        str(cores),
        "--config",
        "mode=" + str(mode),
        "config_path=" + str(config_path),
        "num_threads=" + str(num_threads),
        "deep_mode=" + str(deep_mode),
        "MIN_ALIGN=" + str(MIN_ALIGN),
        "--use-conda",
        "--rerun-incomplete",
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
):
    # run snakemake with specificed gene number and length
    run_snakemake(cores, mode, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN)
    # merging gene trees and mapping files
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
        default="False",
        help="specify if ROADIES will run in deep mode - to capture deeper phylogenetic timescales",
    )
    # assigning argument values to variables
    args = vars(parser.parse_args())
    config_path = args["config"]
    CORES = args["cores"]
    MODE = args["mode"]
    deep_mode = args["deep"]
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
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))
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
