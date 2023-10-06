# converge.py is a script that iteratively runs ROADIES, wherein after each run, the resultant gene trees are concatenated into a master file, bootstrapped, and input into ASTRAL-PRO.
# These bootstrapped trees are then compared with the previous run (iter_dist), and also with the ref if given(ref_dist)
# The program stops after either running converge for ‘MAX_ITER’ or satisfying the distance threshold ‘t’ for ‘STOP_ITER’ consecutive runs for iter_dists_bs

# REQUIREMENTS: Activated conda environment with snakemake and ete3
# USAGE from wga-phylo directory: `python workflow/scripts/converge.py -c {# of cores} --out_dir {converge output directory} --config {config file}`
import os, sys, glob
import argparse
import random
import subprocess
import signal
from ete3 import Tree
from reroot import rerootTree
import yaml
from pathlib import Path
import multiprocessing
from multiprocessing import Pool
import time
import subprocess
import math


# function that finds the average distance between an array of trees and itself
def comp_tree(t1, t2):
    d = t1.compare(t2)
    return d["norm_rf"]


# function to run snakemake with settings and add to run folder
def run_snakemake(cores, mode):
    cmd = [
        "snakemake",
        "--core",
        str(cores),
        "--jobs",
        str(cores),
        "--config",
        "mode="+str(mode),
        "--use-conda",
        "--rerun-incomplete",
    ]
    for i in range(len(cmd)):
        if i == len(cmd) - 1:
            print(cmd[i])
        else:
            print(cmd[i], end=" ")
    subprocess.run(cmd)
    # get the run output in folder
    print("Adding run to converge folder")
    os.system("./workflow/scripts/get_run.sh")


# function for convergence run
def converge_run(cores, mode, ref_exist, trees, roadies_dir):
    # run snakemake with specificed gene number and length
    run_snakemake(cores, mode)
    # merging gene trees and mapping files
    gt = open(roadies_dir + "/genetrees/gene_tree_merged.nwk", "r")
    gene_trees = gt.readlines()
    gt.close()
    t = Tree(roadies_dir + "/roadies.nwk")
    # add species tree to tree list
    trees.append(t)
    if ref_exist:
        print("Rerooting to reference")
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

    parser.add_argument("-c", type=int, default=16, help="number of cores")
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Config file containing global variables",
    )
    parser.add_argument("--mode", default="accurate", help="select modes of operations (fast, accurate, balanced)")
    # assigning argument values to variables
    args = vars(parser.parse_args())
    config_path = args["config"]
    CORES = args["c"]
    MODE = args["mode"]
    # read config.yaml for variables
    config = yaml.safe_load(Path(config_path).read_text())
    ref_exist = False
    if config["REFERENCE"] != None:
        print("Converge has read reference tree {0}".format(config["REFERENCE"]))
        ref_exist = True
        ref = Tree(config["REFERENCE"])
    genomes = config["GENOMES"]
    NUM_GENOMES = len(os.listdir(genomes))
    NUM_GENES = config["GENE_COUNT"]
    LENGTH = config["LENGTH"]
    roadies_dir = config["OUT_DIR"]
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))
    sys.setrecursionlimit(2000)
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
    num_gt = converge_run(CORES, MODE, ref_exist, trees, roadies_dir)
    print("There are {0} gene trees".format(num_gt))
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
