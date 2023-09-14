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
def run_snakemake(cores, out_dir, run, runtime_left):
    cmd = [
        "snakemake",
        "--core",
        str(cores),
        "--jobs",
        str(cores),
        "--use-conda",
        "--rerun-incomplete",
    ]
    for i in range(len(cmd)):
        if i == len(cmd) - 1:
            print(cmd[i])
        else:
            print(cmd[i], end=" ")
    if runtime_left == math.inf:
        print("Running without time limit")
        subprocess.run(cmd)
    else:
        print("Trying to fit in run with time remaining", runtime_left)
        try:
            subprocess.run(cmd, timeout=runtime_left)
        except:
            print("TIME LIMIT EXCEEDED")
            sys.exit(0)
    # get the run output in folder
    print("Adding run to converge folder")
    os.system("./workflow/scripts/get_run.sh {0} {1}".format(out_dir, run))


# function to combine gene trees and mapping files from all iterations
def combine_iter(out_dir, run):
    print("Concatenating run's gene trees and mapping files with master versions")
    os.system(
        "cat {0}/{1}/gene_tree_merged.nwk >> {0}/master_gt.nwk".format(out_dir, run)
    )
    os.system("cp {0}/master_gt.nwk {0}/{1}.gt.nwk")
    os.system("cat {0}/{1}/mapping.txt >> {0}/master_map.txt".format(out_dir, run))
    os.system("cp {0}/master_map.txt {0}/{1}.map.txt")

    # open both files and get lines, each line is a separate gene tree
    os.system(
        "ASTER-Linux/bin/astral-pro -i {0}/master_gt.nwk -o {0}/{1}.nwk -a {0}/master_map.txt".format(
            out_dir, run
        )
    )
    # open both master files and get gene trees and mapping
    gt = open(out_dir + "/master_gt.nwk", "r")
    gene_trees = gt.readlines()
    gt.close()
    return gene_trees


# function for convergence run
def converge_run(
    iteration, cores, out_dir, ref_exist, trees, roadies_dir, runtime_left
):
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))
    run = "run_"
    # allows sorting runs correctly
    if iteration < 10:
        run += "0" + str(iteration)
    else:
        run += str(iteration)
    print("Starting " + run)
    # run snakemake with specificed gene number and length
    run_snakemake(cores, out_dir, run, runtime_left)
    # merging gene trees and mapping files
    gene_trees = combine_iter(out_dir, run)
    t = Tree(out_dir + "/" + run + ".nwk")
    # add species tree to tree list
    trees.append(t)
    if ref_exist:
        print("Rerooting to reference")
        rerootTree(ref, t)
        # print(t)
        t.write(outfile=out_dir + "/" + run + ".rerooted.nwk")
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
    parser.add_argument("--out_dir", default="converge", help="output dir")
    parser.add_argument(
        "--config",
        default="config/config.yaml",
        help="Config file containing global variables",
    )
    # assigning argument values to variables
    args = vars(parser.parse_args())
    config_path = args["config"]
    CORES = args["c"]
    out_dir = args["out_dir"]
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
    ITERATIONS = config["ITERATIONS"]
    roadies_dir = config["OUT_DIR"]
    master_gt = out_dir + "/master_gt.nwk"
    master_map = out_dir + "/master_map.txt"
    os.system("mkdir -p " + out_dir)
    os.system("touch {0}".format(master_gt))
    os.system("touch {0}".format(master_map))
    os.system("mkdir -p " + out_dir + "/tmp")
    sys.setrecursionlimit(2000)
    # initialize lists for runs and distances
    window_dists = []
    time_stamps = []
    if ref_exist:
        ref_dists = []
    # list of roadies trees after each iteration
    trees = []
    start_time = time.time()
    start_time_l = time.asctime(time.localtime(time.time()))
    time_stamps.append(start_time)
    runtime_left = math.inf
    gt_counts = []
    with open(out_dir + "/time_stamps.csv", "a") as t_out:
        t_out.write("Start time: " + str(start_time_l) + "\n")
    for iteration in range(ITERATIONS):
        # returns an array of b bootstrapped trees
        num_gt = converge_run(
            iteration, CORES, out_dir, ref_exist, trees, roadies_dir, runtime_left
        )
        print("There are {0} gene trees after iteration {1}".format(num_gt, iteration))
        gt_counts.append(num_gt)
        curr_time = time.time()
        curr_time_l = time.asctime(time.localtime(time.time()))
        to_previous = curr_time - time_stamps[len(time_stamps) - 1]
        time_stamps.append(curr_time)
        elapsed_time = curr_time - start_time
        with open(out_dir + "/time_stamps.csv", "a") as t_out:
            t_out.write(
                str(num_gt)
                + ","
                + str(curr_time_l)
                + ","
                + str(elapsed_time)
                + ","
                + str(to_previous)
                + "\n"
            )
        # if reference exists get distance between ref and roadies tree
        if ref_exist:
            ref_dist = comp_tree(ref, trees[iteration])
            ref_dists.append(ref_dist)
            with open(out_dir + "/ref_dist.csv", "a") as ref_out:
                ref_out.write(str(num_gt) + "," + str(ref_dist) + "\n")
