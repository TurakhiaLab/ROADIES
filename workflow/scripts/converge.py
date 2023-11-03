# converge.py is a script that iteratively runs ROADIES, wherein after each run, the resultant gene trees are concatenated into a master file, and input into ASTRAL-PRO.
# The program stops after configured number of ITERATIONS

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
import csv


# function that finds the average distance between an array of trees and itself
def comp_tree(t1, t2):
    d = t1.compare(t2)
    return d["norm_rf"]


# function to run snakemake with settings and add to run folder
def run_snakemake(cores, mode, out_dir, run):
    cmd = [
        "snakemake",
        "--core",
        str(cores),
        "--jobs",
        str(cores),
        "--config",
        "mode=" + str(mode),
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
    os.system("./workflow/scripts/get_run.sh {0} {1}".format(out_dir, run))


# function to combine gene trees and mapping files from all iterations
def combine_iter(out_dir, run):
    os.system(
        "cat {0}/{1}/gene_tree_merged.nwk >> {0}/master_gt.nwk".format(out_dir, run)
    )
    os.system("cp {0}/master_gt.nwk {0}/{1}.gt.nwk")
    os.system("cat {0}/{1}/mapping.txt >> {0}/master_map.txt".format(out_dir, run))
    os.system("cp {0}/master_map.txt {0}/{1}.map.txt")

    # open both files and get lines, each line is a separate gene tree
    os.system(
        "ASTER-Linux/bin/astral-pro -t 16 -i {0}/master_gt.nwk -o {0}/{1}.nwk -a {0}/master_map.txt".format(
            out_dir, run
        )
    )
    os.system(
        "ASTER-Linux/bin/astral-pro -t 16 -u 3 -i {0}/master_gt.nwk -o {0}/{1}_stats.nwk -a {0}/master_map.txt".format(
            out_dir, run
        )
    )
    # open both master files and get gene trees and mapping
    gt = open(out_dir + "/master_gt.nwk", "r")
    gene_trees = gt.readlines()
    gt.close()
    return gene_trees


# function for convergence run
def converge_run(iteration, cores, mode, out_dir, ref_exist, roadies_dir):
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))
    run = "iteration_"
    # allows sorting runs correctly
    if iteration < 10:
        run += "0" + str(iteration)
    else:
        run += str(iteration)
    # run snakemake with specificed gene number and length
    run_snakemake(cores, mode, out_dir, run)
    # merging gene trees and mapping files
    gene_trees = combine_iter(out_dir, run)
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
                if value > 0.8:
                    count += 1
                local_pp_values.append(value)

    percent_high_support = (count / (total_rows / 3)) * 100

    return percent_high_support, len(gene_trees), t


# main function
if __name__ == "__main__":
    # taking in arguments, have default values for most; information in README.md
    parser = argparse.ArgumentParser(
        prog="Converge",
        description="Script to continuously run snakemake with a small number of genes combining the gene trees after each run",
    )

    parser.add_argument("--cores", type=int, default=32, help="number of cores")
    parser.add_argument("--out_dir", default="converge", help="converge directory")
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
    # assigning argument values to variables
    args = vars(parser.parse_args())
    config_path = args["config"]
    CORES = args["cores"]
    MODE = args["mode"]
    out_dir = args["out_dir"]
    # read config.yaml for variables
    config = yaml.safe_load(Path(config_path).read_text())
    ref_exist = False
    if config["REFERENCE"] != None:
        ref_exist = True
        ref = Tree(config["REFERENCE"])
    genomes = config["GENOMES"]
    num_itrs = config["ITERATIONS"]
    NUM_GENOMES = len(os.listdir(genomes))
    NUM_GENES = config["GENE_COUNT"]
    LENGTH = config["LENGTH"]
    roadies_dir = config["OUT_DIR"]
    master_gt = out_dir + "/master_gt.nwk"
    master_map = out_dir + "/master_map.txt"
    os.system("mkdir -p " + out_dir)
    os.system("touch {0}".format(master_gt))
    os.system("touch {0}".format(master_map))
    os.system("mkdir -p " + out_dir + "/tmp")
    sys.setrecursionlimit(2000)
    # initialize lists for runs and distances
    time_stamps = []
    if ref_exist:
        ref_dists = []
    high_support_list = []
    iteration = 0
    start_time = time.time()
    start_time_l = time.asctime(time.localtime(time.time()))
    time_stamps.append(start_time)
    with open(out_dir + "/time_stamps.csv", "a") as t_out:
        t_out.write("Start time: " + str(start_time_l) + "\n")
    while True:
        percent_high_support, num_gt, outputtree = converge_run(
            iteration, CORES, MODE, out_dir, ref_exist, roadies_dir
        )
        curr_time = time.time()
        curr_time_l = time.asctime(time.localtime(time.time()))
        to_previous = curr_time - time_stamps[len(time_stamps) - 1]
        time_stamps.append(curr_time)
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

        high_support_list.append(percent_high_support)

        iteration += 1
        if iteration == num_itrs:
            break
