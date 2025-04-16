# converge.py is a script that iteratively runs ROADIES, wherein after each run, the resultant gene trees are concatenated into a master file, and input into ASTRAL-PRO.
# The program stops after configured number of ITERATIONS

# REQUIREMENTS: Activated conda environment with snakemake and ete3
# USAGE: `python workflow/scripts/converge.py -c {# of cores} --out_dir {converge output directory} --config {config file}`

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


# function to run snakemake with settings and add to run folder
def run_snakemake(
    cores, mode, out_dir, run, roadies_dir, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN
):

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
    # get the run output in folder
    os.system(
        "./workflow/scripts/get_run.sh {0} {1} {2}".format(out_dir, run, roadies_dir)
    )


# function to combine gene trees and mapping files from all iterations
def combine_iter(out_dir, run, cores):
    os.system(
        "cat {0}/{1}/gene_tree_merged.nwk >> {0}/master_gt.nwk".format(out_dir, run)
    )
    os.system("cp {0}/master_gt.nwk {0}/{1}.gt.nwk")
    os.system("cat {0}/{1}/mapping.txt >> {0}/master_map.txt".format(out_dir, run))
    os.system("cp {0}/master_map.txt {0}/{1}.map.txt")

    # open both files and get lines, each line is a separate gene tree
    os.system(
        "astral-pro3 -t {2} -i {0}/master_gt.nwk -o {0}/{1}.nwk -a {0}/master_map.txt".format(
            out_dir, run, cores
        )
    )
    os.system(
        "astral-pro3 -t {2} -u 3 -i {0}/master_gt.nwk -o {0}/{1}_stats.nwk -a {0}/master_map.txt".format(
            out_dir, run, cores
        )
    )
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
):
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))
    run = "iteration_"
    # allows sorting runs correctly
    if iteration < 10:
        run += "0" + str(iteration)
    else:
        run += str(iteration)
    # run snakemake with specificed gene number and length
    if iteration >= 2:
        base_gene_count = read_initial_gene_count(
            config_path
        )  # Read initial GENE_COUNT value
        update_config(config_path, base_gene_count)
    run_snakemake(
        cores, mode, out_dir, run, roadies_dir, config_path, fixed_parallel_instances, deep_mode, MIN_ALIGN
    )
    # merging gene trees and mapping files
    gene_trees = combine_iter(out_dir, run, cores)
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
    out_dir = config["ALL_OUT_DIR"]
    NUM_GENOMES = len(os.listdir(genomes))
    NUM_GENES = config["GENE_COUNT"]
    LENGTH = config["LENGTH"]
    MIN_ALIGN = max(4, math.ceil(0.1 * NUM_GENOMES))
    support_thr = config["SUPPORT_THRESHOLD"]
    roadies_dir = config["OUT_DIR"]
    fixed_parallel_instances = config["NUM_INSTANCES"]
    master_gt = out_dir + "/master_gt.nwk"
    master_map = out_dir + "/master_map.txt"
    os.system("rm -r {0}".format(out_dir))
    os.system("mkdir -p " + out_dir)
    os.system("touch {0}".format(master_gt))
    os.system("touch {0}".format(master_map))
    sys.setrecursionlimit(2000)
    os.system("snakemake --unlock")
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
        )
        curr_time = time.time()
        curr_time_l = time.asctime(time.localtime(time.time()))
        to_previous = curr_time - time_stamps[len(time_stamps) - 1]
        time_stamps.append(curr_time)
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
