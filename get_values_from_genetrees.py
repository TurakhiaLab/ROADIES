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

# main function
if __name__ == "__main__":
    ref_exist = True
    if ref_exist:
        ref_dists = []
    # Find all newick files in the directory
    newick_files = [
        f for f in os.listdir(".") if f.endswith(".nwk") and f.startswith("gene_trees_")
    ]

    # Sort files by their numerical ID
    newick_files.sort(key=lambda x: int(x.split("_")[2].split(".")[0]))

    config_path = "config/config.yaml"

    config = yaml.safe_load(Path(config_path).read_text())

    iteration = 0

    output_dir = "output"

    os.system("mkdir {0}".format(output_dir))

    for file in newick_files:
        file_id = file.split("_")[2].split(".")[0]

        # open both files and get lines, each line is a separate gene tree
        os.system(
            "ASTER-Linux/bin/astral-pro -t 16 -i {0} -o output/run_{1}.nwk -a mapping.txt".format(
                file, file_id
            )
        )
        os.system(
            "ASTER-Linux/bin/astral-pro -t 16 -u 3 -i {0} -o output/run_{1}_stats.nwk -a mapping.txt".format(
                file, file_id
            )
        )

        # open both master files and get gene trees and mapping
        gt = open(file, "r")
        gene_trees = gt.readlines()
        gt.close()

        t = Tree("output/run_" + file_id + ".nwk")
        # add species tree to tree list
        if ref_exist:
            reroottree = t
            ref = Tree(config["REFERENCE"])
            rerootTree(ref, reroottree)
            # print(t)
            reroottree.write(outfile="output/run_" + file_id + ".rerooted.nwk")

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
                    if (value) >= 0.7:
                        count += 1
                    local_pp_values.append(value)

        percent_high_support = (count / (total_rows / 3)) * 100

        if ref_exist:
            ref_dist = comp_tree(ref, t)
            ref_dists.append(ref_dist)
            with open("output/ref_dist.csv", "a") as ref_out:
                ref_out.write(
                    str(iteration)
                    + ","
                    + str(len(gene_trees))
                    + ","
                    + str(ref_dist)
                    + ","
                    + str(percent_high_support)
                    + "\n"
                )

        iteration += 1
