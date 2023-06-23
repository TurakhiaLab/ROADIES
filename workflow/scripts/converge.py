# converge.py is a script that iteratively runs ROADIES, wherein after each run, the resultant gene trees are concatenated into a master file, bootstrapped, and input into ASTRAL-PRO.
# These bootstrapped trees are then compared with the previous run (iter_dist), and also with the ref if given(ref_dist)
# The program stops after either running converge for ‘MAX_ITER’ or satisfying the distance threshold ‘t’ for ‘STOP_ITER’ consecutive runs for iter_dists_bs

# REQUIREMENTS: Activated conda environment with snakemake and ete3
# USAGE from wga-phylo directory: `python workflow/scripts/converge.py -c {# of cores} --out_dir {converge output directory} --config {config file}`
import os, sys
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


# function that finds the average distance between an array of trees and itself
def comp_tree(t1, t2):
    d = t1.compare(t2)
    return d["norm_rf"]


# function that finds the average distance between one array of bootstrapped trees with final tree of previous iteration
def comp_bs_final(bs_tree_array, final_tree_prev, idx):
    count = 0
    dist = 0
    for i in range(len(bs_tree_array)):
        count += 1
        print(
            "Comparing run {0} final tree of previous iteration {1} with bootstrap of iteration {2}".format(
                idx, idx - 1, i
            )
        )
        d = final_tree_prev.compare(bs_tree_array[i])
        dist += d["norm_rf"]
    avg = float(dist) / count
    with open(out_dir + "/iter_dist_bs.csv", "a") as w:
        w.write(str(idx) + "," + str(avg) + "\n")
    return avg


class Alarm(Exception):
    pass


def alarm_handler(*args):
    raise Alarm("timeout")


# function to run snakemake with settings and add to run folder
def run_snakemake(
    c, l, k, out_dir, run, roadies_dir, weighted, species_ids, gene_ids, species_lists
):
    if weighted:
        cmd = "snakemake --core {0} --jobs {0} --use-conda --rerun-incomplete --config LENGTH={1} KREG={2} OUTDIR={3}".format(
            c, l, k, roadies_dir
        )
    else:
        cmd = "snakemake --core {0} --jobs {0} --use-conda  --rerun-incomplete --config LENGTH={1} KREG={2} OUTDIR={3} WEIGHTED=0.0 TO_ALIGN=".format(
            c, l, k, roadies_dir
        )
    os.system(cmd)
    # get the run output in folder
    print("Adding run to converge folder")
    os.system(
        "./workflow/scripts/get_run.sh {0} {1} {2} {3} {4}".format(
            out_dir, run, species_ids, gene_ids, species_lists
        )
    )


# function that returns an array of b bootstrapped newick trees
def process_bootstrap(i, out_dir, run, gene_trees):
    print("Creating bootstrapping tree: ", i)
    # path to temp gt
    tmp_path = out_dir + "/tmp/" + run + "." + str(i)
    out = open(tmp_path + ".gt.nwk", "w")
    w = open(tmp_path + ".map.txt", "w")
    leaves = []
    # sample a random line and output to tmp file
    for k in range(len(gene_trees)):
        n = random.randint(0, len(gene_trees) - 1)
        out.write(gene_trees[n])
        # read individual tree in gene trees
        n = Tree(gene_trees[n])
        leaves = n.get_leaf_names()
        for leaf in leaves:
            w.write(leaf + " ")
            s = leaf.split("_")
            for i in range(len(s) - 1):
                if i == len(s) - 2:
                    w.write(s[i] + "\n")
                    # print(s[i])
                else:
                    w.write(s[i] + "_")
                # print(s[i],end='_')
            # print(leaf +' '+ s[0])
    out.close()
    w.close()
    # run astral on bootstrapped tree
    print(
        "ASTER-Linux/bin/astral-pro -i {0}.gt.nwk -o {0}.nwk -a {0}.map.txt".format(
            tmp_path
        )
    )
    os.system(
        "ASTER-Linux/bin/astral-pro -i {0}.gt.nwk -o {0}.nwk -a {0}.map.txt".format(
            tmp_path
        )
    )
    boot_tree = Tree(tmp_path + ".nwk")

    return boot_tree


# function that returns an array of b bootstrapped newick trees
def bootstrap(b, out_dir, run, gene_trees):
    bs_trees = []
    pool = multiprocessing.Pool(b)
    args_list = [(i, out_dir, run, gene_trees) for i in range(b)]
    bs_trees = pool.starmap(process_bootstrap, args_list)
    pool.close()

    return bs_trees


# function to combine gene trees and mapping files from all iterations
def combine_iter(out_dir, run):
    print("Concatenating run's gene trees and mapping files with master versions")
    os.system(
        "cat {0}/{1}/gene_tree_merged.nwk >> {0}/master_gt.nwk".format(out_dir, run)
    )
    os.system("cat {0}/{1}/mapping.txt >> {0}/master_map.txt".format(out_dir, run))
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
    i,
    l,
    k,
    c,
    out_dir,
    b,
    ref_exist,
    trees,
    roadies_dir,
    weighted,
    species_ids,
    gene_ids,
    species_lists,
):
    os.system("rm -r {0}".format(roadies_dir))
    os.system("mkdir {0}".format(roadies_dir))

    run = "run_"
    # allows sorting runs correctly
    if i < 10:
        run += "0" + str(i)
    else:
        run += str(i)
    print("Starting " + run)
    # run snakemake with specificed gene number and length
    run_snakemake(
        c,
        l,
        k,
        out_dir,
        run,
        roadies_dir,
        weighted,
        species_ids,
        gene_ids,
        species_lists,
    )
    # merging gene trees and mapping files
    gene_trees = combine_iter(out_dir, run)
    t = Tree(out_dir + "/" + run + ".nwk")
    # add species tree to tree list
    trees.append(t)
    if ref_exist:
        print("Rerooting to reference")
        rerootTree(ref, t)
        print(t)
        t.write(outfile=out_dir + "/" + run + ".tmp")
        os.system("rm {0}/{1}".format(out_dir, run))
        os.system("mv {0}/{1}.tmp {0}/{1}".format(out_dir, run))
    # create bootstrapping trees
    return bootstrap(b, out_dir, run, gene_trees)


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
    c = args["c"]
    out_dir = args["out_dir"]
    # read config.yaml for variables
    config = yaml.safe_load(Path(config_path).read_text())
    ref_exist = False
    if config["REFERENCE"] != None:
        print("Converge has read reference tree {0}".format(config["REFERENCE"]))
        ref_exist = True
        ref = Tree(config["REFERENCE"])
    k = config["KREG"]
    t = config["DIST_THRESHOLD"]
    l = config["LENGTH"]
    b = config["NUM_BOOTSTRAP"]
    input_gt = config["INPUT_GENE_TREES"]
    input_map = config["INPUT_MAP"]
    max_iter = config["MAX_ITER"]
    stop_iter = config["STOP_ITER"]
    roadies_dir = config["OUT_DIR"]
    species_ids = config["SPECIES_IDS"]
    gene_ids = config["GENE_IDS"]
    species_lists = config["SPECIES_LISTS"]
    quartets = config["QUARTETS"]
    master_gt = out_dir + "/master_gt.nwk"
    master_map = out_dir + "/master_map.txt"
    os.system("mkdir -p " + out_dir)
    if not os.path.isfile(species_ids):
        print("Creating species id file at {0}".format(species_ids))
        os.system("touch {0}".format(species_ids))
    if not os.path.isfile(gene_ids):
        print("Creating gene id file at {0}".format(gene_ids))
        os.system("touch {0}".format(gene_ids))
    if input_gt is None and input_map is None:
        os.system("touch {0}".format(master_gt))
        os.system("touch {0}".format(master_map))
    elif input_gt != master_gt or input_map != master_map:
        print(
            "Reading in reference gene trees {0} and reference mapping {1}".format(
                input_gt, input_map
            )
        )
        os.system("cp {0} {1}".format(input_gt, master_gt))
        os.system("cp {0} {1}".format(input_map, master_map))
    if quartets != "freqQuad.csv":
        print(
            "Input quartet file is different from ./freqQuad.csv Copying quartet file to ./freqQuad.csv"
        )
        if os.path.isfile("freqQuad.csv"):
            os.system("rm freqQuad.csv")
            os.system("cp {0} freqQuad.csv".format(quartets))
        else:
            os.system("cp {0} freqQuad.csv".format(quartets))
        os.system("rm freqQuad.csv")
        os.system("cp {0} freqQuad.csv".format(quartets))
    else:
        if not os.path.isfile("freqQuad.csv"):
            print("Cannot find freqQuad.csv in current directory, so making one ")
            os.system("touch freqQuad.csv")
    os.system("mkdir -p " + out_dir + "/tmp")
    sys.setrecursionlimit(2000)
    # initialize lists for runs and distances
    runs = []
    iter_dists = []
    iter_dists_bs = []
    window_dists = []
    if ref_exist:
        ref_dists = []
    # list of roadies trees after each iteration
    trees = []
    # open files for writing distances
    if ref_exist:
        ref_out = open(out_dir + "/ref_dist.csv", "w")
    iter_out = open(out_dir + "/iter_dist.csv", "w")
    iter_bs = open(out_dir + "/iter_dist_bs.csv", "w")
    avg_out = open(out_dir + "/window_dist.csv", "w")
    # starts main for loop for multiple iterations
    # for max iteration runs; start from 1 index instead of 0
    for i in range(max_iter):
        # returns an array of b bootstrapped trees
        weighted = True
        if i == 0 and input_gt is None:
            weighted = False
        run = converge_run(
            i,
            l,
            k,
            c,
            out_dir,
            b,
            ref_exist,
            trees,
            roadies_dir,
            weighted,
            species_ids,
            gene_ids,
            species_lists,
        )
        runs.append(run)
        # if reference exists get distance between ref and roadies tree
        if ref_exist:
            ref_dist = comp_tree(ref, trees[i])
            ref_dists.append(ref_dist)
            ref_out.write(str(i) + "," + str(ref_dist) + "\n")
        stop_run = False
        iter_flag = False

        # since comparing to previous, hence starts with 1st iteration
        if i >= 1:
            print(len(runs))
            print(runs)
            # get bootstrapped iteration distance
            # iter_dist_bs = comp_runs_bs(runs[i-1],runs[i],i)
            # get avg distance between current bootstrapped trees and previous iter final tree
            print(
                "Getting bootstrapped iteration distances betseen runs {0} and {1}".format(
                    i, i - 1
                )
            )
            iter_dist_bs = comp_bs_final(runs[i], trees[i - 1], i)
            print(
                "Average distance between final tree of iteration {0} and bootstrapped tree of iteration {1} is: {2}".format(
                    i - 1, i, iter_dist_bs
                )
            )
            iter_dists_bs.append(iter_dist_bs)
            iter_bs.write(str(i) + "," + str(iter_dist_bs) + "\n")
            # get single iteration distance
            iter_dist = comp_tree(trees[i], trees[i - 1])
            print(
                "Distance between final trees of iteration {0} and {1} is: {2}".format(
                    i - 1, i, iter_dist
                )
            )
            iter_dists.append(iter_dist)
            iter_out.write(str(i) + "," + str(iter_dist) + "\n")
            # sometimes dist files get updated slowly so create intermediate file for analysis in meantime
            if i % 5 == 0:
                print("Outputting intermediate file")
                iterm = open(out_dir + "/iterm_iter_dist_bs.csv", "w")
                refm = open(out_dir + "/iterm_ref_dist.csv", "w")
                for j in range(len(iter_dists_bs)):
                    iterm.write(str(j) + "," + str(iter_dists_bs[j]) + "\n")
                for j in range(len(ref_dists)):
                    refm.write(str(j) + "," + str(ref_dists[j]) + "\n")
                iterm.close()
                refm.close()
            # checking to see if we should stop script
            # every distance within the window must be lower than threshold to stop
            temp = 0
            temp2 = 0
            # get the difference of average distances between two windows of size stop_iter
            if i > 2 * stop_iter - 1:
                print("Seeing if we should stop")
                for j in range(stop_iter):
                    temp = temp + iter_dists_bs[i - j - 1]
                    temp2 = temp2 + iter_dists_bs[i - j - stop_iter - 1]
                temp = temp / stop_iter
                temp2 = temp2 / stop_iter
                window_diff = abs(temp2 - temp)
                avg_out.write(str(i) + "," + str(window_diff) + "\n")
                print(
                    "The difference between windows of size {0} from run {1} is: {2}".format(
                        stop_iter, i, window_diff
                    )
                )
                window_dists.append(window_diff)
                if window_diff < t:
                    print("crossed threshold")
                    stop_run = True
                else:
                    print("did not cross threshold")
                    stop_run = False
                    iter_flag = True
            print("Distance to previous Tree")
            for j in range(len(iter_dists)):
                print("Run: " + str(j) + ": " + str(iter_dists[j]))
            print(
                "Distances to previous iterations final trees with bootstrapped trees so far:"
            )
            for j in range(len(iter_dists_bs)):
                print("Run: " + str(j) + ": " + str(iter_dists_bs[j]))
            if ref_exist:
                print("Distance to reference so far:")
                for j in range(len(ref_dists)):
                    print("Run: " + str(j) + ": " + str(ref_dists[j]))
            if i > 2 * stop_iter - 1:
                print(
                    "Difference between average of windows of size {0} so far:".format(
                        stop_iter
                    )
                )
                for j in range(len(window_dists)):
                    print(
                        "Run: " + str(j + (2 * stop_iter)) + ": " + str(window_dists[j])
                    )
            if stop_run == True and iter_flag == False:
                break
        if stop_run == True and iter_flag == False:
            break
    if ref_exist:
        ref_out.close()
    iter_out.close()
