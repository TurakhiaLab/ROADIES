# converge.py is a script that iteratively runs ROADIES, wherein after each run, the resultant gene trees are concatenated into a master file, bootstrapped, and input into ASTRAL-PRO.
# These bootstrapped trees are then compared within the same run (self_dist), with the previous run (iter_dist), and also with the ref if given(ref_dist)
# The program stops after either running converge for ‘MAX_ITER’ or satisfying the distance threshold ‘t’ for ‘STOP_ITER’ consecutive runs for both self_dists and iter_dists_bs

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


"""
def comp_self(run,idx):
    dist = 0
    count = 0
    print(run)
    for i in range(len(run)-1):
        for j in range(i+1,len(run)):
            count +=1
            print("Self-comparing run {0} bootstrap trees {1} and {2}".format(idx,i,j))
            d = run[i].compare(run[j])
            dist += d['norm_rf']
    avg = float(dist)/count
    with open(out_dir+'/self_dist.csv','a') as w:
        w.write(str(idx)+','+str(avg)+'\n')
    return avg

#function that finds the average distance between two arrays of trees
def comp_runs_bs(run1,run2,idx):
    count = 0
    dist = 0
    for i in range(len(run1)):
        for j in range(len(run2)):
            count += 1
            print("Comparing run {0} bootstrap {1} to prev bootstrap {2}".format(idx,i,j))
            d = run1[i].compare(run2[j])
            dist += d['norm_rf']
    avg = float(dist)/count
    with open(out_dir+'/iter_dist_bs.csv','a') as w:
        w.write(str(idx)+','+str(avg)+'\n')
    return avg

"""


# function that finds the average distance between one array of bootstrapped trees with final tree of current iteration
def comp_self_final(bs_tree_array, final_tree_cur, idx):
    dist = 0
    count = 0
    for i in range(len(bs_tree_array)):
        count += 1
        print(
            "Self-comparing run {0} final tree of current iteration {1} to bootstrap trees {2}".format(
                idx, idx, i
            )
        )
        d = final_tree_cur.compare(bs_tree_array[i])
        dist += d["norm_rf"]
    avg = float(dist) / count
    with open(out_dir + "/self_dist.csv", "a") as w:
        w.write(str(idx) + "," + str(avg) + "\n")
    return avg


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
def run_snakemake(c, l, k, out_dir, run, roadies_dir):
    cmd = "snakemake --core {0} --use-conda --rerun-incomplete --config LENGTH={1} KREG={2} OUTDIR={3}".format(
        c, l, k, roadies_dir
    )
    os.system(cmd)
    # get the run output in folder
    print("Adding run to converge folder")
    os.system("./workflow/scripts/get_run.sh {0}/{1}".format(out_dir, run))


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

    print(boot_tree)

    return boot_tree


# function that returns an array of b bootstrapped newick trees
def bootstrap(b, out_dir, run, gene_trees):
    bs_trees = []

    pool = multiprocessing.Pool(b)

    args_list = [(i, out_dir, run, gene_trees) for i in range(b)]

    print(args_list)

    bs_trees = pool.starmap(process_bootstrap, args_list)

    pool.close()

    print(bs_trees)

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
def converge_run(i, l, k, c, out_dir, b, ref_exist, trees, roadies_dir):
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
    run_snakemake(c, l, k, out_dir, run, roadies_dir)
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
    os.system("rm -r " + out_dir)
    os.system("mkdir -p " + out_dir)
    os.system("mkdir " + out_dir + "/tmp")
    sys.setrecursionlimit(2000)

    # if there are input gene trees use that to build on or else make empty file
    master_gt = out_dir + "/master_gt.nwk"
    master_map = out_dir + "/master_map.txt"
    if input_gt is None or input_map is None:
        os.system("touch {0}".format(master_gt))
        os.system("touch {0}".format(master_map))
    else:
        os.system("cp {0} {1}".format(input_gt, master_gt))
        os.system("cp {0} {1}".format(input_map, master_map))

    # initialize lists for runs and distances
    runs = []
    self_dists = []
    iter_dists = []
    iter_dists_bs = []
    if ref_exist:
        ref_dists = []
    # list of roadies trees after each iteration
    trees = []
    # open files for writing distances
    if ref_exist:
        ref_out = open(out_dir + "/ref_dist.csv", "w+")
    iter_out = open(out_dir + "/iter_dist.csv", "w+")
    self_out = open(out_dir + "/self_dist.csv", "w+")
    iter_bs = open(out_dir + "/iter_dist_bs.csv", "w+")

    # starts main for loop for multiple iterations
    # for max iteration runs; start from 1 index instead of 0
    for i in range(max_iter):
        # returns an array of b bootstrapped trees
        run = converge_run(i, l, k, c, out_dir, b, ref_exist, trees, roadies_dir)
        runs.append(run)
        # get bootstrapped distance to itself
        # self_dist = comp_self(run,i)
        self_dist = comp_self_final(run, trees[i], i)
        print("Run " + str(i) + " average distance to itself: ", self_dist)
        self_dists.append(self_dist)
        self_out.write(str(i) + "," + str(self_dist) + "\n")
        # if reference exists get distance between ref and roadies tree
        if ref_exist:
            ref_dist = comp_tree(ref, trees[i])
            ref_dists.append(ref_dist)
            ref_out.write(str(i) + "," + str(ref_dist) + "\n")
        print("Average distance to self per iter so far: ", self_dists)
        for i in range(len(self_dists)):
            print("Run: " + str(i) + ": " + str(self_dists[i]))
        if ref_exist:
            for i in range(len(ref_dists)):
                print("Run: " + str(i) + ": " + str(ref_dists[i]))

        stop_run = False
        iter_flag = False

        # since comparing to previous, hence starts with 1st iteration
        if i >= 1:
            print(len(runs))
            print(runs)
            # get bootstrapped iteration distance
            # iter_dist_bs = comp_runs_bs(runs[i-1],runs[i],i)
            # get avg distance between current bootstrapped trees and previous iter final tree
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

            if i % 10 == 0:
                iterm = open(out_dir + "/iterm_iter_dist_bs.csv", "w")
                refm = open(out_dir + "/iterm_ref_dist.csv", "w")
                for j in range(len(iter_dists_bs)):
                    iterm.write(str(j) + "," + str(iter_dists_bs[j]) + "\n")
                for j in range(len(ref_dists)):
                    refm.write(str(j) + "," + str(ref_dists[j]) + "\n")
                iterm.close()
                refm.close()

            for i in range(len(iter_dists_bs)):
                print("Run: " + str(i) + ": " + str(iter_dists_bs[i]))
            print(
                "Distances to previous iterations final trees with bootstrapped trees so far: ",
                iter_dists_bs,
            )
            for i in range(len(iter_dists)):
                print("Run: " + str(i) + ": " + str(iter_dists[i]))
            # checking to see if we should stop script
            # every distance within the window must be lower than threshold to stop
            temp = 0
            temp2 = 0
            if i > 2 * stop_iter - 1:
                print("Seeing if we should stop")
                for j in range(stop_iter):
                    temp = temp + iter_dists_bs[i - j - 1]
                    temp2 = temp2 + iter_dists_bs[i - j - stop_iter - 1]

                temp = temp / stop_iter
                temp2 = temp2 / stop_iter
                if (temp2 - temp) < t:
                    print("crossed threshold")
                    # print("stop iteration number", j)
                    stop_run = True
                    # print("stop run value is", stop_run)
                else:
                    print("did not cross threshold")
                    stop_run = False
                    iter_flag = True
                    # print("stop run value is", stop_run)
                    # print("iter flag value is", iter_flag)
                    # break

            if stop_run == True and iter_flag == False:
                break
        if stop_run == True and iter_flag == False:
            break
    if ref_exist:
        ref_out.close()
    iter_out.close()
