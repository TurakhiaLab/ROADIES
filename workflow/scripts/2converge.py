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


# function that finds the average distance between one array of bootstrapped trees with final tree of previous iteration
def comp_bs_final(bs_tree_array, final_tree_prev, iteration):
    count = 0
    dist = 0
    for i in range(len(bs_tree_array)):
        count += 1
        print(
            "Comparing run {0} final tree of previous iteration {1} with bootstrap of iteration {2}".format(
                iteration, iteration - 1, i
            )
        )
        d = final_tree_prev.compare(bs_tree_array[i])
        dist += d["norm_rf"]
    avg = float(dist) / count
    with open(out_dir + "/iter_dist_bs.csv", "a") as w:
        w.write(str((iteration+1)*NUM_GENES) + "," + str(avg) + "\n")
    return avg



# function to run snakemake with settings and add to run folder
def run_snakemake(
    cores, out_dir, run,  weighted,runtime_left
):
    if weighted:
        cmd = ['snakemake','--core',str(cores),'--jobs',str(cores),'--use-conda','--rerun-incomplete']
    else:
        cmd = ['snakemake','--core',str(cores),'--jobs',str(cores),'--use-conda','--rerun-incomplete','--config','WEIGHTED=0']
    for i in range(len(cmd)):
        if i == len(cmd)-1:
            print(cmd[i])
        else:
            print(cmd[i],end=" ")
    if runtime_left == math.inf:
        print("Running without time limit")
        subprocess.run(cmd)
    else:
        print("Trying to fit in run with time remaining",runtime_left)
        try:
            subprocess.run(cmd,timeout=runtime_left)
        except:
            print("TIME LIMIT EXCEEDED")
            sys.exit(0)
    # get the run output in folder
    print("Adding run to converge folder")
    os.system(
        "./workflow/scripts/get_run.sh {0} {1}".format(
            out_dir, run
        )
    )


# function that returns an array of b bootstrapped newick trees
def process_bootstrap(iteration, out_dir, run, gene_trees):
    print("Creating bootstrapping tree: ", iteration)
    # path to temp gt
    tmp_path = out_dir + "/tmp/" + run + "." + str(iteration)
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
def bootstrap(num_bootstrap, out_dir, run, gene_trees):
    bs_trees = []
    pool = multiprocessing.Pool(num_bootstrap)
    args_list = [(i, out_dir, run, gene_trees) for i in range(num_bootstrap)]
    bs_trees = pool.starmap(process_bootstrap, args_list)
    pool.close()

    return bs_trees


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
def converge_run(iteration,cores,out_dir,num_bootstrap,ref_exist,trees,roadies_dir,weighted,runtime_left):
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
    run_snakemake(cores, out_dir, run,  weighted, runtime_left)
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
    return bootstrap(num_bootstrap, out_dir, run, gene_trees),len(gene_trees)


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
    genomes= config["GENOMES"]
    NUM_GENOMES = len(os.listdir(genomes))
    gene_mult = config["GENE_MULT"]
    NUM_GENES = gene_mult * NUM_GENOMES
    LENGTH = config["LENGTH"]
    NUM_BOOTSTRAP = config["NUM_BOOTSTRAP"]
    input_gt = config["INPUT_GENE_TREES"]
    input_map = config["INPUT_MAP"]
    MIN_ITER = config["MIN_ITER"]
    MAX_RUNTIME = config["MAX_RUNTIME"]
    roadies_dir = config["OUT_DIR"]
    species_ids = "species_ids.csv"
    gene_ids = "gene_ids.csv"
    species_lists = "species_lists.csv"
    
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
    os.system("mkdir -p " + out_dir + "/tmp")
    sys.setrecursionlimit(2000)
    # initialize lists for runs and distances
    runs = []
    iter_dists = []
    iter_dists_bs = []
    window_dists = []
    time_stamps = []
    if ref_exist:
        ref_dists = []
    # list of roadies trees after each iteration
    trees = []
    # open files for writing distances
    
    # starts main for loop for multiple iterations
    # for max iteration runs; start from 1 index instead of 0
    iteration = 0
    start_time = time.time()
    start_time_l = time.asctime(time.localtime(time.time()))
    time_stamps.append(start_time)
    runtime_left = math.inf
    gt_counts = []
    with open(out_dir+"/time_stamps.csv",'a') as t_out:
        t_out.write("Start time: "+str(start_time_l)+"\n")
    while(True):
        # returns an array of b bootstrapped trees
        weighted = True
        if iteration == 0 and input_gt is None:
            weighted = False
        run,num_gt = converge_run(iteration,CORES,out_dir,NUM_BOOTSTRAP,ref_exist,trees,roadies_dir,weighted,runtime_left)
        print("There are {0} gene trees after iteration {1}".format(num_gt,iteration))
        runs.append(run)
        gt_counts.append(num_gt)
        # if reference exists get distance between ref and roadies tree
        if ref_exist:
            ref_dist = comp_tree(ref, trees[iteration])
            ref_dists.append(ref_dist)
            with open(out_dir + "/ref_dist.csv", "a") as ref_out:
                ref_out.write(str(num_gt) + "," + str(ref_dist) + "\n")
        stop_run = False
        iter_flag = False

        # since comparing to previous, hence starts with 1st iteration
        if iteration >= 1:
            print(len(runs))
            print(runs)
            # get bootstrapped iteration distance
            # iter_dist_bs = comp_runs_bs(runs[i-1],runs[i],i)
            # get avg distance between current bootstrapped trees and previous iter final tree
            print(
                "Getting bootstrapped iteration distances betseen runs {0} and {1}".format(
                    iteration, iteration - 1
                )
            )
            iter_dist_bs = comp_bs_final(runs[iteration], trees[iteration - 1], iteration)
            print(
                "Average distance between final tree of iteration {0} and bootstrapped tree of iteration {1} is: {2}".format(
                    iteration - 1, iteration, iter_dist_bs
                )
            )
            iter_dists_bs.append(iter_dist_bs)
            
            
            # get single iteration distance
            iter_dist = comp_tree(trees[iteration], trees[iteration - 1])
            print(
                "Distance between final trees of iteration {0} and {1} is: {2}".format(
                    iteration - 1, iteration, iter_dist
                )
            )
            iter_dists.append(iter_dist)
            with open(out_dir + "/iter_dist.csv", "a") as iter_out:
                iter_out.write(str(num_gt) + "," + str(iter_dist) + "\n")
            # sometimes dist files get updated slowly so create intermediate file for analysis in meantime
    
            # checking to see if we should stop script
            # every distance within the window must be lower than threshold to stop
            temp = 0
            temp2 = 0
            # get the difference of average distances between two windows of size stop_iter
            print("Distance to previous Tree")
            for j in range(len(iter_dists)):
                print("Run: " + str(j) + " ("+str((j+1)*NUM_GENES)+"): " + str(iter_dists[j]))
            print(
                "Distances to previous iterations final trees with bootstrapped trees so far:"
            )
            for j in range(len(iter_dists_bs)):
                print("Run: " + str(j) +" ("+str((j+1)*NUM_GENES)+"): " + str(iter_dists_bs[j]))
            if ref_exist:
                print("Distance to reference so far:")
                for j in range(len(ref_dists)):
                    print("Run: " + str(j) + " ("+str((j+1)*NUM_GENES)+"): " + str(ref_dists[j]))
        curr_time = time.time()
        curr_time_l = time.asctime(time.localtime(time.time()))
        to_previous = curr_time - time_stamps[len(time_stamps)-1]
        time_stamps.append(curr_time)
        elapsed_time = curr_time - start_time
        with open(out_dir+"/time_stamps.csv",'a') as t_out:
            t_out.write(str(num_gt)+','+str(curr_time_l)+','+str(elapsed_time)+','+str(to_previous)+"\n")
        if iteration+1 >= MIN_ITER:
            if elapsed_time >= MAX_RUNTIME:
                print("Current elapsed time",elapsed_time)
                print("Max Runtime",MAX_RUNTIME)
                break
            else:
                runtime_left = MAX_RUNTIME-elapsed_time
                print("Current elapsed time",elapsed_time)
                print("Max Runtime",MAX_RUNTIME)
                print("Runtime left",runtime_left)
        
        iteration+=1

          
