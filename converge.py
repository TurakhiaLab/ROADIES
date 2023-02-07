#requirements ETE and snakemake only
import os
import argparse
import random
import subprocess
import signal
from ete3 import Tree 

#function that finds the average distance between an array of trees and itself
def comp_ref(run,ref,i):
    dist = 0
    print(run)
    for i in range(len(run)):
        d = run[i].compare(ref)
        dist += d['norm_rf']
    avg = float(dist)/len(run)
    with open(out_dir+'/ref_dist.csv','a') as w:
        w.write(str(i)+','+str(avg)+'\n')
    return avg
def comp_self(run,i):
    dist = 0
    count = 0
    print(run)
    for i in range(len(run)-1):
        for j in range(i+1,len(run)):
            count +=1
            d = run[i].compare(run[j])
            dist += d['norm_rf']
    avg = float(dist)/count
    with open(out_dir+'/self_dist.csv','a') as w:
        w.write(str(i)+','+str(avg)+'\n')
    return avg

#function that finds the average distance between two arrays of trees
def comp_runs(run1,run2,i):
    count = 0
    dist = 0
    for i in range(len(run1)):
        for j in range(len(run2)):
            count += 1
            d = run1[i].compare(run2[j])
            dist += d['norm_rf']
    avg = float(dist)/count
    with open(out_dir+'/iter_dist.csv','a') as w:
        w.write(str(i)+','+str(avg)+'\n')
    return avg
class Alarm(Exception):
    pass
def alarm_handler(*args):
    raise Alarm("timeout")
#function to run snakemake with settings and add to run folder
def run_snakemake(c,l,k,out_dir,run):

    cmd ='snakemake --core {0} --use-conda --rerun-incomplete --config LENGTH={1} KREG={2}'.format(c,l,k)
    os.system(cmd)
    #get the run output in folder
    print("Adding run to converge folder")
    os.system('./workflow/scripts/get_run.sh {0}/{1}'.format(out_dir,run))

#function that returns an array of b bootstrapped newick trees 
def bootstrap(b,out_dir,run,gene_trees):
    bs_trees = []
    for i in range(b):
        
        print("Creating bootstrapping tree: ",i)
        #path to temp gt
        tmp_path = out_dir+'/tmp/'+run+'.'+str(i)
        out = open(tmp_path+'.gt.nwk','w')
        leaves = []
        w = open(tmp_path+'.map.txt','w')
        #sample a random line and output to tmp file
        for k in range(len(gene_trees)):
            n = random.randint(0,len(gene_trees)-1)
            out.write(gene_trees[n])
            n = Tree(gene_trees[n])
            leaves = n.get_leaf_names()
            for leaf in leaves:
                w.write(leaf+' ')
                s = leaf.split('_')
                w.write(s[0]+'\n')
                print(leaf +' '+ s[0])
        out.close()

        #run astral on bootstrapped tree
        os.system('ASTER-Linux/bin/astral-pro -i {0}.gt.nwk -o {0}.nwk -a {0}.map.txt'.format(tmp_path))
        boot_tree=Tree(tmp_path+'.nwk')
        bs_trees.append(boot_tree)
    return bs_trees

def combine_iter(out_dir,run):
    print("Concatenating run's gene trees and mapping files with master versions")
    os.system('cat {0}/{1}/gene_tree_merged.newick >> {0}/master_gt.nwk'.format(out_dir,run))
    os.system('cat {0}/{1}/mapping.txt >> {0}/master_map.txt'.format(out_dir,run))
    #open both files and get lines, each line is a separate gene tree
    os.system('ASTER-Linux/bin/astral-pro -i {0}/master_gt.nwk -o {0}/{1}.nwk -a {0}/master_map.txt'.format(out_dir,run))
    #open both master files and get gene trees and mapping
    gt = open(out_dir+'/master_gt.nwk','r')
    m = open(out_dir+'/master_map.txt','r')
    gene_trees = gt.readlines()
    mapping = m.readlines()
    gt.close()
    m.close()
    return gene_trees,mapping

def converge_run(i,l,k,c,out_dir,b):
    os.system('rm -r results')
    os.system('mkdir results')
    run = "run_"
    #allows sorting runs correctly
    if i < 10:
        run += "0" + str(i)
    else:
        run += str(i)
    print("Starting " +run)
    #run snakemake with specificed gene number and length
    run_snakemake(c,l,k,out_dir,run)
    #merging gene trees and mapping files
    gene_trees,mapping = combine_iter(out_dir,run)
    #create bootstrapping trees
    return bootstrap(b,out_dir,run,gene_trees)

if __name__=="__main__":
#taking in arguments, have default values for most
    parser = argparse.ArgumentParser(
        prog = 'Converge',
        description = 'Script to continuously run snakemake with a small number of genes combining the gene trees after each run'
    )
    parser.add_argument('--input_gt',default=None,help='can put in previously generated gene trees')
    parser.add_argument('--input_map',default=None,help='can put in previously generated mapping')
    parser.add_argument('--ref',default='trees/cn48.nwk',help='reference tree (input as .nwk or .newick')
    parser.add_argument('-c',type=int,default=16,help='number of cores')
    parser.add_argument('-k',type=int,default=100,help='number of genes')
    parser.add_argument('-t',type=float,default=0.1,help=' maximum TreeDistance threshold')
    parser.add_argument('-l',type=int,default=500,help='length of genes')
    parser.add_argument('--bootstrap',type=int,default=10,help='number of trees for bootstrapping when comparing')
    parser.add_argument('--max_iter',type=int,default=5,help='maximum number of runs before stopping')
    parser.add_argument('--stop_iter',type=int,default=3,help='number of runs satisfying threshold before halt')
    parser.add_argument('--out_dir',default='converge',help='output dir')
    parser.add_argument('--smk_dir',default='results',help='Snakemake output directory')
    #assigning argument values to variables
    args = vars(parser.parse_args())
    ref = Tree(args['ref'])
    c = args['c']
    k = args['k']
    t = args['t']
    l = args['l']
    b= args['bootstrap']
    input_gt = args['input_gt']
    input_map = args['input_map']
    max_iter = args['max_iter']
    stop_iter= args['stop_iter']
    out_dir = args['out_dir']
    smk_dir=args['smk_dir']
    os.system('rm -r '+out_dir)
    os.system('mkdir -p '+out_dir)
    os.system('mkdir '+out_dir+'/tmp')
    # if there are input gene trees use that to build on or else make empty file
    master_gt = out_dir+'/master_gt.nwk'
    master_map = out_dir+'/master_map.txt'
    if input_gt is None or input_map is None:
        os.system('touch {0}'.format(master_gt))
        os.system('touch {0}'.format(master_map))
    else:
        os.system('cp {0} {1}'.format(input_gt,master_gt))
        os.system('cp {0} {1}'.format(input_map,master_map))
    #initialize lists for runs and distances
    runs= []
    self_dists=[]
    iter_dists=[]
    ref_dists = []
    #for max iteration runs; start from 1 index instead of 0
    for i in range(3):
        
        #returns an array of b bootstrapped trees
        run = converge_run(i,l,k,c,out_dir,b)
        runs.append(run)
        self_dist = comp_self(run,i)
        print("Average distance to itself: ",self_dist)
        self_dists.append(self_dist)
        ref_dist = comp_ref(run,ref,i)
        print("Average distance to itself: ",ref_dist)

        ref_dists.append(ref_dist)
        print("Average distance to self per iter: ",self_dists)
        print("Average distance to ref per iter: ",ref_dists)

        stop_run = False

        if i >= 1:
            print(len(runs))
            print(runs)
            #compare between iterations, indexes due to starting from 1 index
            iter_dist = comp_runs(runs[i-1],runs[i],i)
            print("Average distance between iteration {0} and {1} is: {2}".format(i-1,i,iter_dist))
            iter_dists.append(iter_dist)
            if i > stop_iter-1:
                print("Seeing if We should stop")
             
                for j in range(stop_iter):
                    if self_dists[i-j] < t or iter_dists[i-j-1] < t:
                        stop_run = True
                        break
                if stop_run:
                    break
                print("Average distance to prev per iter: ",iter_dists)

            if stop_run:
                break
        if stop_run:
            break
