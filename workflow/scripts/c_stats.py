#libraries for plotting data
import glob, os, sys
from collections import OrderedDict
import matplotlib.pyplot as plt
from datetime import datetime
import subprocess
#input directory
dir = sys.argv[1]
#length of genes
length = sys.argv[2]
geneiter = sus.argv[3]
#function to parse stat file for x and y values
def get_run(statfile):
    run_dict=OrderedDict()
    with open(dir+'/'+statfile,'r') as f:
        lines=f.readlines()
        for i in range(len(lines)):
            if i%2==0:
                line = lines[i].strip()
                #print(line)
                cs = line.split(',')
                cs2 = cs[1].split('/')
                #print(cs2)
                idx = 0
                for j in range(len(cs2)):
                    if 'run' in cs2[j]:
                        idx = j
                cs3 = cs2[idx].replace('.nwk','')
                #print(cs3,idx)
                cs3 = cs3.replace('run_','')
                #print(cs3)
                id = int(cs3)*int(geneiter)
                run_dict[id] = float(lines[i+1])
    return run_dict

for f in ['refdist.txt','iterdist.txt','refdist2.txt','iterdist2.txt']:
    run_dict = get_run(f)
    lengths = list(run_dict.keys())
    #print(order)
    x = []
    y = []
    print(run_dict)
    for run in run_dict:
        x.append(run)
        y.append(run_dict[run])
    plt.plot(x,y,label='length='+str(length))
    plt.xlabel('Number of Genes')
    if '2' in f:
        plt.ylabel('RF Distance')
    else:
        plt.ylabel('Normalized Distance')
    if 'ref' in f:
        plt.title('Distance to Reference')
    else:
        plt.title('Distance to Previous Iter')
    plt.legend()
    plt.tight_layout()
    plt.savefig(dir+'/'+f.replace('txt','png'))
    plt.clf()
    
