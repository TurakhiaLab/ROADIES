#input directory from either get_stats or cstats
#python combined.py [stat dir] [outdir] [auto/converge]
import glob, os, sys
from re import L
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import matplotlib.dates as mdates
path = sys.argv[1]
outdir = sys.argv[2]
conv = True
if sys.argv[3] == 'autorun':
    conv = False
    print('autorun')
elif sys.argv[3] == 'converge':
    print('converge')
else:
    exit(1)
def get_runs(statfile):
    run_dict=OrderedDict()
    for subdir,dirs, files in os.walk(path):
        for file in files:
            dir = os.path.join(subdir,file)
            #print(file)
            if file == statfile:
                #print(os.path.join(subdir, file))
                with open(dir,'r') as f:
                    split = subdir.split('/')
                    idx = 0
                    for i in range(len(split)):
                        if 'length' in split[i]:
                            idx =i
                    length = int(split[i].replace('length',''))
                    #print(length)
                    run_dict[length] = []
                    lines = f.readlines()
                    if statfile == 'runtimes.txt':
                        for l in lines:
                            split = l.strip().split(',')
                            k = int(split[0])
                            time = split[1].strip()
                            print(time)
                            format_time = "%H:%M:%S"
                            dt = datetime.strptime(time,format_time)
                            run_dict[length].append((k,dt))

                    else:
                        for i in range(len(lines)):
                            if i%2==0:
                                line = lines[i].strip()
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
                                if not conv:
                                    s = cs3.split('.')
                                    cs3 = s[1].replace('K','')
                                    cs3 = cs3.replace('A','')
                                #print(cs3)
                                id = int(cs3)
                                if conv:
                                    id = id *200
                                run_dict[length].append((id,float(lines[i+1])))
                           # print(id,lines[i+1])
    #print(run_dict)
    return run_dict


for f in ['refdist.txt','iterdist.txt','refdist2.txt','iterdist2.txt','runtimes.txt']:
    print(f)
    if f == 'runtimes.txt' and conv:
        continue
    run_dict = get_runs(f)
    lengths = list(run_dict.keys())
    order = np.argsort(lengths)
    #print(order)
    for run in run_dict:
        #print(run_dict[run])
        x = []
        y = []
        for i in range(len(run_dict[run])):
            #print(run_dict[run][i])
            x.append(run_dict[run][i][0])
            y.append(run_dict[run][i][1])
        plt.plot(x,y,label='length='+str(run))
    plt.xlabel('Number of Genes')
    ax=plt.gca()
    if '2' in f:
        plt.ylabel('RF Distance')
    elif 'run' in f:
        plt.ylabel('Total Time')
    else:
        plt.ylabel('Normalized Distance')
    if 'ref' in f:
        plt.title('Distance to Reference')
    elif 'run' in f:
        plt.title("Runtimes")
        yfmt = mdates.DateFormatter('%H:%M:%S')
        ax.yaxis.set_major_formatter(yfmt)
    else:
        plt.title('Distance to Previous Iter')
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
    plt.tight_layout()
    plt.savefig(outdir+'/'+f.replace('txt','png'))
    plt.clf()
