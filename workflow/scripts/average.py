#python script to calculate averages from stat files with 2 columns
import sys
import numpy as np
stats = []
with open(sys.argv[1]) as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        stat = float(split[1])
        stats.append(stat)
print('mean = ' + "{:.2f}".format(np.mean(stats))+', med = '+ "{:.2f}".format(np.median(stats))+', std = '+"{:.2f}".format(np.std(stats)))

