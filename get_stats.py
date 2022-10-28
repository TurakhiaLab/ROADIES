#import need
import glob, os, sys
import matplotlib.pyplot as plt
from datetime import datetime
import numpy as np
import matplotlib.dates as mdates
#argument for path to directory and output dir
path = sys.argv[1]
dir = sys.argv[2]
iter = []
runs = []
times = []
#command = "Rscript"
#path2script = "dist.R"
for subdir,dirs, files in os.walk(path):
    for file in files:
        #print(os.path.join(subdir, file))
        if file == "runtime.txt":
            with open(os.path.join(subdir,file),'r') as f:
                print(subdir,file)
                #print(subdir)
                split = subdir.split('/')
                run = split[2]
                runs.append(run)
                lines = f.readlines()
                total = lines[len(lines)-1]
                ts = total.split(',')
                time = ts[1].strip()
                print(time)
                format_time = "%H:%M:%S"
                dt = datetime.strptime(time,format_time)
                times.append(dt)
                #print(total)

ds = dict(zip(runs,times))
run_dict = {}
print(ds)
#print(ds)
for d in ds:
    split = d.split('_')
    params = split[1]
    ps = params.split('.')
    #print(ps)
    print(ps[0],ps[1])
    l = int(ps[0].replace('L',''))
    k=0
    if 'A' in ps[1]:
        k = int(ps[1].replace('KA',''))
    else:
        k = int(ps[1].replace('K',''))
    if l not in run_dict:
        run_dict[l] = [(k,ds[d])]
    else:
        run_dict[l].append((k,ds[d]))
for r in run_dict:
    run_dict[r]=sorted(run_dict[r])
    #print(run_dict[r])
    x = [x[0] for x in run_dict[r]]
    y = [y[1] for y in run_dict[r]]
    with open(dir+'/runtimes.txt','w') as w:
        for i in range(len(x)):
            s= str(y[i]).split()
            time = s[1]
            w.write(str(x[i])+','+time+'\n')
    #print(x)
    #print(y)
    plt.plot(x,y,label='L='+str(r))
lengths = list(run_dict.keys())
order = np.argsort(lengths)
handles, labels = plt.gca().get_legend_handles_labels()
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order])
plt.show()
plt.title('Runtime in hours')
ax=plt.gca()
yfmt = mdates.DateFormatter('%H:%M:%S')
ax.yaxis.set_major_formatter(yfmt)
plt.xlabel('# genes')
plt.ylabel("minutes")
plt.savefig(dir+'/runtimes.png')
plt.clf()
run_dict = {}
with open(dir+"/refdist.txt",'r') as f:
    lines=f.readlines()
    for i in range(len(lines)):
        if i%2==0:
            params=lines[i].strip()
            dist = float(lines[i+1].strip())
            ps = params.split('/')
            print(ps[2])
            print(dist)
            split = ps[2].split('_')
            params = split[1]
            ps = params.split('.')
            print(ps)
            l = int(ps[0].replace('L',''))
            k=0
            if 'A' in ps[1]:
                k = int(ps[1].replace('KA',''))
            else:
                k = int(ps[1].replace('K',''))
            if l not in run_dict:
                run_dict[l] = [(k,dist)]
            else:
                run_dict[l].append((k,dist))
print(run_dict)
for r in run_dict:
    run_dict[r]=sorted(run_dict[r])
    #print(run_dict[r])
    x = [x[0] for x in run_dict[r]]
    y = [y[1] for y in run_dict[r]]
    #print(x)
    #print(y)
    plt.plot(x,y,label='L='+str(r))
plt.legend()
plt.title('Distance from TimeTree')
plt.show()
plt.xlabel('# genes')
plt.ylabel("distance")
plt.savefig(dir+"/refdist.png")
plt.clf()
run_dict = {}
with open(dir+"/refdist2.txt",'r') as f:
    lines=f.readlines()
    for i in range(len(lines)):
        if i%2==0:
            params=lines[i].strip()
            dist = int(lines[i+1].strip())
            ps = params.split('/')
            print(ps[2])
            print(dist)
            split = ps[2].split('_')
            params = split[1]
            ps = params.split('.')
            print(ps)
            l = int(ps[0].replace('L',''))
            k=0
            if 'A' in ps[1]:
                k = int(ps[1].replace('KA',''))
            else:
                k = int(ps[1].replace('K',''))
            if l not in run_dict:
                run_dict[l] = [(k,dist)]
            else:
                run_dict[l].append((k,dist))
print(run_dict)
for r in run_dict:
    run_dict[r]=sorted(run_dict[r])
    #print(run_dict[r])
    x = [x[0] for x in run_dict[r]]
    y = [y[1] for y in run_dict[r]]
    #print(x)
    #print(y)
    plt.plot(x,y,label='L='+str(r))
plt.legend()
plt.title('Distance from TimeTree RF')
plt.show()
plt.xlabel('# genes')
plt.ylabel("distance")
plt.savefig(dir+"/refdist2.png")
