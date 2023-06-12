#logparser.py reads in a snakemake log file and parses the runtime for each rule
#REQUIRES: pandas, matplotlib, seaborn
#USAGE: `python workflow/scripts/logparser.py [log_file] [output]`
import sys
import datetime
log_file = sys.argv[1]
output = sys.argv[2]
rules = []
times = []
pasta = {}
iqtree = {}
parallel = ['pasta', 'iqtree']
def to_seconds(time):
    s = time.split(':')
    secs = int(s[2])+(60*int(s[1]))+(3600*int(s[0]))
    return secs
with open(log_file,'r') as f:
    lines = f.readlines()
    for i in range(len(lines)):
        if 'rule ' in lines[i]:
            the_rule = lines[i].strip().replace(':','')
            rule = the_rule.split()[1]
            #print(rule)
            if rule in parallel:
                jobid_line = lines[i+3]
                s = jobid_line.split()
                id = s[1]
                idx = 0
                if '[' in lines[i-1]:
                    idx= i-1
                else:
                    idx = i-2
                time_line = lines[idx]
                ts = time_line.split()
                time = ts[3]
                if rule == 'pasta':
                    pasta[id] = [to_seconds(time)]
                else:
                    iqtree[id] = [to_seconds(time)]
            else:
                if rule in rules:
                    continue
                else: 
                    #print("rule: ",rule,i)
                    rules.append(rule)
                    time_line = lines[i-1]
                    #print(time_line)
                    ts = time_line.split()
                    time = to_seconds(ts[3])
                    #print(time)
                    times.append(time)
        elif 'Finished' in lines[i]:
            s = lines[i].split()
            num = s[2].replace('.','')
            if num in pasta:
                #print('end pasta')
                time_line = lines[i-1]
                #print(time_line)
                ts = time_line.split()
                time = ts[3]
                [num].append(to_seconds(time))
                #print(mafft[num])
            elif num in iqtree:
                time_line = lines[i-1]
                #print(time_line)
                ts = time_line.split()
                time = ts[3]
                iqtree[num].append(to_seconds(time))
                #print(iqtree[num])
t = list(pasta.values())
print(t)
s = t[0][0]
e = t[0][1]
for j in t:
    s1 = j[0]
    e2 = j[1]
    if s1 < s:
        s = s1
    if e2 > e:
        e = e2
t2 = list(iqtree.values())
rules.insert(4,'pasta')
times.insert(4,s)
rules.insert(5,'iqtree')
times.insert(5,e)
x = rules[:len(rules)-1]
runtimes= []
with open(output+'/statistics/runtime.txt','w') as w:
    w.write(log_file+'\n')
    for i in range(len(rules)-1):
        start = times[i]
        end = times[i+1]
        rule = rules[i]
        total = end -start
        time = str(datetime.timedelta(seconds=(end-start)))
        print(rule+': '+time)
        w.write(rules[i]+', '+time+'\n')
        runtimes.append(total)
    start = times[0]
    end = times[len(times)-1]
    time = str(datetime.timedelta(seconds=(end-start)))
    print('total_time: '+time)
    w.write('total_time, '+time+'\n')
    x.append('total')
    runtimes.append(end-start)
