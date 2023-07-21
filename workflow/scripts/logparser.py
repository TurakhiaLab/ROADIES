# logparser.py reads in a snakemake log file and parses the runtime for each rule
# REQUIRES: pandas, matplotlib, seaborn
# USAGE: `python workflow/scripts/logparser.py [log_file] [output]`
import sys
import datetime
import time
log_file = sys.argv[1]
output = sys.argv[2]
rules = []
times = []
pasta = {}
iqtree = {}
parallel = ["pasta", "iqtree"]
short_rules = ["filtermsa"]
def to_seconds(time_line):
    without_bars = time_line.replace('[','').replace(']','')
    to_datetime = datetime.datetime.strptime(without_bars.strip(),"%a %b %d %H:%M:%S %Y")
    epoch_time = datetime.datetime(1970, 1, 1)
    delta = (to_datetime-epoch_time)
    seconds= delta.total_seconds()
    return seconds

with open(log_file, "r") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        if "rule " in lines[i]:
            the_rule = lines[i].strip().replace(":", "")
            rule = the_rule.split()[1]
            if rule in short_rules:
                continue
            # print(rule)
            if rule in parallel:
                jobid_line = lines[i + 3]
                s = jobid_line.split()
                id = s[1]
                idx = 0
                if "[" in lines[i - 1]:
                    idx = i - 1
                else:
                    idx = i - 2
                time_line = lines[idx]
                if rule == "pasta":
                    pasta[id] = [to_seconds(time_line)]
                else:
                    iqtree[id] = [to_seconds(time_line)]
            else:
                if rule in rules:
                    continue
                else:
                    # print("rule: ",rule,i)
                    rules.append(rule)
                    time_line = lines[i - 1]
                    # print(time_line)
                    the_time = to_seconds(time_line)
                    # print(time)
                    times.append(the_time)
                   # print(rule,the_time)
        elif "Finished" in lines[i]:
            s = lines[i].split()
            num = s[2].replace(".", "")
            if num in pasta:
              #  print('end pasta')
                time_line = lines[i - 1]
                # print(time_line)
               
                pasta[num].append(to_seconds(time_line))
                # print(mafft[num])
            elif num in iqtree:
                time_line = lines[i - 1]
                iqtree[num].append(to_seconds(time_line))
                # print(iqtree[num])
#print("before",times)
t = list(pasta.values())
print(t)
s = t[0][0]
e = t[0][1]
for j in t:
    s1 = j[0]
    e2 = j[1]
    if s1 < s:
        print("s2s",s1,s)
        s = s1
    if e2 > e:
        print("e2e",e2,e)
        e = e2
#print(s,e)
t2 = list(iqtree.values())
rules.insert(5, "treebuilding")
times.insert(5, s)
for i in range(len(rules)):
    print(i,rules[i],times[i])
#print(times)
x = rules[: len(rules) - 1]
runtimes = []
with open(output, "w") as w:
    w.write(log_file + "\n")
    for i in range(len(rules) - 1):
        start = times[i]
        end = times[i + 1]
        rule = rules[i]
        total = end - start
        the_time = str(datetime.timedelta(seconds=(end - start)))
        print(start,end,total,the_time)
        print(rule + ": " + the_time)
        w.write(rules[i] + ", " + the_time + "\n")
        runtimes.append(total)
    start = times[0]
    end = times[len(times) - 1]
    the_time = str(datetime.timedelta(seconds=(end - start)))
    print("total_time: " + the_time)
    w.write("total_time, " + the_time + "\n")
    x.append("total")
    runtimes.append(end - start)