import sys
import re
t1 = sys.argv[1]
t2 = sys.argv[2]
regex = re.compile('[^a-zA-Z_ ]')
tips = []
tips2 = []
with open(t1,'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        print(split)
        for s in split:
            s = regex.sub('', s)
            if len(s) > 1:
                tips.append(s)
print(tips)
with open(t2,'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        print(split)
        for s in split:
            s = regex.sub('', s)
            if len(s) > 1:
                tips2.append(s)
print(len(tips))
print(len(tips2))
for tip in tips:
    if tip not in tips2:
        print("not found in ",t1,tip)
for tip in tips2:
    if tip not in tips:
        print("not found in tree1",t2,tip)


