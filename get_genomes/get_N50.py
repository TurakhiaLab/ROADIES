import sys
from collections import OrderedDict
input = sys.argv[1]
mapfile = sys.argv[2]
datarank = sys.argv[3]
out = sys.argv[4]
names = []
name_convert = {}
with open(input,'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        sn = split[0]
        cn = split[1]
        names.append(sn)
        name_convert[sn] = cn
mapping = {}
mapping2 = {}
with open(mapfile,'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        GCA = split[0]
        sn = split[1]
        mapping[sn] = GCA
        mapping2[GCA] = sn
GCAS = []
for name in names:
    if name in mapping:
        GCAS.append((mapping[name],name))
    else:
        print("Not mapped: ",name)
N50s={}
with open(datarank,'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        id = split[0]
        N50 = split[1]
        N50s[id] = N50
orderbyN50 = {}
for g in GCAS:
    if g[0] in N50s:
        orderbyN50[g[0]] = int(N50s[g[0]])
orderbyN50 = {k: v for k, v in sorted(orderbyN50.items(), key=lambda item: item[1],reverse=True)}
print(orderbyN50)
with open(out,'w') as w:
    for g in orderbyN50:
        w.write(g + ',' + mapping2[g] + ','+ name_convert[mapping2[g]] + ',' + str(orderbyN50[g])+'\n')

        
