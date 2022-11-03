import sys, os
mapping = sys.argv[1]
path = sys.argv[2]
to = sys.argv[3]
m = {}
with open(mapping,'r') as f:
    lines = f.readlines()
    for l in lines:
        split = l.strip().split(',')
        m[split[0]] = split[1].strip()
print(m)
for filename in os.listdir(path):
    print(filename)
    name = filename.replace('.fa','')
    if name in m:
        print('found')
        print(m[name])
        os.link(path+'/'+filename,to+'/'+m[name]+'.fa')