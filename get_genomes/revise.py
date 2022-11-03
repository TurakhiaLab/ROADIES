import sys
import glob,os
path = sys.argv[1]
mapfile = sys.argv[2]
mapping = {}
with open(mapfile,'r') as f:
    lines = f.readlines()
    for line in lines:
        split = line.strip().split(',')
        mapping[split[0]]=split[1]
    print(mapping)
    for filename in os.listdir(path):
        print(filename)
        name = filename.replace('.fa','')
        if name in mapping:
            print('found')
            os.rename(path+'/'+filename,path+'/'+mapping[name]+'.fa')