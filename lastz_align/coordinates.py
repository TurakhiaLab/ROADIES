import re
import os, glob
import sys
p = sys.argv[1]
o = sys.argv[2]
output = "Coordinates"
path = os.path.join(o,output)
print(p,o,output,path)
if not os.path.exists(path):
    os.mkdir(path)
for filename in glob.glob(os.path.join(p,'*.maf')):
    with open(os.path.join(os.getcwd(),filename),'r') as f:
        s = filename.split('/')
        s2 = s[1].split('.')
        genome = s2[0]+'.'+s2[1]
        lines = f.readlines()
        for l in range(15,len(lines)):
            if (l-15)%4 == 0:
                split = lines[l].split()
                seq = split[len(split)-1]
                seq = re.sub('-','',seq)
                chrnum = split[1]
                index = split[2]
                name = split[1]+":"+index
                split2 = lines[l+1].split()
                seq2 = split2[1]
                w = open(path+'/'+seq2+'.m.fa','a')
                w.write('>'+genome+'.'+name+'\n')
                w.write(seq +'\n')
