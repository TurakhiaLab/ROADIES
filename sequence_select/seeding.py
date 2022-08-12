import argparse, os
from collections import OrderedDict
import random
def dir_path(str):
    if os.path.isdir(str):
        return str
    else:
        raise NotADirectoryError(str)
parser = argparse.ArgumentParser(description='Seed genomes')
parser.add_argument('-l',type=int,default = 1000)
parser.add_argument('-k',type=int,default=300)
parser.add_argument('--path',type=dir_path,default="./../genomes")
parser.add_argument('--output',default="./index.csv")
args = parser.parse_args()
length = args.l
num = args.k
directory = args.path
dir_list = os.listdir(directory)
num_genomes = len(dir_list)
# print("There are "+ str(len(dir_list))+" genomes")
od = OrderedDict()
for i in range(num):
    index = random.randint(0,num_genomes-1)
    if index in od.keys():
        od[index] = od[index]+1
    else:
        od[index] = 1
od = OrderedDict(sorted(od.items()))
f = open(args.output,'w')
for key,value in od.items():
    # print(str(key)+ ", "+dir_list[key] +", " +str(value))
    f.write(dir_list[key]+","+str(value)+'\n')

