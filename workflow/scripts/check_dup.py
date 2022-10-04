import glob, os
import sys
path = sys.argv[1]
for filename in glob.glob(os.path.join(path,'*.fa')):
    with open(os.path.join(os.getcwd(),filename),'r') as f:
        lines=f.readlines()
        nodup=True
        isin = []
        for i in range(len(lines)):
            if i %2 == 0:
                line = lines[i]
                if line in isin:
                    print("found dup")
                    nodup=False
                else:
                    isin.append(line)
        if nodup:
            print("No Duplicates")