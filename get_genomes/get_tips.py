import sys
path = sys.argv[1]
out = sys.argv[2]
with open(path,'r') as f:
    with open(out,'w') as w:
        lines = f.readlines()
        print(len(lines))
        for i in range(len(lines)):
            print(lines[i])
            line = lines[i]
            line = line.replace('(','')
            line = line.replace(')','')
            line = line.replace(':','')
            line = line.replace('.','')
            line = line.replace("'",'')
            line = ''.join([i for i in line if not i.isdigit()])
            print(line)
            split = line.split(',')
            print(split)
            for s in split:
                w.write(s+'\n')

