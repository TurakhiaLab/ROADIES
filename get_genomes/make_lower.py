import sys
with open(sys.argv[1],'r') as f:
    with open(sys.argv[2],'w') as w:
        lines = f.readline()
        for l in lines:
            l = l.lower()
            w.write(l)
        