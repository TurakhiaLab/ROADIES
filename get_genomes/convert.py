import sys
path = sys.argv[1]
mapping = sys.argv[2]
out = sys.argv[3]
init = int(sys.argv[4])
fin = int(sys.argv[5])
with open(path,'r') as f:
    with open(mapping,'r') as m:
        with open(out,'w') as w:
            with open("names.txt",'w') as w2:
                mapping = {}
                lines = m.readlines()
                lines = lines[1:]
                for l in lines:
                    print(l)
                    split=l.strip().split(',')
                    print(split)
                    fr = split[init].strip()
                    to = split[fin].strip()
                    mapping[fr] = to
                print(mapping)
                nwk = f.readlines()
                for l in nwk:
                    for m in mapping:
                        if m in l:
                            l = l.replace(m,mapping[m])
                            w2.write(m+', '+mapping[m]+'\n')
                    w.write(l)
                    
            

        