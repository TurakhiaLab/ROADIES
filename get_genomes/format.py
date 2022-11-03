import sys
with open(sys.argv[1],'r') as f:
    with open(sys.argv[2],'w') as w: 
        lines = f.readlines()
        for l in lines:
            l = l.replace("'",'')
            l = l.replace(' ','_')
            w.write(l)
