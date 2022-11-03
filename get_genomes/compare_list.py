import sys
curr = sys.argv[1]
need = sys.argv[2]
current = []
with open(curr,'r') as f:
    lines = f.readlines()
    for l in lines:
        if len(l)>1:
            current.append(l.strip().replace('.fa',''))
        
needs = []
with open(need,'r') as f:
    lines = f.readlines()
    for l in lines:
        if len(l)>1:
            needs.append(l.strip())
for n in needs:
    if n not in current:
        print('needs',n)
        