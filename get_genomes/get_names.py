from pickle import FALSE
import pandas as pd
import sys
import glob,os
#path = sys.argv[1]
def get_id(id):
    s = id.split('_')
    return s[0] +'_' +s[1]
def remove_space(name):
    return name.replace(' ','_')
def format_name(name):
    s = ""
    print("unformatted:",name)
    if '(' in name or ')' in name:
        split = name.split()
        print(split)
        within_p = False
        for i in range(len(split)):
            #print(sp[i])
            if '(' in split[i]:
                within_p = True
                continue
            if ')' in split[i]:
                within_P = False
                continue
            if within_p:
                continue
            if len(s) == 0:
                s = s + split[i]
            else:
                s = s+ "_"+split[i]
        return s
    else:
        return name.replace(' ','_')

url = 'https://hgdownload.soe.ucsc.edu/hubs/birds/index.html'
tables = pd.read_html(url)
table = tables[1]
#print(table.head)
common_name = table.iloc[:,1]
common_name = common_name.apply(format_name)
#print(common_name)
science_name = table.iloc[:,2]
science_name = science_name.apply(remove_space)
ids = table.iloc[:,3]
ids = ids.apply(get_id)
sn = science_name.tolist()
id = ids.tolist()
cn = common_name.tolist()
names = list(zip(sn,cn))
#print(names)
mapping = dict(zip(id,names))
print(mapping)
with open("species_map.txt",'w') as w:
    w.write("Accession Number, scientific name, common name\n")
    for map in mapping:
        w.write(map+','+mapping[map][0]+','+ mapping[map][1]+'\n')
#for filename in glob.glob(os.path.join(path,'*.fa')):
 #   s = filename.split('/')
  #  ext = s[0]+'/'+s[1] + '/'
   # f = s[2]
    #id = f.replace('.fa','')
    #new_loc = ext+mapping[id]+'.fa'
    #print(new_loc)
    #os.rename(filename,new_loc)