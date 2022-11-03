import sys
import re
path = sys.argv[1]
science = sys.argv[2]
common = sys.argv[3]
both = sys.argv[4]
regex = re.compile('[^a-zA-Z ]')
mapping = {}
with open(path,'r') as f:
    with open(science,'w') as w:
        with open(common,'w') as w2:
            with open(both,'w')as w3:
                lines = f.readlines()
                for l in lines:
                    split = l.split(',')
                    print(split)
                    for s in split:
                        s = regex.sub('', s)
                        print(s)
                        s2 = s.split()
                        isScience=True
                        science_name = s2[0]
                        common_name = ""
                        for i in range(1,len(s2)):
                            if s2[i][0].islower() and isScience:
                                science_name +='_'+s2[i]
                            else:
                                isScience = False
                                if i == len(s2)-1:
                                    common_name += s2[i]
                                else:
                                    common_name += s2[i] +'_'
                        mapping[s] = science_name
                        print("science name",science_name)
                        w.write(science_name+'\n')
                        print("common name",common_name)
                        w2.write(common_name+'\n')
                        w3.write(science_name+', '+common_name+'\n')

                        print(s)
                #w.write(s+'\n')
print(mapping)    
with open(path,'r') as f:
    with open('science_names2014.tre','w')as w:
        lines = f.readlines()
        for l in lines:
            for m in mapping:
                l = l.replace(m,mapping[m])
                l = l.replace("'",'')
            w.write(l)

            

                    
            

        