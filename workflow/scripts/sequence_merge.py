#sequence_merge.py concats all the genes together and gathers stats
#REQUIRES: Seaborn and Matplotlib
#USAGE: `python workflow/scripts/sequence_merge.py [gene directory] [output file] [plot directory]`
import shutil
import glob,os,sys
import seaborn as sns 
import matplotlib.pyplot as plt
directory = sys.argv[1]
output = sys.argv[2]
plotdir = sys.argv[3]
#stats = sys.argv[3]
species_count = {}
with open(output, 'wb') as outfile:
    for filename in glob.glob(os.path.join(directory,'*.fa')):
        if filename == output:
            # don't want to copy the output into the output
            continue
        f = open(filename,'r')
        t = f.read()
        num = t.count('>')
        s = filename.split('/')
        #print(s[len(s)-1])
        n = s[len(s)-1].split('_')
        name = n[0]+'_'+n[1]
        species_count[name] = num
        with open(filename, 'rb') as readfile:
            shutil.copyfileobj(readfile, outfile)
print(species_count)
#with open(stats,'w') as stat:
 #   for species in species_count:
  #      print(species,species_count[species])
   #     stat.write(species+': '+str(species_count[species])+'\n')
x = list(species_count.keys())
y = list(species_count.values())
ax = sns.barplot(x=x,y=y)
ax.set_xticklabels(ax.get_xticklabels(),rotation=90,ha='right')
ax.set_title('Number of Reference Genes from Each Genome')
plt.tight_layout()
fig = ax.get_figure()
fig.savefig(plotdir)