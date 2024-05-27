# sequence_merge.py concats all the genes together and gathers stats
# REQUIRES: Seaborn and Matplotlib
# USAGE: `python workflow/scripts/sequence_merge.py [gene directory] [output file] [plot directory]`
import shutil
import glob, os, sys
import seaborn as sns
import matplotlib.pyplot as plt

directory = sys.argv[1]
output = sys.argv[2]
plotdir = sys.argv[3]
statdir = sys.argv[4]
species_count = {}

with open(output, "wb") as outfile:
    for filename in glob.glob(os.path.join(directory, "*.fa")):
        if filename == output or "temp" not in filename:
            # don't want to copy the output into the output
            continue
        f = open(filename, "r")
        t = f.read()
        num = t.count(">")
        # s = filename.split("/")
        # # print(s[len(s)-1])
        # n = s[len(s) - 1].split("_")
        # name = n[0] + "_" + n[1]
        s = filename.split("/")
        # Extracts the file name part from the path
        file_name = s[-1]
        name = file_name[:-len("_temp.fa")]
        species_count[name] = num
        with open(filename, "rb") as readfile:
            shutil.copyfileobj(readfile, outfile)

x = list(species_count.keys())
y = list(species_count.values())

with open(statdir + "/sampling.csv", "w") as f:
    f.write("Species name,number of genes sampled"+"\n")
    for i in range(len(x)):
        f.write(x[i] + "," + str(y[i]) + "\n")
        
plt.figure(figsize=(40, 24))
ax = sns.barplot(x=x, y=y)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
ax.set_title("Number of genes sampled from each input genome",fontsize=18)
ax.set_ylabel("Number of sampled genes", fontsize=14)
plt.tight_layout()
fig = ax.get_figure()
fig.savefig(plotdir,dpi=300)
