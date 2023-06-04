#USAGE:: `python workflow/scripts/plot_dists.py [converge dir] [out_file]`
import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
converge_dir = sys.argv[1]
iter = pd.read_csv(sys.argv[1]+"/iter_dist.csv",header=None)
iter_bs = pd.read_csv(sys.argv[1]+"/iter_dist_bs.csv",header=None)
ref = pd.read_csv(sys.argv[1]+"/ref_dist.csv",header=None)
iter = iter.rename(columns={0:'Run #',1:'iter_dist'})
iter_bs = iter_bs.rename(columns={0:'Run #',1:'bs_dist'})
ref = ref.rename(columns={0:'Run #',1:'ref_dist'})
iter = iter.merge(iter_bs,on="Run #")
combined = iter.merge(ref,on='Run #')
print(combined)
title = input("Please enter plot title: ")
ax = sns.lineplot(data=combined[['iter_dist','bs_dist','ref_dist']])
ax.set_title(title)
plt.tight_layout()
fig = ax.get_figure()
fig.savefig(sys.argv[2])
