import seaborn as sns
import matplotlib.pyplot as plt
import glob, os, sys
import pandas as pd
import argparse

# benchmark file path
parser = argparse.ArgumentParser(
    description="Process benchmark files and plot results."
)
parser.add_argument("--path", help="Directory containing benchmark files.")
parser.add_argument("--out_dir", help="Output directory for plots and CSV files.")
args = parser.parse_args()

path = args.path
out_dir = args.out_dir

sampling = []
alignment = []
msa = []
tree_build = []
step_names = ["sampling", "lastz", "pasta", "raxml"]
steps = [sampling, alignment, msa, tree_build]
extensions = ["*.sample.txt", "*.lastz.txt", "*.pasta.txt", "*.raxml.txt"]

for i in range(len(extensions)):
    for filename in glob.glob(os.path.join(path, extensions[i])):
        with open(os.path.join(os.getcwd(), filename), "r") as f:
            lines = f.readlines()
            s = lines[1].strip().split()
            seconds = float(s[0]) / 60
            s2 = filename.split("/")
            fn = s2[len(s2) - 1].replace(extensions[i][1:], "")
            cpu_time = float(s[len(s) - 1]) / 60
            steps[i].append((fn, seconds, cpu_time))

step_totals = []
cpu_totals = []

for i in range(len(steps)):
    sec_total = 0
    cpu_total = 0
    for j in range(len(steps[i])):
        sec_total += steps[i][j][1]
        cpu_total += steps[i][j][2]
    step_totals.append((step_names[i], sec_total, cpu_total))

tops = []

for i in range(len(steps)):
    if len(steps[i]) > 50:
        top = sorted(steps[i], key=lambda x: x[1], reverse=True)[:50]
        print(top)
        tops.append(top)
    else:
        tops.append(steps[i])

with open(out_dir + "/step_avg.csv", "w") as w:
    for s in step_totals:
        w.write(s[0] + "," + str(s[1]) + "," + str(s[2]) + "\n")
# with open("astral.txt",'r') as f:
# lines = f.readlines()
# s = lines[1].strip().split()
# seconds = float(s[0])
# cpu = float(s[len(s)-1])
# step_totals.append(("ASTRAL",seconds,cpu))

for i in range(len(tops)):
    print(step_names[i])
    fig, ax = plt.subplots(figsize=(20, 10))
    df = pd.DataFrame(
        {
            "job": [x[0] for x in tops[i]],
            "minutes": [x[1] for x in tops[i]],
            "cpu-time": [x[2] for x in tops[i]],
        }
    )
    print(df.head())
    m = df.melt(id_vars="job")
    print(m.head())
    sns.barplot(data=m, x="job", y="value", hue="variable", ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
    ax.set_title(step_names[i])
    plt.tight_layout()
    fig = ax.get_figure()
    fig.savefig(out_dir + "/" + step_names[i] + ".png")

fig, ax = plt.subplots(figsize=(20, 10))

ax = sns.barplot(x=[k[0] for k in step_totals], y=[k[1] for k in step_totals])
ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
ax.set_title(step_names[i])

plt.tight_layout()
fig = ax.get_figure()
fig.savefig(out_dir + "/step_totals.png")
