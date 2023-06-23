# python3 resampleLines.py [input file path] [output file path] [output lines#]

import sys
import random

inFilePath = sys.argv[1]
outFilePath = sys.argv[2]
if len(sys.argv) > 3:
    outLines = int(sys.argv[3])
else:
    outLines = 0

with open(inFilePath, "r") as inFile:
    data = inFile.read().split("\n")
    if len(data[-1]) == 0:
        data.pop(-1)
    if outLines == 0:
        outLines = len(data)

    resampled = []
    for i in range(outLines):
        resampled.append(random.choice(data))

    with open(outFilePath, "w") as outFile:
        outFile.write("\n".join(resampled))
