import numpy as np
import sys
import gzip

allS = open(sys.argv[1])

allPeaks = {}
allScores = {}
scores = []

for line in allS:
    parts = line.strip().split()
    name = parts[0] + ":" + parts[1]
    score = float(parts[3])
    scores.append(score)
    allScores[name] = score
allS.close()

pcs = []
for x in range(1, 101, 1):
    pcs.append(np.percentile(scores, x))

allS = open(sys.argv[1])
for line in allS:
    parts = line.strip().split()
    name = parts[0] + ":" + parts[1]
    v = allScores[name]
    mid = int(parts[1]) + ((int(parts[2]) - int(parts[1])) / 2)
    i = 0
    while v >= pcs[i]:
        if i == len(pcs) - 1:
            break
        i += 1

    print "chr"+parts[0] + "\t" + parts[1] +"\t"+ parts[2] + "\t" + name + "\t" + parts[3] + "\t+\t" + str(mid) + "\t" + str(i)
