import sys
import gzip

allS = gzip.open(sys.argv[2])
allPeaks = {}
otherPeaks = {}
otherPercs = {}
otherIntervals = {}

for i in range(1, 23):
    allPeaks["chr%d" % i] = []
    otherPeaks["chr%d" % i] = []

# read in all autosomal peaks
for line in allS:
    parts = line.strip().split()
    chrom = parts[0]
    if chrom == "chrX": continue
    if chrom == "chrY": continue
    if chrom == "X": continue
    if chrom == "Y": continue
    center = parts[6]
    allPeaks[chrom].append(int(center))

# read in DHSes
otherS = open(sys.argv[1])
sampleName= sys.argv[1].split("/")[-1]
sampleName = sampleName.split(".")[0]

for line in otherS:
    parts = line.strip().split()
    chrom = parts[0]
    name = chrom + ":" + parts[1]
    # name = parts[3]
    if chrom == "chrX": continue
    if chrom == "chrY": continue
    if chrom == "X": continue
    if chrom == "Y": continue

    center = int(parts[1]) + ((int(parts[2]) - int (parts[1])) / 2)
    otherPeaks[chrom].append(int(center))
    otherPercs[chrom+":"+str(center)] = parts[7]
    otherIntervals[chrom+ ":"+str(center)] = chrom+":"+parts[1]+":"+parts[2]
otherS.close()

for chrom, v in otherPeaks.items():
    otherPeaks[chrom] = sorted(v)

# find the minimum distance to the peak upstream of the center, 
# and the peak downstream of the center, and then find the 
# distance between those two peaks:  ...^...|.........^...
#                                       ==============  this distance

for count in range(1, 23):
    chrom = "chr%d" % count
    otheri = 0
    allj = 0
    counter = 0

    while otheri < len(otherPeaks[chrom]):
        minneg = 999999
        minpos = 999999
        counter += 1
        if allj >= len(allPeaks[chrom]):
            allj -= 1
        mymin = otherPeaks[chrom][otheri] - allPeaks[chrom][allj]
        if mymin < 0: minneg = abs(mymin)
        elif mymin > 0: minpos = mymin
        else: # mymin == 0
            # not sure what the "right" way to handle this is
            minpos = 0

        nextj = allj
        
        # now try decreasing j
        origj = allj
        allj -= 1
        while allj >= 0:
            newdist = otherPeaks[chrom][otheri] - allPeaks[chrom][allj]
            # if distance is negative, keep decreasing j
            if newdist < 0:
                if abs(newdist) < abs(minneg):
                    minneg = abs(newdist)
                allj -= 1
                continue
            
            # otherwise, distance is positive
            if newdist < minpos:
                minpos = newdist
                nextj = allj
                allj -=1
            else:
                # stop
                break
            
        # now look at negative distances
        allj = origj + 1
        while allj < len(allPeaks[chrom]):
            newdist = otherPeaks[chrom][otheri] - allPeaks[chrom][allj]
            if newdist > 0: # gotta keep going
                if newdist < minpos:
                    minpos = newdist
                allj += 1
                continue
            
            # otherwise, distance is negative
            if abs(newdist) < abs(minneg):
                minneg = abs(newdist)
                nextj = allj
                allj += 1

            else:
                # the distance is increasing
                break

        # okay, now calculate the distance between the two
        totaldistance = minpos + minneg - 1
        
                
        print "%s\t%s\t%s\t%d\t%d\t%d\t%s" % (chrom, otherIntervals[chrom+":"+str(otherPeaks[chrom][otheri])].split(":")[1], otherIntervals[chrom+":"+str(otherPeaks[chrom][otheri])].split(":")[2], minneg, minpos, totaldistance, otherPercs[chrom+":"+str(otherPeaks[chrom][otheri])])
        allj = nextj
        otheri += 1

