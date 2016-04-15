# compute peak density in 100kb bins

import numpy as np
import sys
import gzip

# min and max distances to consider
maxdist = 280
mindist = 120

allPeaks = {}

for i in range(1, 23):
    allPeaks["chr%d" % i] = []
allPeaks["chrX"] = []
allPeaks["chrY"] = []

# read in all peaks
fpeaks = gzip.open(sys.argv[1])
for line in fpeaks:
    parts = line.strip().split()
    chrom = parts[0]
    center = parts[6]
    allPeaks[chrom].append(int(center))

# finished reading peaks, now compute bins
lastchrom = ""
i = 0

# length of each reference chromosome
f1 = open("chrNameLength.txt")
for line in f1:
	parts = line.strip().split()
	chrom = "chr" + parts[0]
	if not chrom == lastchrom:
		i = 0
	lastchrom = chrom
	length = int(parts[1])
	if chrom[3:] in map(str, xrange(23)) or chrom == "chrX" or chrom == "chrY":
		start = 0
		stop = 100000
                i = 0
		while stop < length:
                    lastpeak = None
                    mylist = []
                    name = chrom + ":" + str(start) + "-" + str(stop)
                    for peak in allPeaks[chrom][i:]:
                        if peak <= stop:
                            if lastpeak is not None:
                                dist = peak - lastpeak
                                if dist <= maxdist and dist >= mindist:
                                    mylist.append(dist)
                                # else do nothing
                            lastpeak = peak
                            i += 1
                    print "%s\t%d\t%d\t%s\t%f" % (chrom, start, stop, name, np.median(mylist))
                    start += 100000
                    stop += 100000

		# now handle the last bin
		stop = length
                lastpeak = None
                mylist = []
                name = chrom + ":" + str(start) + "-" + str(stop)
		for peak in allPeaks[chrom][i:]:
                    if start <= peak and peak <= stop:
                        if lastpeak is not None:
                            dist = peak - lastpeak
                            if dist <= maxdist and dist >= mindist:
                                mylist.append(dist)
                                # else do nothing
                        lastpeak = peak
                print "%s\t%d\t%d\t%s\t%f" % (chrom, start, stop, name, np.median(mylist))
