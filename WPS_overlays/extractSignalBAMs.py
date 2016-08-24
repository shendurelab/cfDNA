#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *17.02.2016
"""

import sys, os
from optparse import OptionParser
import pysam
import random
from collections import defaultdict 

import gzip
from bx.intervals.intersection import Intersecter, Interval
import random

def isSoftClipped(cigar):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch
  for (op,count) in cigar:
    if op in [4,5,6]: return True
  return False

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def readIterator(filenames,chrom,start,end):
  for bamfile in filenames:
    if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
      input_file = pysam.Samfile( bamfile, "rb" )
      for read in input_file.fetch(chrom,start-1,end):
        yield read
      input_file.close()

def calculate_score(chrom, start, end,ctype="WPS"):
  filteredReads = Intersecter()
  posRange = defaultdict(int)
  for read in readIterator(args,chrom,start,end):
    if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
    if isSoftClipped(read.cigar): continue
    
    if read.is_paired:
      if read.mate_is_unmapped: continue
      if read.rnext != read.tid: continue
      if read.is_read1 or (read.is_read2 and read.pnext+read.qlen < start):
        if read.isize == 0: continue
        if options.downsample != None and random.random() >= options.downsample: continue
        rstart = min(read.pos,read.pnext)+1 # 1-based
        rend = rstart+abs(read.isize)-1 # end included
        rlength = rend-rstart+1
        if options.minLength <= rlength <= options.maxLength:
          filteredReads.add_interval(Interval(rstart,rend))
          if ctype == "COV":
            for i in range(rstart,rend+1):
              if i >= start and i <= end:
                posRange[i]+=1
          elif ctype == "STARTS":
            if rstart >= start and rstart <= end:
              posRange[rstart]+=1
            if rend >= start and rend <= end:
              posRange[rend]+=1
    else:
      if options.downsample != None and random.random() >= options.downsample: continue
      rstart = read.pos+1 # 1-based
      rend = rstart+aln_length(read.cigar)-1 # end included
      rlength = rend-rstart+1
      if options.minLength <= rlength <= options.maxLength:
        filteredReads.add_interval(Interval(rstart,rend))
        if ctype == "COV":
          for i in range(rstart,rend+1):
            if i >= start and i <= end:
              posRange[i]+=1
        elif ctype == "STARTS":
          if rstart >= start and rstart <= end:
            posRange[rstart]+=1
          if rend >= start and rend <= end:
            posRange[rend]+=1

  if ctype == "WPS":
    protection = options.protection//2
    for pos in xrange(start,end+1):
      rstart,rend = pos-protection,pos+protection
      gcount,bcount = 0,0
      for read in filteredReads.find(rstart,rend):
        if (read.start > rstart) or (read.end < rend): bcount +=1
        else: gcount +=1
      posRange[pos]+=gcount-bcount

  res = []
  for pos in xrange(start,end+1):
    res.append(posRange[pos])
  return res


parser = OptionParser()
parser.add_option("-t", "--inputType", dest="inputType", help="Input coordinate type (def: offset, alt: sites,named,nameFirst)",default="offset")
parser.add_option("-o", "--outprefix", dest="outprefix", help="Output file prefix (def ./)",default="./")
parser.add_option("-w", "--window", dest="window", help="Window to extract around central position (def 5000)",default=5000,type="int")
parser.add_option("-r", "--random", dest="random", help="Number of random bins to split data into (def 10)",default=10,type="int")
parser.add_option("", "--name", dest="name", help="Overwrite name obtained from annotation (def None)",default=None)
parser.add_option("-l","--lengthSR", dest="lengthSR", help="Length of full reads (default 76)",default=76,type="int")
parser.add_option("","--merged", dest="merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
parser.add_option("","--trimmed", dest="trimmed", help="Assume reads are trimmed (default Off)",default=False,action="store_true")
parser.add_option("-p","--protection", dest="protection", help="Base pair protection assumed for elements (default 120)",default=120,type="int")
parser.add_option("","--minLength", dest="minLength", help="Minimum fragment length to include (default 120)",default=120,type="int")
parser.add_option("","--maxLength", dest="maxLength", help="Minimum fragment length to include (default 180)",default=180,type="int")
parser.add_option("--downsample", dest="downsample", help="Ratio to down sample reads (default OFF)",default=None,type="float")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-s","--scoreType", dest="scoreType", help="Score type (i.e WPS, COV, STARTS, default WPS)",default="WPS",action="store_true")
(options, args) = parser.parse_args()

if options.scoreType not in ["WPS","COV","STARTS"]: 
  sys.stderr.write("Error: undefined type for score calculation!\n")
  sys.exit()

outprefix = options.outprefix #.rstrip("/")

def initDicts(name,count,length):
  global sites
  sites[name]=[]
  for i in range(count):
    sites[name].append([0.0]*length)

sites = {}

window = options.window
halfwindow = options.window//2

for line in sys.stdin:
  if line.startswith("#"): continue
  chrom,center,strand,name = None,None,None,None
  if options.inputType == "offset":
    chrom,start,end,name,offset,strand = line.rstrip().split("\t")
    start = int(start)
    end = int(end)
    offset = int(offset)
    center = start+offset
  elif options.inputType == "sites":
    chrom,pos,strand = line.rstrip().split(":")
    center = int(pos)
  elif options.inputType == "named":
    chrom,start,end,name,strand = line.rstrip().split(":")
    start = int(start)
    end = int(end)
    center = start+(end-start)//2
  elif options.inputType == "nameFirst":
    name,chrom,start,end,strand = line.rstrip().split(":")
    start = int(start)
    end = int(end)
    center = start+(end-start)//2
  else:
    sys.stderr.write("Unknown file type specificed\n")
    sys.exit()  
  
  if options.name != None: name = options.name
  selBin = random.randint(0,9)
  res1 = []
  
  res1 = calculate_WPS(chrom, center-halfwindow, center+halfwindow-1)
  if (len(res1) != 2*halfwindow): 
    sys.stderr.write("Length does not match %d != %d\n"%(len(res1),2*halfwindow))
    continue
  
  if (strand == "-"): res1=res1[::-1]
  if name not in sites: initDicts(name,options.random,len(res1))
  
  # Sum signal
  sites[name][selBin] = map(lambda (a,b): a+b,zip(sites[name][selBin],res1))

  
for name,values in sites.iteritems():
  outfile = open("%s%s.tsv"%(outprefix,name.replace("/","-").replace("::","-")),'w')
  outfile.write("#Pos\t%s\n"%("\t".join(map(lambda x:"Sites-%d"%(x+1),range(options.random)))))
  for ind,val in enumerate(values[0]):
    outfile.write("%d\t%s\n"%(ind+1,"\t".join(map(lambda x:"%d"%values[x][ind],range(options.random)))))
  outfile.close()

