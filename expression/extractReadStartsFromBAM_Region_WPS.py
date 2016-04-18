#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *03.06.2014
"""

import sys, os
from optparse import OptionParser
import gzip
import pysam
import random

from collections import defaultdict
from bx.intervals.intersection import Intersecter, Interval

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

parser = OptionParser()
parser.add_option("-i","--input", dest="input", help="Use regions transcript file (def transcriptAnno.tsv)",default="transcriptAnno.tsv")
parser.add_option("-l","--lengthSR", dest="lengthSR", help="Length of full reads (default 76)",default=76,type="int")
parser.add_option("-m","--merged", dest="merged", help="Assume reads are merged (default Off)",default=False,action="store_true")
parser.add_option("-t","--trimmed", dest="trimmed", help="Assume reads are trimmed (default Off)",default=False,action="store_true")
parser.add_option("-w","--protection", dest="protection", help="Base pair protection assumed for elements (default 120)",default=120,type="int")
parser.add_option("-o","--outfile", dest="outfile", help="Outfile prefix (def 'block_%s.tsv.gz')",default='block_%s.tsv.gz') # reserve atleast 6 digits 
parser.add_option("-e","--empty", dest="empty", help="Keep files of empty blocks (def Off)",default=False,action="store_true")
parser.add_option("--minInsert", dest="minInsSize", help="Minimum read length threshold to consider (def None)",default=-1,type="int")
parser.add_option("--maxInsert", dest="maxInsSize", help="Minimum read length threshold to consider (def None)",default=-1,type="int")
parser.add_option("--max_length", dest="max_length", help="Assumed maximum insert size (default 1000)",default=1000,type="int")
parser.add_option("--downsample", dest="downsample", help="Ratio to down sample reads (default OFF)",default=None,type="float")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
(options, args) = parser.parse_args()

minInsSize,maxInsSize = None,None
if options.minInsSize > 0 and options.maxInsSize > 0 and options.minInsSize < options.maxInsSize:
  minInsSize = options.minInsSize
  maxInsSize = options.maxInsSize
  sys.stderr.write("Using min/max length cutoffs: %d/%d\n"%(minInsSize,maxInsSize))
  
options.outfile = options.outfile.strip("""\'""")

protection = options.protection//2

validChroms = set(map(str,range(1,23)+["X","Y"]))

if os.path.exists(options.input):
  infile = open(options.input)
  for line in infile:
    cid,chrom,start,end,strand = line.split() # positions should be 0-based and end non-inclusive
    if chrom not in validChroms: continue
    
    regionStart,regionEnd = int(start),int(end)
    
    if regionStart < 1: continue
    
    #outchrom = chrom.replace("chr","")
    #if outchrom.startswith('gi|'):
      #NCfield = outchrom.split("|")[-2]
      #if NCfield.startswith("NC"):
        #outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))

    posRange = defaultdict(lambda:[0,0])
    filteredReads = Intersecter()

    for bamfile in args:
      if options.verbose: sys.stderr.write("Reading %s\n"%bamfile)
      bamfile = bamfile.strip("""\'""")
      if os.path.exists(bamfile) and (os.path.exists(bamfile.replace(".bam",".bai")) or os.path.exists(bamfile+".bai")):
        input_file = pysam.Samfile( bamfile, "rb" )
        prefix = ""
        for tchrom in input_file.references:
          if tchrom.startswith("chr"): 
            prefix = "chr"
            break

        for read in input_file.fetch(prefix+chrom,regionStart-protection-1,regionEnd+protection+1):
          if read.is_duplicate or read.is_qcfail or read.is_unmapped: continue
          if isSoftClipped(read.cigar): continue
          
          if read.is_paired:
            if read.mate_is_unmapped: continue
            if read.rnext != read.tid: continue
            if read.is_read1 or (read.is_read2 and read.pnext+read.qlen < regionStart-protection-1):
              if read.isize == 0: continue
              if options.downsample != None and random.random() >= options.downsample: continue
              rstart = min(read.pos,read.pnext)+1 # 1-based
              lseq = abs(read.isize)
              rend = rstart+lseq-1 # end included
              if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)): continue
              
              filteredReads.add_interval(Interval(rstart,rend))
              #print read.qname,rstart,rend,rend-rstart,abs(read.isize)
              for i in range(rstart,rend+1):
                if i >= regionStart and i <= regionEnd:
                  posRange[i][0]+=1
              if rstart >= regionStart and rstart <= regionEnd:
                posRange[rstart][1]+=1
              if rend >= regionStart and rend <= regionEnd:
                posRange[rend][1]+=1
          else:
            if options.downsample != None and random.random() >= options.downsample: continue
            rstart = read.pos+1 # 1-based
            lseq = aln_length(read.cigar)
            rend = rstart+lseq-1 # end included
            if minInsSize != None and ((lseq < minInsSize) or (lseq > maxInsSize)): continue
            
            filteredReads.add_interval(Interval(rstart,rend))
            #print read.qname,rstart,rend,rend-rstart,aln_length(read.cigar)
            for i in range(rstart,rend+1):
              if i >= regionStart and i <= regionEnd:
                posRange[i][0]+=1
            if ((options.merged or read.qname.startswith('M_')) or ((options.trimmed or read.qname.startswith('T_')) and read.qlen <= options.lengthSR-10)):
              if (rstart >= regionStart and rstart <= regionEnd):
                posRange[rstart][1]+=1
              if rend >= regionStart and rend <= regionEnd:
                posRange[rend][1]+=1
            elif read.is_reverse:
              if rend >= regionStart and rend <= regionEnd:
                posRange[rend][1]+=1
            else:
              if (rstart >= regionStart and rstart <= regionEnd):
                posRange[rstart][1]+=1

    filename = options.outfile%cid
    outfile = gzip.open(filename,'w')
    cov_sites = 0
    outLines = []
    for pos in range(regionStart,regionEnd+1):
      rstart,rend = pos-protection,pos+protection
      gcount,bcount = 0,0
      for read in filteredReads.find(rstart,rend):
        if (read.start > rstart) or (read.end < rend): bcount +=1
        else: gcount +=1
      covCount,startCount = posRange[pos]
      cov_sites += covCount
      outLines.append("%s\t%d\t%d\t%d\t%d\n"%(chrom,pos,covCount,startCount,gcount-bcount))

    if strand == "-": outLines = outLines[::-1]
    for line in outLines: outfile.write(line)
    outfile.close()

    if cov_sites == 0 and not options.empty:
      os.remove(filename)