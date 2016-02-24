#!/usr/bin/env python

"""
:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *14.04.2015
"""

import sys, os
from optparse import OptionParser
import pysam
import string
from collections import defaultdict

table = string.maketrans('TGCA','ACGT') # COMPLEMENT DNA

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

def alnLength(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def valid(sequence):
  check = True
  for b in sequence:
    if b == "A" or b == "G" or b == "T" or b == "C": continue
    else:
      check = False
      break
  return check

parser = OptionParser()
parser.add_option("--kmerLE", dest="kmerSizeLE", help="Size of kmers to count (def 8)",default=8,type="int")
parser.add_option("--kmerRE", dest="kmerSizeRE", help="Size of kmers to count (def 5)",default=5,type="int")
parser.add_option("-r","--region", dest="region", help="Region to be looked up (def 12:34,443,233-34,453,733)",default="12:34,443,233-34,453,733")
parser.add_option("-m","--max", dest="maxreads", help="Maximum number of reads to consider (default 0 = Off)",default=0,type="int")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("--outfileLE", dest="outfileLE", help="Write output to file (def leftEnd[_f/_r.tsv])",default="leftEnd")
parser.add_option("--outfileRE", dest="outfileRE", help="Write output to file (def rightEnd[_f/_r.tsv])",default="rightEnd")
(options, args) = parser.parse_args()

kmerFrequencyLEp = defaultdict(int)
kmerFrequencyLEm = defaultdict(int)
kmerFrequencyREp = defaultdict(int)
kmerFrequencyREm = defaultdict(int)

chrom,start,end = None,None,None
outchrom = None
options.region = options.region.strip("""\'""").strip()
if options.region.upper() == "ALL": options.region = ""

if len(options.region) != 0:
  try:
    chrom = options.region.split(':')[0]
    start,end = map(int,options.region.split(':')[1].replace(",","").split("-"))
    outchrom = chrom
    outchrom = outchrom.replace("chr","")
    if outchrom.startswith('gi|'): # gi|224589803|ref|NC_000012.11|
      NCfield = outchrom.split("|")[-2]
      if NCfield.startswith("NC"):
        outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))
    if options.verbose: 
      sys.stderr.write("Coordinates parsed: Chrom %s Start %d End %d\n"%(chrom,start,end))
  except:
    sys.stderr.write("Error: Region not defined in a correct format!\n")
    sys.exit()

validChromosomes = set(map(str,range(1,23))+["X","Y"])

maxKmerLen = max(options.kmerSizeLE,options.kmerSizeRE)

readCount = 0
for bamfile in args:
  if not os.path.exists(bamfile): continue
  if (chrom != None) and not os.path.exists(bamfile.replace(".bam",".bai")) and not os.path.exists(bamfile+".bai"): continue
  if options.verbose: sys.stderr.write("Reading %s\n"%bamfile)

  input_file = pysam.Samfile( bamfile, "rb" )
  BAMreferences = dict(enumerate(input_file.references))
  
  if chrom == None:
    BAMiter = input_file
  else:
    BAMiter = input_file.fetch( chrom,start,end )
  
  for alignment in BAMiter:
    if alignment.is_duplicate or alignment.is_qcfail: continue
    if alignment.cigar == None: continue
    if isSoftClipped(alignment.cigar): continue
    if len(alignment.seq) < maxKmerLen: continue
    
    if chrom == None:
      try: 
        outchrom = BAMreferences[alignment.tid]
        outchrom = outchrom.replace("chr","")
      except: continue
    if outchrom not in validChromosomes: continue
    
    if alignment.is_paired:
      if alignment.mate_is_unmapped: continue
      if alignment.rnext != alignment.tid: continue
      
      if alignment.is_reverse: 
        if alignment.is_read1:
          ckmer = alignment.seq[-options.kmerSizeLE:]
          if valid(ckmer): kmerFrequencyLEm[ckmer]+=1
        else:
          ckmer = alignment.seq[-options.kmerSizeRE:]
          if valid(ckmer): kmerFrequencyREm[ckmer]+=1
      else: 
        if alignment.is_read1:
          ckmer = alignment.seq[:options.kmerSizeLE]
          if valid(ckmer): kmerFrequencyLEp[ckmer]+=1
        else:
          ckmer = alignment.seq[:options.kmerSizeRE]
          if valid(ckmer): kmerFrequencyREp[ckmer]+=1
    else:
      if alignment.is_reverse:
        ckmer = alignment.seq[:options.kmerSizeRE]
        if valid(ckmer): kmerFrequencyREp[ckmer]+=1
        ckmer = alignment.seq[-options.kmerSizeLE:]
        if valid(ckmer): kmerFrequencyLEm[ckmer]+=1
      else:
        ckmer = alignment.seq[:options.kmerSizeLE]
        if valid(ckmer): kmerFrequencyLEp[ckmer]+=1
        ckmer = alignment.seq[-options.kmerSizeRE:]
        if valid(ckmer): kmerFrequencyREm[ckmer]+=1
        
    readCount+=1

    if (options.maxreads > 0) and (readCount > options.maxreads): 
      sys.stderr.write("Maximum number of reads reached.\n")
      break 
  input_file.close()

outfile = open(options.outfileLE+"_f.tsv","w")
for key,value in sorted(kmerFrequencyLEp.iteritems()):
  outfile.write("%s\t%s\n"%(key,value))
outfile.close()
outfile = open(options.outfileLE+"_r.tsv","w")
for key,value in sorted(kmerFrequencyLEm.iteritems()):
  outfile.write("%s\t%s\n"%(key,value))
outfile.close()

outfile = open(options.outfileRE+"_f.tsv","w")
for key,value in sorted(kmerFrequencyREp.iteritems()):
  outfile.write("%s\t%s\n"%(key,value))
outfile.close()
outfile = open(options.outfileRE+"_r.tsv","w")
for key,value in sorted(kmerFrequencyREm.iteritems()):
  outfile.write("%s\t%s\n"%(key,value))
outfile.close()
