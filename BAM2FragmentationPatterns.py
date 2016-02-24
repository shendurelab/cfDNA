#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *09.01.2015
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

parser = OptionParser()
parser.add_option("-r","--region", dest="region", help="Region to be looked up (def 12:34,443,233-34,453,733)",default="12:34,443,233-34,453,733")
parser.add_option("-f","--fasta", dest="reference", help="Fasta index reference genome (default /net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa)",default="/net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa")
parser.add_option("-m","--max", dest="maxreads", help="Maximum number of reads to consider (default 0 = Off)",default=0,type="int")
parser.add_option("-d","--dinucleotides", dest="dinucleotides", help="Additionaly extract dinucleotide patterns (def Off)",default=False,action="store_true")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
parser.add_option("-o","--outfile", dest="outfile", help="Write output to file (def STDOUT)",default="")
(options, args) = parser.parse_args()

if not os.path.exists(options.reference) or not os.path.exists(options.reference+".fai"):
  sys.stderr.write("Fasta indexed reference genome is not available.\n")
  sys.exit()

reference = pysam.Fastafile(options.reference)

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

data = {  True:{'A':defaultdict(int),'C':defaultdict(int),'G':defaultdict(int),'T':defaultdict(int)}, 
         False:{'A':defaultdict(int),'C':defaultdict(int),'G':defaultdict(int),'T':defaultdict(int)} }

dinucleotides = { True:{}, False:{} }
if options.dinucleotides:
  for base1 in ['A','C','G','T']:
    for base2 in ['A','C','G','T']:
      dinucleotides[True][base1+base2]=defaultdict(int)
      dinucleotides[False][base1+base2]=defaultdict(int)


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
    if alignment.is_unmapped: continue
    if alignment.is_duplicate or alignment.is_qcfail: continue
    if alignment.cigar == None: continue
    if isSoftClipped(alignment.cigar): continue
    
    if chrom == None:
      try:
        outchrom = BAMreferences[alignment.tid]
        outchrom = outchrom.replace("chr","")
      except:
        continue
    
    posList = []
    if alignment.is_paired:
      if alignment.mate_is_unmapped: continue
      if alignment.rnext != alignment.tid: continue
      if alignment.is_reverse: 
        posList.append((alignment.pos+alnLength(alignment.cigar),True,alignment.is_read1))
      else: 
        posList.append((alignment.pos,False,alignment.is_read1)) 
    else:
      posList.append((alignment.pos,False,not(alignment.is_reverse)))
      posList.append((alignment.pos+alnLength(alignment.cigar),True,alignment.is_reverse))
    
    readCount+=1
    for (pos,orient,threePrime) in posList:
      try:
        seq = reference.fetch(outchrom,pos-20,pos+21)
      except:
        continue
      
      if orient:
        rseq = seq.translate(table)[::-1]
        for ind,base in enumerate(rseq[1:]):
          try: data[threePrime][base][ind]+=1
          except: pass
          
          if options.dinucleotides:
            try: dinucleotides[threePrime][base+seq[ind+1]][ind]+=1
            except: pass

      else:
        for ind,base in enumerate(seq[:-1]):
          try: data[threePrime][base][ind]+=1
          except: pass

          if options.dinucleotides:
            try: dinucleotides[threePrime][base+seq[ind+1]][ind]+=1
            except: pass

      ## Sanity checks
      #if not orient:
        #print seq,orient,threePrime
        #print 20*" "+alignment.seq[:20]
      #else:
        #print seq,orient,threePrime
        #print alignment.seq[-20:]
        #print rseq
        #print 20*" "+alignment.seq[-20:].translate(table)[::-1]

      #print
      #for base in ['A','C','G','T']:
        #print sorted(data[True][base].iteritems())
      #print
      #for base in ['A','C','G','T']:
        #print sorted(data[False][base].iteritems())

    if (options.maxreads > 0) and (readCount > options.maxreads): 
      sys.stderr.write("Maximum number of reads reached.\n")
      break 
  input_file.close()

#print
#for base in ['A','C','G','T']:
  #print sorted(data[True][base].iteritems())
#print
#for base in ['A','C','G','T']:
  #print sorted(data[False][base].iteritems())


if options.outfile != "":
  outfile = open(options.outfile,"w")
else:
  outfile = sys.stdout

outfile.write("Type\tBase\t"+"\t".join(map(str,range(-20,0)+range(1,21)))+"\n")
for threePrime in [True,False]:
  helper = map(lambda (x,y):y,sorted(data[threePrime]['A'].iteritems()))
  for base in ['C','G','T']:
    for ind,val in data[threePrime][base].iteritems():
      helper[ind]+=val
  
  for base,subDict in sorted(data[threePrime].iteritems()):
    outfile.write(("ThreePrime" if threePrime else "FivePrime")+"\t"+base+"\t"+"\t".join(map(lambda (ind,x): "NA" if helper[ind] == 0 else "%.4f"%(x/float(helper[ind])),sorted(subDict.iteritems())))+"\n")

if options.dinucleotides:
  for threePrime in [True,False]:
    helper = map(lambda (x,y):y,sorted(dinucleotides[threePrime]['AA'].iteritems()))
    for base in dinucleotides[threePrime].keys():
      if base == "AA": continue
      for ind,val in dinucleotides[threePrime][base].iteritems():
        helper[ind]+=val

    for base,subDict in sorted(dinucleotides[threePrime].iteritems()):
      outfile.write(("ThreePrime" if threePrime else "FivePrime")+"\t"+base+"\t"+"\t".join(map(lambda (ind,x): "NA" if helper[ind] == 0 else "%.4f"%(x/float(helper[ind])),sorted(subDict.iteritems())))+"\n")

if options.outfile != "":
  outfile.close()