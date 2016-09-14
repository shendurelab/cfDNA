#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *12.04.2015
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 
from sortedcontainers import SortedDict
import random
import mmap
import gzip
import pysam

def read_genome_fastaindex(faifile):
  #READ FASTA INDEX TO MEMORY
  fastaindex = {}
  if os.path.exists(faifile):
    infile = open(faifile)
    for line in infile:
      fields = line.split()
      if len(fields) == 5:
        cname,length,start,line,cline = fields
        fastaindex[cname]=int(length),int(start),int(line),int(cline)
      else:
        sys.stderr.write('Error: Unexpected line in fasta index file: %s\n'%(line.strip()))
        sys.exit()
    infile.close()
  else:
    sys.stderr.write('Error: Fasta index file not available: %s\n'%faifile)
    sys.exit()
  return fastaindex

def get_DNA_file(chrom,start,end,fastafile,fastaindex):
  try:
    length,bstart,bline,cline = fastaindex[chrom]
    bases = ""
    for pos in xrange(start,end+1):
      if pos <= length and pos >= 1:
        hpos = pos -1
        npos = (hpos // bline)*cline+(hpos % bline)
        fastafile.seek(bstart+npos)
        bases+=fastafile.read(1)
    return bases
  except:
    return None

def DNAiterator(chrom,start,fastamap,fastaindex):
  try:
    length,bstart,bline,cline = fastaindex[chrom]
    pos = start
    while (pos <= length and pos >= 1):
      hpos = pos -1
      npos = (hpos // bline)*cline+(hpos % bline)
      fastamap.seek(bstart+npos)
      yield fastamap.read(1)
      pos += 1
  except:
    pass
  raise StopIteration

def writeBAMentry(chrom,pos,seq,length,strand):
  global BAMoutfile,rcounter
  global options

  chromTID = BAMoutfile.gettid(chrom)
  if chromTID == -1: return
  
  qual = "I"*len(seq)
  
  if length > options.readlength:
    forward = pysam.AlignedRead()
    forward.qname = "SIM-%s-%09d"%(chrom,rcounter)
    forward.is_paired = True
    forward.is_proper_pair = True
    forward.is_read1 = True
    forward.tid = chromTID
    forward.rnext = chromTID

    reverse = pysam.AlignedRead()
    reverse.qname = "SIM-%s-%09d"%(chrom,rcounter)
    reverse.is_paired = True
    reverse.is_proper_pair = True
    reverse.is_read2 = True
    reverse.tid = chromTID
    reverse.rnext = chromTID
    
    if strand:
      forward.seq = seq[:options.readlength]
      forward.qual = qual[:options.readlength]
      forward.pos = pos-1
      forward.mpos = pos+length-options.readlength-1
      forward.isize = -length
      forward.mate_is_reverse = True

      reverse.seq = seq[-options.readlength:]
      reverse.qual = qual[-options.readlength:]
      reverse.is_reverse = True
      reverse.pos = pos+length-options.readlength-1
      reverse.mpos = pos-1
      reverse.isize = length
    
    else:
      forward.seq = seq[-options.readlength:]
      forward.qual = qual[-options.readlength:]
      forward.is_reverse = True
      forward.pos = pos-options.readlength
      forward.mpos = pos-length
      forward.isize = -length
      
      reverse.seq = seq[:options.readlength]
      reverse.qual = qual[:options.readlength]
      reverse.pos = pos-length
      reverse.mpos = pos-options.readlength
      reverse.isize = length
      reverse.mate_is_reverse = True

    forward.mapq = 255
    forward.cigar = [(0,options.readlength)]
    reverse.mapq = 255
    reverse.cigar = [(0,options.readlength)]
    BAMoutfile.write(forward)
    BAMoutfile.write(reverse)
    rcounter += 1
  else:
    forward = pysam.AlignedRead()
    forward.is_reverse = not strand
    forward.qname = "M_SIM-%s-%09d"%(chrom,rcounter)
    forward.seq = seq
    forward.qual = qual
    forward.tid = chromTID
    if strand: forward.pos = pos-1
    else: forward.pos = pos-length
    forward.mapq = 255
    forward.cigar = [(0,length)]
    forward.mrnm = -1
    forward.mpos = -1
    
    BAMoutfile.write(forward)
    rcounter += 1

  #print chrom,pos,length,seq,len(seq)


parser = OptionParser()
parser.add_option("-o","--outfile", dest="outfile", help="Name of output file (def output.bam)",default="output.bam")
parser.add_option("-r","--region", dest="region", help="Region to be looked up (def chr12:34439500-34470500)",default="chr12:34439500-34470500")
parser.add_option("-f","--fasta", dest="reference", help="Fasta index reference genome (default /net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa)",default="/net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa")
parser.add_option("-p","--pipe", dest="pipe", help="Do not fetch, stream data (def Off)",default=False,action="store_true")
parser.add_option("-d","--lengthDist", dest="lengthDist", help="Length distribution (default P12.19.17_Normal_lenDist.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/P12.19.17_Normal_lenDist.tsv")
parser.add_option("--fwdKMerGenome", dest="fwdKMerGenome", help="Frequency distribution of KMers for left reads ends in the genome (default grch37_regChroms_8mers.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/grch37_regChroms_8mers.tsv")
parser.add_option("--revKMerGenome", dest="revKMerGenome", help="Frequency distribution of KMers for right reads ends in the genome (default grch37_regChroms_5mers.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/grch37_regChroms_5mers.tsv")
parser.add_option("--fwdPKMers", dest="fwdPKMers", help="Frequency distribution of KMers at left reads ends on plus strand (default P12.19.17_Normal_left_f.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/P12.19.17_Normal_left_f.tsv")
parser.add_option("--fwdMKMers", dest="fwdMKMers", help="Frequency distribution of KMers at left reads ends on minus strand (default P12.19.17_Normal_left_r.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/P12.19.17_Normal_left_r.tsv")
parser.add_option("--revPKMers", dest="revPKMers", help="Frequency distribution of KMers at right reads ends on plus strand (default P12.19.17_Normal_right_f.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/P12.19.17_Normal_right_f.tsv")
parser.add_option("--revMKMers", dest="revMKMers", help="Frequency distribution of KMers at right reads ends on minus strand (default P12.19.17_Normal_right_r.tsv)",default="/net/shendure/vol1/home/mkircher/nucleosomes/simulation/v2/P12.19.17_Normal_right_r.tsv")
#parser.add_option("--fwdPKMers", dest="fwdPKMers", help="Frequency distribution of KMers at left reads ends on plus strand (default testLeft_f.tsv)",default="testLeft_f.tsv")
#parser.add_option("--fwdMKMers", dest="fwdMKMers", help="Frequency distribution of KMers at left reads ends on minus strand (default testLeft_r.tsv)",default="testLeft_r.tsv")
#parser.add_option("--revPKMers", dest="revPKMers", help="Frequency distribution of KMers at right reads ends on plus strand (default testRight_f.tsv)",default="testRight_f.tsv")
#parser.add_option("--revMKMers", dest="revMKMers", help="Frequency distribution of KMers at right reads ends on minus strand (default testRight_r.tsv)",default="testRight_r.tsv")
parser.add_option("-l", "--readlength", dest="readlength", help="Read length (default 45)", default=45, type="int")
parser.add_option("-s", "--sample", dest="sample", help="Number of times to sample from a position (default 30)", default=30, type="int")
parser.add_option("-c", "--correction", dest="correction", help="Likelihood correction factor to increase chances of piking a putative alignment (default 1.0)", default=1.0, type="float")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
(options, args) = parser.parse_args()

#################
# PARSE REGION
#################

chromosomes = []
options.region = options.region.replace('"',"").replace("'","").strip()

if len(options.region) > 0:
  try:
    chrom = options.region.split(':')[0]
    start,end = map(int,options.region.split(':')[1].replace(",","").split("-"))
    outchrom = chrom
    outchrom = outchrom.replace("chr","")
    if outchrom.startswith('gi|'): # e.g. gi|224589803|ref|NC_000012.11|
      NCfield = outchrom.split("|")[-2]
      if NCfield.startswith("NC"):
        outchrom = "%d"%(int(NCfield.split(".")[0].split("_")[-1]))
    if options.verbose: 
      sys.stderr.write("Coordinates parsed: Chrom %s Start %d End %d\n"%(chrom,start,end))
    chromosomes.append((outchrom,start,end))
  except:
    sys.stderr.write("Error: Region not defined in a correct format!\n")
    sys.exit()

#############################
# READ LENGTH DISTRIBUTION
#############################

lengthDist = SortedDict(23) # load -> math.sqrt(500)
if os.path.exists(options.lengthDist):
  infile = open(options.lengthDist)
  line = infile.readline() # Skip header
  total = 0
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      length,count = map(int,fields)
      total += count
  infile.close()
  infile = open(options.lengthDist)
  line = infile.readline() # Skip header
  total = float(total)
  rsum = 0
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      length,count = map(int,fields)
      rsum+=count
      lengthDist[rsum/total] = length
  infile.close()
else:
  sys.stderr.write("Error: Could not open length distribution!\n")
  sys.exit()

#print lengthDist
#sel = random.random()
#key = lengthDist.bisect_left(sel)
#print sel,key,lengthDist.iloc[key]

minLength = lengthDist[lengthDist.iloc[0]]
maxLength = lengthDist[lengthDist.iloc[-1]]

#############################
# READ KMER DISTRIBUTIONS
#############################

forwardGenome = defaultdict(int)
reverseGenome = defaultdict(int)
forwardLength = None
reverseLength = None

if os.path.exists(options.fwdKMerGenome):
  infile = open(options.fwdKMerGenome)
  total = 0 
  minCount,maxCount = None,None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
      if forwardLength == None:
        forwardLength = len(seq)
      total += count
  infile.close()
  if options.verbose: print "total/min/max:",total,minCount,maxCount
  total = float(maxCount)
  infile = open(options.fwdKMerGenome)
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      forwardGenome[seq] = count/total
  infile.close()
else:
  sys.stderr.write("Error: Could not open left read end kmer distribution of genome!\n")
  sys.exit()

if os.path.exists(options.revKMerGenome):
  infile = open(options.revKMerGenome)
  total = 0 
  minCount,maxCount = None,None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
      if reverseLength == None:
        reverseLength = len(seq)
      total += count
  infile.close()
  if options.verbose: print "total/min/max:",total,minCount,maxCount
  total = float(maxCount)
  infile = open(options.revKMerGenome)
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      reverseGenome[seq] = count/total
  infile.close()
else:
  sys.stderr.write("Error: Could not open right read end kmer distribution of genome!\n")
  sys.exit()

forward = { True:defaultdict(int), False:defaultdict(int) }
reverse = { True:defaultdict(int), False:defaultdict(int) }

if os.path.exists(options.fwdPKMers):
  infile = open(options.fwdPKMers)
  total = 0
  minCount,maxCount = None,None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
  infile.close()
  if options.verbose: print "total/min/max:",total,minCount,maxCount
  total = float(maxCount)
  infile = open(options.fwdPKMers)
  maxCount = None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      freq = (count/total)/forwardGenome[seq]
      if maxCount == None: 
        maxCount = freq
      if maxCount < freq: maxCount = freq
      forward[True][seq] = freq
  infile.close()
  if options.verbose: print maxCount 
  maxCount = float(maxCount)
  for key,value in forward[True].iteritems():
    forward[True][key]=value/maxCount
else:
  sys.stderr.write("Error: Could not open left read end kmer [plus strand] distribution!\n")
  sys.exit()
  
if os.path.exists(options.fwdMKMers):
  infile = open(options.fwdMKMers)
  total = 0 
  minCount,maxCount = None,None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
  infile.close()
  if options.verbose: print "total/min/max:",total,minCount,maxCount
  total = float(maxCount)
  infile = open(options.fwdMKMers)
  maxCount = None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      freq = (count/total)/forwardGenome[seq]
      if maxCount == None: 
        maxCount = freq
      if maxCount < freq: maxCount = freq
      forward[False][seq] = freq
  infile.close()
  if options.verbose: print maxCount 
  maxCount = float(maxCount)
  for key,value in forward[False].iteritems():
    forward[False][key]=value/maxCount
else:
  sys.stderr.write("Error: Could not open left read end kmer [minus strand] distribution!\n")
  sys.exit()

if os.path.exists(options.revPKMers):
  infile = open(options.revPKMers)
  total = 0
  minCount,maxCount = None,None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
  infile.close()
  if options.verbose: print "total/min/max:",total,minCount,maxCount
  total = float(maxCount)
  infile = open(options.revPKMers)
  maxCount = None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      freq = (count/total)/reverseGenome[seq]
      if maxCount == None: 
        maxCount = freq
      if maxCount < freq: maxCount = freq
      reverse[True][seq] = freq
  infile.close()
  if options.verbose: print maxCount 
  maxCount = float(maxCount)
  for key,value in reverse[True].iteritems():
    reverse[True][key]=value/maxCount
else:
  sys.stderr.write("Error: Could not open right read end kmer [plus strand] distribution!\n")
  sys.exit()
  
if os.path.exists(options.revMKMers):
  infile = open(options.revMKMers)
  total = 0 
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      total += count
      if minCount == None:
        minCount = count
      if maxCount == None:
        maxCount = count
      if minCount > count: minCount = count
      if maxCount < count: maxCount = count
  infile.close()
  if options.verbose: print "total/min/max:",total,minCount,maxCount
  total = float(maxCount)
  infile = open(options.revMKMers)
  maxCount = None
  for line in infile:
    fields = line.split()
    if len(fields) == 2:
      seq,count = fields[0],int(fields[-1])
      freq = (count/total)/reverseGenome[seq]
      if maxCount == None: 
        maxCount = freq
      if maxCount < freq: maxCount = freq
      reverse[False][seq] = freq
  infile.close()
  if options.verbose: print maxCount 
  maxCount = float(maxCount)
  for key,value in reverse[False].iteritems():
    reverse[False][key]=value/maxCount
else:
  sys.stderr.write("Error: Could not open right read end kmer [minus strand] distribution!\n")
  sys.exit()

if options.verbose:
  print 
  print "Forward kmer length:",forwardLength
  print "Reverse kmer length:",reverseLength
  print
  print "Initialized kmer sampling rates:"
  print "FOW+",str(forward[True])[:100],sum(forward[True].values())
  print "FOW-",str(forward[False])[:100],sum(forward[False].values())
  print "REV+",str(reverse[True])[:100],sum(reverse[True].values())
  print "REV-",str(reverse[False])[:100],sum(reverse[False].values())

genome_index = read_genome_fastaindex(options.reference+".fai")
genome_map = open(options.reference,'r')
#genome_map = get_fasta_mmap(options.reference)
#genome = pysam.Fastafile(options.reference)

rcounter = 1
BAMoutfile = None
BAMoutfile = pysam.Samfile(options.outfile, "wb", header={'SQ': [{'LN': 249250621, 'SN': '1'}, {'LN': 243199373, 'SN': '2'}, {'LN': 198022430, 'SN': '3'}, {'LN': 191154276, 'SN': '4'}, {'LN': 180915260, 'SN': '5'}, {'LN': 171115067, 'SN': '6'}, {'LN': 159138663, 'SN': '7'}, {'LN': 146364022, 'SN': '8'}, {'LN': 141213431, 'SN': '9'}, {'LN': 135534747, 'SN': '10'}, {'LN': 135006516, 'SN': '11'}, {'LN': 133851895, 'SN': '12'}, {'LN': 115169878, 'SN': '13'}, {'LN': 107349540, 'SN': '14'}, {'LN': 102531392, 'SN': '15'}, {'LN': 90354753, 'SN': '16'}, {'LN': 81195210, 'SN': '17'}, {'LN': 78077248, 'SN': '18'}, {'LN': 59128983, 'SN': '19'}, {'LN': 63025520, 'SN': '20'}, {'LN': 48129895, 'SN': '21'}, {'LN': 51304566, 'SN': '22'}, {'LN': 155270560, 'SN': 'X'}, {'LN': 59373566, 'SN': 'Y'}, {'LN': 16569, 'SN': 'MT'}, {'LN': 4262, 'SN': 'GL000207.1'}, {'LN': 15008, 'SN': 'GL000226.1'}, {'LN': 19913, 'SN': 'GL000229.1'}, {'LN': 27386, 'SN': 'GL000231.1'}, {'LN': 27682, 'SN': 'GL000210.1'}, {'LN': 33824, 'SN': 'GL000239.1'}, {'LN': 34474, 'SN': 'GL000235.1'}, {'LN': 36148, 'SN': 'GL000201.1'}, {'LN': 36422, 'SN': 'GL000247.1'}, {'LN': 36651, 'SN': 'GL000245.1'}, {'LN': 37175, 'SN': 'GL000197.1'}, {'LN': 37498, 'SN': 'GL000203.1'}, {'LN': 38154, 'SN': 'GL000246.1'}, {'LN': 38502, 'SN': 'GL000249.1'}, {'LN': 38914, 'SN': 'GL000196.1'}, {'LN': 39786, 'SN': 'GL000248.1'}, {'LN': 39929, 'SN': 'GL000244.1'}, {'LN': 39939, 'SN': 'GL000238.1'}, {'LN': 40103, 'SN': 'GL000202.1'}, {'LN': 40531, 'SN': 'GL000234.1'}, {'LN': 40652, 'SN': 'GL000232.1'}, {'LN': 41001, 'SN': 'GL000206.1'}, {'LN': 41933, 'SN': 'GL000240.1'}, {'LN': 41934, 'SN': 'GL000236.1'}, {'LN': 42152, 'SN': 'GL000241.1'}, {'LN': 43341, 'SN': 'GL000243.1'}, {'LN': 43523, 'SN': 'GL000242.1'}, {'LN': 43691, 'SN': 'GL000230.1'}, {'LN': 45867, 'SN': 'GL000237.1'}, {'LN': 45941, 'SN': 'GL000233.1'}, {'LN': 81310, 'SN': 'GL000204.1'}, {'LN': 90085, 'SN': 'GL000198.1'}, {'LN': 92689, 'SN': 'GL000208.1'}, {'LN': 106433, 'SN': 'GL000191.1'}, {'LN': 128374, 'SN': 'GL000227.1'}, {'LN': 129120, 'SN': 'GL000228.1'}, {'LN': 137718, 'SN': 'GL000214.1'}, {'LN': 155397, 'SN': 'GL000221.1'}, {'LN': 159169, 'SN': 'GL000209.1'}, {'LN': 161147, 'SN': 'GL000218.1'}, {'LN': 161802, 'SN': 'GL000220.1'}, {'LN': 164239, 'SN': 'GL000213.1'}, {'LN': 166566, 'SN': 'GL000211.1'}, {'LN': 169874, 'SN': 'GL000199.1'}, {'LN': 172149, 'SN': 'GL000217.1'}, {'LN': 172294, 'SN': 'GL000216.1'}, {'LN': 172545, 'SN': 'GL000215.1'}, {'LN': 174588, 'SN': 'GL000205.1'}, {'LN': 179198, 'SN': 'GL000219.1'}, {'LN': 179693, 'SN': 'GL000224.1'}, {'LN': 180455, 'SN': 'GL000223.1'}, {'LN': 182896, 'SN': 'GL000195.1'}, {'LN': 186858, 'SN': 'GL000212.1'}, {'LN': 186861, 'SN': 'GL000222.1'}, {'LN': 187035, 'SN': 'GL000200.1'}, {'LN': 189789, 'SN': 'GL000193.1'}, {'LN': 191469, 'SN': 'GL000194.1'}, {'LN': 211173, 'SN': 'GL000225.1'}, {'LN': 547496, 'SN': 'GL000192.1'}, {'LN': 171823, 'SN': 'NC_007605'}, {'LN': 35477943, 'SN': 'hs37d5'}], 'HD': {'SO': 'unknown', 'VN': '1.4'}})

if len(chromosomes) == 0:
  for chrom,(length,pos,lenLine1,lenLine2) in genome_index.iteritems():
    if chrom.startswith("GL") or chrom.startswith("NC_") or chrom in ["MT","hs37d5"]: continue
    chromosomes.append((chrom,1,length))
chromosomes.sort()

# Empirically identified factor for required iterations to achieve required coverage
targetCov = abs(options.sample)
#if options.correction < 1.0 and abs(int(options.sample*options.correction)) > 1:
 #targetCov = abs(int(options.sample*options.correction)+1)
 #if targetCov < options.sample*options.correction:
   #options.correction = 1.0

correctForward = options.correction
if ((4**forwardLength)/float(4**reverseLength) > options.correction and options.correction > 1):
  correctReverse = 1.0
else:
  correctReverse = options.correction/((4**forwardLength)/float(4**reverseLength))

if options.verbose:
  print "Sampling, FwdCorrect, RevCorrect:",targetCov,correctForward,correctReverse

count = 0
for chrom,start,end in chromosomes:
  if options.verbose:
    sys.stderr.write("Reading %s:%d-%d...\n"%(chrom,start,end))
  context = get_DNA_file(chrom,max(1,start-maxLength),min(genome_index[chrom][0],start+maxLength),genome_map,genome_index)
  posInWindow = maxLength
  if (start-maxLength) < 1:
    posInWindow = maxLength-(maxLength-start)-1
  #print posInWindow,len(context),context[posInWindow-25:posInWindow+25]
  for pos in xrange(start,end+1):
    #print chrom,pos,context[posInWindow],posInWindow,get_DNA_file(chrom,pos,pos,genome_map,genome_index)

    posStrand = []
    # CAN ONLY EXTRACT PLUS STRAND
    if (posInWindow < maxLength) and (posInWindow+maxLength <= len(context)): 
      posStrand = [True]
    # CAN ONLY EXTRACT MINUS STRAND
    elif (posInWindow >= maxLength) and (posInWindow+maxLength > len(context)): 
      posStrand = [False]
    # CAN THEORETICALLY EXTRACT FROM BOTH STRANDS
    else: posStrand = [True,False]
    
    for strand in posStrand:
      kmer = ""
      if strand: kmer = context[posInWindow:posInWindow+forwardLength]
      else: kmer = context[posInWindow-forwardLength+1:posInWindow+1]
      
      for i in xrange(targetCov//len(posStrand)):
        # Pick a threshold
        rval = random.random()
        if rval <= forward[strand][kmer]*correctForward:
          # Pick a length
          lval = random.random()
          key = lengthDist.bisect_left(lval)
          selLength = lengthDist[lengthDist.iloc[key]]
        
          # Check the kmer at the other end
          rkmer = ""
          if strand: rkmer = context[posInWindow+selLength-reverseLength:posInWindow+selLength]
          else: rkmer = context[posInWindow-selLength+1:posInWindow-selLength+reverseLength+1]

          rval = random.random()
          if rval <= reverse[not strand][rkmer]*correctReverse:
            # Extract sequence and write BAM record
            selSeq = ""
            if strand: selSeq = context[posInWindow:posInWindow+selLength]
            else: selSeq = context[posInWindow-selLength+1:posInWindow+1]
            writeBAMentry(chrom,pos,selSeq,selLength,strand)
            #print strand,kmer,rkmer,selSeq,len(selSeq),selLength
      #else:
        #print rval,forward[strand][kmer]
    # ADVANCE IN WINDOW
    nbase = get_DNA_file(chrom,pos+maxLength+1,pos+maxLength+1,genome_map,genome_index)
    if nbase != None:
      if len(context) == 2*maxLength+1:
        context = context[1:]+nbase
      else:
        context += nbase
        posInWindow += 1
    else:
      posInWindow += 1
      
    count += 1
    #if count > 5: sys.exit()
    if options.verbose and (count % 100000 == 0): 
      print "CurPos:", chrom,pos,context[posInWindow],posInWindow,get_DNA_file(chrom,pos+1,pos+1,genome_map,genome_index),"Sim:",rcounter

if BAMoutfile != None:
  BAMoutfile.close()