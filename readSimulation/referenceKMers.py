#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *12.04.2015
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 
import mmap
import gzip

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

def get_fasta_mmap(fastafile, length=0):
  cmap = None
  if os.path.exists(fastafile):
    #OPEN FASTA WITH MMAP
    f = open(fastafile, "r")
    cmap = mmap.mmap(f.fileno(), length=length, access=mmap.ACCESS_READ)
  else:
    sys.stderr.write('Error: Fasta file not available: %s\n'%fastafile)
    sys.exit()
  return cmap

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

def valid(sequence):
  check = True
  for b in sequence:
    if b == "A" or b == "G" or b == "T" or b == "C": continue
    else:
      check = False
      break
  return check

parser = OptionParser()
parser.add_option("-k","--kmer", dest="kmerSize", help="Size of kmers to count (def 8)",default=8,type="int")
parser.add_option("-o","--outfile", dest="outfile", help="Name of output file (def output.tsv)",default="output.tsv")
parser.add_option("-r","--region", dest="region", help="Region to be looked up, empty string for all regular chromosomes (def chr12:34439500-34470500)",default="chr12:34439500-34470500")
parser.add_option("-f","--fasta", dest="reference", help="Fasta index reference genome (default /net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa)",default="/net/shendure/vol1/home/mkircher/sequencedb/genome/grch37_1000g_phase2/whole_genome.fa")
parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
(options, args) = parser.parse_args()

#################
# PARSE REGION
#################

chromosomes = []
options.region = options.region.replace('"','').replace("'","").strip()
if len(options.region) > 0:
  chrom,start,end = None,None,None
  outchrom = None
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
  except:
    sys.stderr.write("Error: Region not defined in a correct format!\n")
    sys.exit()
  chromosomes.append((outchrom,start,end))

if options.verbose:
  sys.stderr.write("Opening genome file...\n")

genome_index = read_genome_fastaindex(options.reference+".fai")
genome_map = get_fasta_mmap(options.reference)
#print genome_index

if len(chromosomes) == 0:
  for chrom,(length,pos,lenLine1,lenLine2) in genome_index.iteritems():
    if chrom.startswith("GL") or chrom.startswith("NC_") or chrom in ["MT","hs37d5"]: continue
    chromosomes.append((chrom,1,length))

kSize = options.kmerSize
chromosomes.sort()
kmerFrequency = defaultdict(int)
for chrom,start,end in chromosomes:
  kmer = ""
  if options.verbose:
    sys.stderr.write("Reading %s:%d-%d...\n"%(chrom,start,end))
  for pos,base in enumerate(DNAiterator(chrom,start,genome_map,genome_index)): 
    if start+pos > end: break
    #Update kmer
    if len(kmer) < kSize:
      kmer=kmer+base.upper()
    else:
      kmer=kmer[1:]+base.upper()
    if len(kmer) == kSize and valid(kmer):
      kmerFrequency[kmer]+=1

if options.verbose:
  sys.stderr.write("Writing output file... Observed %d Kmers\n"%(len(kmerFrequency)))
output = open(options.outfile,'w')
for key,value in sorted(kmerFrequency.iteritems()):
  output.write("%s\t%s\n"%(key,value))
output.close()
