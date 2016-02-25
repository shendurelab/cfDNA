#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *17.06.2015
"""

import sys, os
from optparse import OptionParser
import pysam
import random

def initDicts(name,count,length):
  global sites
  sites[name]=[]
  for i in range(count):
    sites[name].append([0.0]*length)

parser = OptionParser("%prog [options]")
parser.add_option("-t", "--inputType", dest="inputType", help="Input coordinate type (def: offset, alt: sites,named,nameFirst)",default="offset")
parser.add_option("-o", "--outprefix", dest="outprefix", help="Output file prefix (def ./)",default="./")
parser.add_option("-w", "--window", dest="window", help="Window to extract around central position (def 5000)",default=5000,type="int")
parser.add_option("-r", "--random", dest="random", help="Number of random bins to split data into (def 10)",default=10,type="int")
parser.add_option("-n", "--name", dest="name", help="Overwrite name obtained from annotation (def None)",default=None)
parser.add_option("-i", "--individual", dest="individual", help="File with normalized score of the individual (def 'normedWPS_CA01.tsv.gz')",default="normedWPS_CA01.tsv.gz")
(options, args) = parser.parse_args()

outprefix = options.outprefix.rstrip("/")

sites = {}

window = options.window
halfwindow = options.window//2

individual = pysam.Tabixfile(options.individual)
chromSeg = set(individual.contigs)

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
  if chrom not in chromSeg: continue
  selBin = random.randint(0,9)
  res1 = []
  
  try:
    res1 = map(lambda x:float(x.split("\t")[4]),individual.fetch(chrom, center-halfwindow, center+halfwindow))
  except:
    continue
  
  #print "\n".join(map(lambda x:x,individual.fetch(chrom, center-halfwindow, center+halfwindow)))
  #print res1
  #print len(res1), 2*halfwindow
  #print
  if (len(res1) != 2*halfwindow): continue
  
  if (strand == "-"): res1=res1[::-1]
  if name not in sites: initDicts(name,options.random,len(res1))
  
  # Only sum positive signal
  #sites[name][selBin] = map(lambda (a,b): a+b if b > 0 else a,zip(sites[name][selBin],res1))

  # Sum signal
  sites[name][selBin] = map(lambda (a,b): a+b,zip(sites[name][selBin],res1))
  
for name,values in sites.iteritems():
  outfile = open("%s%s.tsv"%(outprefix,name.replace("/","-").replace("::","-")),'w')
  outfile.write("#Pos\t%s\n"%("\t".join(map(lambda x:"Sites-%d"%(x+1),range(options.random)))))
  for ind,val in enumerate(values[0]):
    outfile.write("%d\t%s\n"%(ind+1,"\t".join(map(lambda x:"%d"%values[x][ind],range(options.random)))))
  outfile.close()
