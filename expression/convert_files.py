#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *05.07.2014

"""

import sys, os
import pysam
import gzip
from optparse import OptionParser
from collections import defaultdict

parser = OptionParser()
parser.add_option("-i","--individual", dest="individual", help="Individual to use (comma-separate, def 'A11')",default="A11")
parser.add_option("-p","--project", dest="project", help="Project name to use (def 'TSS')",default="TSS")
parser.add_option("-t","--temp", dest="temp", help="Temp directory with blocks (def '/net/shendure/vol1/home/mkircher/nucleosomes/expression')",default="/net/shendure/vol1/home/mkircher/nucleosomes/expression")
parser.add_option("-r","--root", dest="root", help="Root directory for outputs (def '/net/shendure/vol1/home/mkircher/nucleosomes/expression')",default="/net/shendure/vol1/home/mkircher/nucleosomes/expression")
parser.add_option("-a","--annotation", dest="annotation", help="Use regions transcript file (def transcriptAnno.tsv)",default="transcriptAnno.tsv")
parser.add_option("-s","--suffix", dest="suffix", help="Suffix for fft_summaries folder (def '')",default="")
(options, args) = parser.parse_args()

rootDir = options.root
tempDir = options.temp

project = options.project
ind = options.individual

outfileCov = gzip.open("%s/%s/fft_summaries%s/fft_%s_cov.tsv.gz"%(rootDir,project,options.suffix,ind),'w')
outfileStarts = gzip.open("%s/%s/fft_summaries%s/fft_%s_starts.tsv.gz"%(rootDir,project,options.suffix,ind),'w')
outfileWPS = gzip.open("%s/%s/fft_summaries%s/fft_%s_WPS.tsv.gz"%(rootDir,project,options.suffix,ind),'w')

geneIDs = set()
if os.path.exists(options.annotation):
  infile = open(options.annotation)
  for line in infile:
    cid,chrom,start,end,strand = line.split() # positions are 1-based and inclusive
    geneIDs.add(cid)
  infile.close()

wrote_header = False
for cid in sorted(geneIDs):
  if os.path.exists("%s/%s/%s/fft%s/block_%s.tsv.gz"%(tempDir,project,ind,options.suffix,cid)):
    if not wrote_header:
      header = ['#Region']
      infile = gzip.open("%s/%s/%s/fft%s/block_%s.tsv.gz"%(tempDir,project,ind,options.suffix,cid))
      infile.readline()
      for line in infile:
        fields = line.split()
        header.append(fields[0])
      infile.close()
      outfileCov.write("\t".join(header)+"\n")
      outfileStarts.write("\t".join(header)+"\n")
      outfileWPS.write("\t".join(header)+"\n")
      wrote_header=True
    infile = gzip.open("%s/%s/%s/fft%s/block_%s.tsv.gz"%(tempDir,project,ind,options.suffix,cid))  
    covs = [cid]
    starts = [cid]
    wps = [cid]
    infile.readline()
    for line in infile:
      fields = line.split()
      covs.append(fields[1])
      starts.append(fields[2])
      wps.append(fields[3])
    infile.close()
    outfileCov.write("\t".join(covs)+"\n")
    outfileStarts.write("\t".join(starts)+"\n")
    outfileWPS.write("\t".join(wps)+"\n")

outfileCov.close()
outfileStarts.close()
outfileWPS.close()