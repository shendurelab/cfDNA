#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *02.03.2015
"""

import sys, os
import gzip
from optparse import OptionParser

parser = OptionParser()
(options, args) = parser.parse_args()

IDs = "ensemblv75_canonicalTranscriptIDs.lst.gz"
# zcat ensemblv75_canonicalTranscriptIDs.lst.gz | sort | uniq | wc -l
# 64102

GTFfile = "Homo_sapiens.GRCh37.75.gtf.gz"

filterIDs = set()
for line in gzip.open(IDs):
  cid = line.strip()
  filterIDs.add(cid)

sys.stdout.write("\t".join(["#Chrom","Start","End","TranscriptID","Strand","BioType","GeneID","GeneName"])+"\n")
for line in gzip.open(GTFfile):
  if line.startswith("#"): continue
  fields = line.split()
  chrom = fields[0]
  start = fields[3]
  end = fields[4]
  ltype = fields[2]
  strand = fields[6]
  biotype = fields[1]
  geneID = fields[9].replace("\"","").replace(";","")
  transcriptID = fields[11].replace("\"","").replace(";","")
  geneName = fields[13].replace("\"","").replace(";","")
  #print ltype,chrom,start,end,transcriptID,strand,biotype,geneID,geneName
  if (ltype == "transcript") and (transcriptID in filterIDs):
    sys.stdout.write("\t".join([chrom,start,end,transcriptID,strand,biotype,geneID,geneName])+"\n")
