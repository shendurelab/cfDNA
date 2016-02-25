#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *04.06.2015
"""

import sys, os
from optparse import OptionParser
import gzip

parser = OptionParser()
parser.add_option("--min", dest="minLength", help="Minimum length of exons to be reported (default 0)",default=0,type="int")
parser.add_option("--max", dest="maxLength", help="Maximum length of exons to be reported (default -1)",default=-1,type="int")
parser.add_option("--geneIDs", dest="geneIDs", help="List of gene IDs (default '' )",default="")
parser.add_option("--transcriptIDs", dest="transcriptIDs", help="List of transcript IDs (default ensemblv75_canonicalTranscriptIDs.protein_coding.lst.gz)",default="ensemblv75_canonicalTranscriptIDs.protein_coding.lst.gz")
(options, args) = parser.parse_args()

IDlist = options.transcriptIDs
transcriptIDs = set()
infile = gzip.open(IDlist)
for line in infile:
  transcriptIDs.add(line.strip())
infile.close()

filterIDs = set()
if (options.geneIDs != "") and os.path.exists(options.geneIDs):
  infile = open(options.geneIDs)
  for line in infile:
    filterIDs.add(line.strip())
  infile.close()

infile = gzip.open("Homo_sapiens.GRCh37.75.gtf.gz")
for line in infile:
  if line.startswith("#"): continue
  fields = line.rstrip().split()
  lfield = ""
  tID = None
  gID = None
  for field in fields:
    if lfield == "gene_id":
      gID = field.rstrip(";").strip('"')
    if lfield == "transcript_id":
      tID = field.rstrip(";").strip('"')
      break
    lfield = field
  if tID in transcriptIDs and ((len(filterIDs) == 0) or (gID in filterIDs)):
    exonSize = int(fields[4])-int(fields[3])+1
    if exonSize < options.minLength: continue
    if options.maxLength > 0 and exonSize > options.maxLength: continue
    
    if fields[2] == "exon" and fields[6] == "+":
      print ":".join([fields[0],fields[3],fields[6]])
    elif fields[2] == "exon" and fields[6] == "-":
      print ":".join([fields[0],fields[4],fields[6]])
infile.close()

