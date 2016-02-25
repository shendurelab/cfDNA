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
#parser.add_option("-v","--verbose", dest="verbose", help="Turn debug output on",default=False,action="store_true")
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
    if fields[2] == "stop_codon" and fields[6] == "+":
      print ":".join([fields[0],fields[3],fields[6]])
    elif fields[2] == "stop_codon" and fields[6] == "-":
      print ":".join([fields[0],fields[4],fields[6]])
infile.close()
