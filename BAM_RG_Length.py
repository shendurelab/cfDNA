#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *13.04.2012
"""

import sys, os
from optparse import OptionParser
from collections import defaultdict 
import pysam

parser = OptionParser("%prog [options]")
parser.add_option("-p","--prefix", dest="prefix", help="Prefix for output filenames (default Length)",default="Length")
parser.add_option("-l","--library", dest="library", help="Use library name from RG read header rather than the sample ID",default=False,action="store_true")
parser.add_option("--max_length", dest="max_length", help="Maximum length considered for the output (default 1000)",default=1000,type="int")
parser.add_option("--noRG", dest="noRG", help="Ignore read group information (output name = output prefix + .tsv )",default=False,action="store_true")
(options, args) = parser.parse_args()

if options.library: options.all=True
have_XP = False

rgroups = {}

for filename in args:
  if os.path.exists(filename):
    print "Reading %s..."%filename
    cbamfile = pysam.Samfile(filename, "rb" )
    id2lib = {}
    if options.library and 'RG' in cbamfile.header:
      for rgroup in cbamfile.header['RG']:
        if 'LB' in rgroup and 'ID' in rgroup:
          id2lib[rgroup['ID']] = rgroup['LB']
    for read in cbamfile:
        library,count = '',1
        for (key,value) in read.tags:
          if key == "RG":
            if value in id2lib: library = id2lib[value]
            else: library = value
          elif key == "XP":
            have_XP = True
            count = value
        if options.noRG: library = ''
        if library not in rgroups: rgroups[library] = [defaultdict(int),defaultdict(int)]
        if not read.is_paired:
          rgroups[library][0][len(read.seq)]+=1
          rgroups[library][1][len(read.seq)]+=count
        elif read.is_read1:
          length = min(options.max_length,abs(read.isize))
          if length == 0: length = options.max_length
          rgroups[library][0][length]+=1
          rgroups[library][1][length]+=count

for library in rgroups:
  if library != '': outfile = open("%s_%s.tsv"%(options.prefix.rstrip("_"),library),'w')
  else: outfile = open("%s.tsv"%(options.prefix),'w')
  if have_XP:
    outfile.write('Length\tCounts\tInclDuplicates\n')
    for length in range(min(rgroups[library][0].keys()),max(rgroups[library][0].keys())+1):
      outfile.write("%d\t%d\t%d\n"%(length,rgroups[library][0][length],rgroups[library][1][length]))
  else:
    outfile.write('Length\tCounts\n')
    for length in range(min(rgroups[library][0].keys()),max(rgroups[library][0].keys())+1):
      outfile.write("%d\t%d\n"%(length,rgroups[library][0][length]))
  outfile.close()
