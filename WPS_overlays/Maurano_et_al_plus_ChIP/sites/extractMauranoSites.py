#!/usr/bin/env python

"""
:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *27.06.2015
"""

import sys, os
from optparse import OptionParser
import gzip

parser = OptionParser()
(options, args) = parser.parse_args()

#CTCF -> CTCF_Core/2 -> CTCF
#AP-2s -> AP-2 AP-2/2 -> TFAP2A,TFAP2C
#E2F-2 -> E2F/2 -> E2F1,E2F4
#EBF1 -> EBF1 -> EBF1
#EBOX -> Ebox Ebox/CACCTG -> TCF3,TAL1 TCF3
#ESR1 -> ESR1 -> ESR1
#ETS -> ETS -> ETS1 ELF1 ELK1 ELK4 STAT1 STAT2 STAT3 GABPA SPI1
#IRF* -> IRF IRF/2 IRF/3 -> IRF1,STAT1,PRDM1 IRF3 IRF3
#MAFK -> MAFK -> MAFK
#MEF2A -> MEF2A MEF2A/2 -> MEF2A
#MYC/MAX -> MYC/MAX -> MYC,MAX,USF1,USF2,SREBP1
#PAX5 -> PAX5/2 -> PAX5
#RUNX* -> RUNX2 RUNX/AML -> RUNX3
#STAF-2 -> STAF/2 -> ZNF143
#TCF-LEF -> TCF/LEF -> TCF7L2,TCF3
#YY1 -> YY1 -> YY1

factors = set(['CTCF_Core/2','AP-2','AP-2/2','E2F/2','EBF1','Ebox','Ebox/CACCTG','ESR1','ETS','IRF','IRF/2','IRF/3','MAFK','MEF2A','MEF2A/2','MYC/MAX','PAX5/2','RUNX2','RUNX/AML','STAF/2','TCF/LEF','YY1'])

outfiles = {}
for factor in factors:
  outfiles[factor] = open("%s.sites"%(factor.replace("/","-")),'w')

infile = gzip.open("data/Maurano_et_al_func_var/hg19.taipale.recluster.bed.gz")
for line in infile:
  fields = line.split("\t")
  if fields[3] in factors:
    outfiles[fields[3]].write(line)

for factor in factors:
  outfiles[factor].close()
