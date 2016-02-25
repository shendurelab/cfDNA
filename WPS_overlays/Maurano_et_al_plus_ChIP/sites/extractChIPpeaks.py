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
#ETS -> ETS -> ETS1,ELF1,ELK1,ELK4,STAT1,STAT2,STAT3,GABPA,SPI1
#IRF* -> IRF IRF/2 IRF/3 -> IRF1,STAT1,PRDM1 IRF3 IRF3
#MAFK -> MAFK -> MAFK
#MEF2A -> MEF2A MEF2A/2 -> MEF2A
#MYC/MAX -> MYC/MAX -> MYC,MAX,USF1,USF2,SREBP1
#PAX5 -> PAX5/2 -> PAX5
#RUNX* -> RUNX2 RUNX/AML -> RUNX3
#STAF-2 -> STAF/2 -> ZNF143
#TCF-LEF -> TCF/LEF -> TCF7L2,TCF3
#YY1 -> YY1 -> YY1

factorAssignment = { 'CTCF':['CTCF_Core/2'], 'TFAP2A':['AP-2','AP-2/2'], 'TFAP2C':['AP-2','AP-2/2'], 'E2F1':['E2F/2'], 'E2F4':['E2F/2'],  'EBF1':['EBF1'], 'TCF3':['Ebox','Ebox/CACCTG','TCF/LEF'], 'TAL1':['Ebox'], 'ESR1':['ESR1'], 'ETS1':['ETS'], 'ELF1':['ETS'], 'ELK1':['ETS'], 'ELK4':['ETS'], 'STAT1':['ETS'], 'STAT2':['ETS'], 'STAT3':['ETS'], 'GABPA':['ETS'], 'SPI1':['ETS'], 'IRF1':['IRF'], 'STAT1':['IRF'], 'PRDM1':['IRF'], 'IRF3':['IRF/2','IRF/3'], 'MAFK':['MAFK'], 'MEF2A':['MEF2A','MEF2A/2'], 'MYC':['MYC/MAX'], 'MAX':['MYC/MAX'], 'USF1':['MYC/MAX'], 'USF2':['MYC/MAX'], 'SREBP1':['MYC/MAX'], 'PAX5':['PAX5/2'], 'RUNX3':['RUNX2','RUNX/AML'], 'ZNF143':['STAF/2'], 'TCF7L2':['TCF/LEF'], 'YY1':['YY1'] }

outfiles = {}
for key,value in factorAssignment.iteritems():
  for factor in value:
    if factor not in outfiles: outfiles[factor] = open("%s.peaks"%(factor.replace("/","-")),'w')

infile = gzip.open("data/ENCODE_TfbsClusteredV3/tfbs.v3.tsv.gz")
for line in infile:
  fields = line.split("\t")
  if fields[3] in factorAssignment:
    for elem in factorAssignment[fields[3]]:
      outfiles[elem].write(line)

for key,value in outfiles.iteritems():
  value.close()
