#!/usr/bin/env python

"""

:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *10.03.2015
"""

import sys, os
from optparse import OptionParser

parser = OptionParser("%prog [options]")
parser.add_option("--PWM_mapping", dest="PWM_mapping", help="PWM mapping file (def pwm.gene.mapping.V3.withclusters.txt)",default="pwm.gene.mapping.V3.withclusters.txt")
#parser.add_option("--PWM_cluster", dest="PWM_cluster", help="Quality score character to be set for adapter bases (def TableS11_PWM.bycluster.txt)",default="TableS11_PWM.bycluster.txt")
(options, args) = parser.parse_args()

## pwm.gene.mapping.V3.withclusters.txt
#Motif_Database  PWM_Name        Gene_Name       Cluster       Strand  Offset  Cluster_number
#CTCF_upstream   CTCF_Core       CTCF    c251-CTCF_Core/2        +       1       251
#CTCF_upstream   CTCF_Upstream   CTCF    c138-CTCF_Upstream      +       0       138
#CTCF_upstream   CTCF_UpstreamP1 CTCF    c12-CTCF_UpstreamP1     +       0       12
#JASPAR          MA0002.1-RUNX1  RUNX1   c188-RUNX/AML           +       0       188
#JASPAR          MA0002.2-RUNX1  RUNX1   c188-RUNX/AML           +       0       188
#JASPAR          MA0003.1-TFAP2A TFAP2A  c203-AP-2/2             +       0       203
#JASPAR          MA0004.1-ARNT   ARNT    c90-MYC/MAX             +       -2      90

## TableS11_PWM.bycluster.txt
#Cluster Cluster.num     Motifs
#1.P53        1       V_P53_05
#2.RORA1      2       V_RORA1_01
#3.PR         3       V_PR_02
#4.PAX5       4       V_PAX5_01
#5.PAX4       5       MA0068.1-Pax4
#7.HNF1       7       V_HNF1_Q6_01
#8.PAX4/2     8       V_PAX4_04

## INPUT
#chr1    10096   10114   ESRRG_nuclearreceptor_1         6.25347e-05     -       TAACCCTAACCCAACCC|T
#chr1    10106   10120   Rxra.mouse_nuclearreceptor_1    9.41875e-05     -       CCAACCCTAACCCT
#chr1    10106   10120   Rxrb.mouse_nuclearreceptor_1    5.84544e-05     -       CCAACCCTAACCCT
#chr1    10131   10151   TBX1_TBX_1                      9.39743e-05     -       TAACCCTAACCCTAACCCCT
#chr1    10142   10157   NR2F1_nuclearreceptor_1         1.76882e-05     -       CTAACCCCTAACCCT
#chr1    10143   10157   NR2C2_nuclearreceptor_1         1.25498e-05     -       TAACCCCTAACCCT
#chr1    10143   10157   NR2F6_nuclearreceptor_2         1.83942e-05     -       TAACCCCTAACCCT
#chr1    10143   10157   NR2F6_nuclearreceptor_3         3.45994e-05     -       TAACCCCTAACCCT
#chr1    10143   10157   Nr2f6.mouse_nuclearreceptor_2   3.27350e-05     -       TAACCCCTAACCCT
#chr1    10142   10159   ZNF282_C2H2_1                   1.24456e-05     +       CTAACCCCTAACCCTAA

## INPUT PWM
#Motif_Database  PWM_Name        Gene_Name       Cluster Strand  Offset  Cluster_number
#taipale_selex   ESRRG_nuclearreceptor_1         ESRRG   c64.ESRRA/2     +       1       64
#taipale_selex   Rxra.mouse_nuclearreceptor_1    Rxra    c85-NR/2        +       0       85
#taipale_selex   Rxrb.mouse_nuclearreceptor_1    Rxrb    c85-NR/2        +       0       85
#taipale_selex   TBX1_TBX_1                      TBX1    c91-T-box       +       -4      91
#taipale_selex   NR2F1_nuclearreceptor_1         NR2F1   c85-NR/2        +       0       85
#taipale_selex   NR2C2_nuclearreceptor_1         NR2C2   c85-NR/2        +       0       85
#taipale_selex   NR2F6_nuclearreceptor_2         NR2F6   c85-NR/2        +       0       85
#taipale_selex   NR2F6_nuclearreceptor_3         NR2F6   c85-NR/2        +       0       85
#taipale_selex   Nr2f6.mouse_nuclearreceptor_2   Nr2f6   c85-NR/2        +       0       85
#taipale_selex   ZNF282_C2H2_1                   ZNF282  c165.ZNF282     -       0       165

## OUTPUT
#1       10096   10114   ESRRA/2 -       17
#1       10106   10120   NR/2    -       14
#1       10131   10151   T-box   -       24
#1       10142   10157   NR/2    -       15
#1       10142   10159   ZNF282  +       17

## SECOND CASE

## INPUT
#chr1    857025  857048  TBX1_TBX_4      3.71544e-05     -       GCACGGACAGTACAGGTGTGAGC
#chr1    857038  857046  TBX1_TBX_3      1.31552e-05     +                    AGGTGTGA

## INPUT PWM
#taipale_selex  TBX1_TBX_3     TBX1           c91-T-box      +              -4             91
#taipale_selex  TBX1_TBX_4     TBX1           c91-T-box      -              9              91

#chr1    857025  857048  TBX1_TBX_4      3.71544e-05     -       GCACGGACAGTACAGGTGTGAGC
                                                                         #*
                                                               #'01234567890123456789012'
                                                                         #*
#chr1    857038  857046  TBX1_TBX_3      1.31552e-05     +                    AGGTGTGA

## OUTPUT
#1       857025  857048  T-box   -       9


#def unifyLineBuffers(lines):
  #global curChrom,curPos
  #lines.sort()
  #cCluster,cChrom,cStart,cEnd,cStrand,cOffset = None,"",0,0,"",0
  #results = []
  #for (PWMcluster,chrom,start,end,strand,PWMoffset) in lines:
    #if (PWMcluster == cCluster) and (chrom == cChrom) and (min(end,cEnd)-max(start,cStart)>0):
      #if start < cStart:
        #cStart = start
        #cOffset = PWMoffset
      #if end > cEnd:
        #cEnd = end
    #else:
      #if cCluster != None:
        #results.append((cChrom,cStart,cEnd,cCluster,cStrand,cOffset))
      #cCluster,cChrom,cStart,cEnd,cStrand,cOffset = PWMcluster,chrom,start,end,strand,PWMoffset
  #if cCluster != None:
    #results.append((cChrom,cStart,cEnd,cCluster,cStrand,cOffset))
  #results.sort()
  #lines = []
  #for (cChrom,cStart,cEnd,cCluster,cStrand,cOffset) in results:
    #if (cChrom != curChrom) or ((cChrom == curChrom) and (curPos > cEnd+50)):
      #print "%s\t%d\t%d\t%s\t%d\t%s"%(cChrom,cStart,cEnd,cCluster,cOffset,cStrand)
    #else:
      #lines.append((cCluster,cChrom,cStart,cEnd,cStrand,cOffset))
  #return lines


def unifyLineBuffers(lines):
  global curChrom,curPos
  lines.sort()
  cCluster,cChrom,cStart,cEnd,cStrand = None,"",0,0,""
  results = []
  for (PWMcluster,chrom,start,end,strand) in lines:
    if (PWMcluster == cCluster) and (chrom == cChrom) and (min(end,cEnd)-max(start,cStart)>0):
      if start < cStart:
        cStart = start
      if end > cEnd:
        cEnd = end
    else:
      if cCluster != None:
        results.append((cChrom,cStart,cEnd,cCluster,cStrand))
      cCluster,cChrom,cStart,cEnd,cStrand = PWMcluster,chrom,start,end,strand
  if cCluster != None:
    results.append((cChrom,cStart,cEnd,cCluster,cStrand))
  results.sort()
  lines = []
  for (cChrom,cStart,cEnd,cCluster,cStrand) in results:
    if (cChrom != curChrom) or ((cChrom == curChrom) and (curPos > cEnd+50)):
      print "%s\t%d\t%d\t%s\t%s"%(cChrom,cStart,cEnd,cCluster,cStrand)
    else:
      lines.append((cCluster,cChrom,cStart,cEnd,cStrand))
  return lines

PWMinfos = {}
motifClusters = {}

infile = open(options.PWM_mapping)
infile.readline()
for line in infile:
  motifDB,PWM,gene,cluster,strand,offset,clusterNr = line.rstrip().split("\t")
  gene=gene.upper()
  PWMinfos[PWM] = strand,int(offset)
  #motifClusters[PWM] = cluster[(len(clusterNr)+2):]
  motifClusters[PWM] = gene
infile.close()

#print len(PWMinfos),len(motifClusters),len(set(motifClusters.values()))

lineBuffer = []

curChrom,curPos = "",0
for line in sys.stdin:
  chrom,start,end,PWM,pVal,strand,seq = line.rstrip().split("\t")
  chrom = chrom.replace("chr","")
  start,end = int(start),int(end)
  curChrom,curPos=chrom,start

  if PWM in motifClusters: # PWM in PWMinfos and 
    PWMcluster = motifClusters[PWM]
    #PWMstrand,PWMoffset = PWMinfos[PWM]
    #if PWMstrand != strand:
      #PWMoffset = (end-start)-PWMoffset
    #lineBuffer.append((PWMcluster,chrom,start,end,strand,PWMoffset))
    lineBuffer.append((PWMcluster,chrom,start,end,strand))
  else:
    sys.stderr.write("Error: motif not found\n")

  if len(lineBuffer) > 100:
    lineBuffer = unifyLineBuffers(lineBuffer)

lineBuffer = unifyLineBuffers(lineBuffer)      
#for (PWMcluster,chrom,start,end,strand,PWMoffset) in lineBuffer:
  #print "%s\t%d\t%d\t%s\t%d\t%s"%(chrom,start,end,PWMcluster,PWMoffset,strand)
for (PWMcluster,chrom,start,end,strand) in lineBuffer:
  print "%s\t%d\t%d\t%s\t%s"%(chrom,start,end,PWMcluster,strand)
