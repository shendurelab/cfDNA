#!/usr/bin/env python

"""
:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *30.04.2015
"""

import sys, os
from optparse import OptionParser
from sortedcontainers import SortedList
import math

parser = OptionParser()
parser.add_option("-w","--WPS", dest="WPS", help="Window protection score file (def WPS.wig)",default="WPS.wig")
parser.add_option("-c","--Coverage", dest="Coverage", help="Coverage score file (def Coverage.wig)",default="Coverage.wig")
parser.add_option("-v","--verbose", dest="verbose", help="Turn on verbose messages on stderr",default=False,action="store_true")
(options, args) = parser.parse_args()

windowSize = 1000
halfWindow = windowSize//2
windowValuesSorted = SortedList(load=int(math.sqrt(windowSize)))
windowValues = []
posInWindow = 0

def medianS(sortedValues):
  if len(sortedValues) % 2 == 0:
    return 0.5*(sortedValues[len(sortedValues)//2]+sortedValues[len(sortedValues)//2-1])
  else:
    return sortedValues[len(sortedValues)//2]

chrom = "1"
pos = 0
step = 1
cstart,cend = None,None
ivalue = 0

if not os.path.exists(options.Coverage) or not os.path.exists(options.WPS):
  sys.stderr.write("ERROR: Could not access input files\n")
  sys.exit()

#count = 0

sys.stdout.write("#Chrom\tPos\tnmWPS\tmWPS\tWPS\tCov\n")

cov = open(options.Coverage)
wps = open(options.WPS)
cline = cov.readline()
wpsline = wps.readline()
while cline != '' or wpsline != '':
  # Skip empty lines
  while (len(wpsline.strip()) == 0) and (wpsline != ''):
    wpsline = wps.readline()
  while (len(cline.strip()) == 0) and (cline != ''):
    cline = cov.readline()

  #count += 1
  #if count % 100000 == 0:
    #sys.stderr.write("Read %d lines; %d %d\n"%(count,len(windowValuesSorted),len(windowValues)))

  # A new WIG region starts
  if wpsline.startswith("fixedStep") and cline.startswith("fixedStep"):
    #Decide wheather we need to clean-up or whether this block continues 
    fields = wpsline.split()
    nchrom = fields[1].split("=")[1].replace("chr","")
    npos = int(fields[2].split("=")[1])-1
    nstep = int(fields[3].split("=")[1])
    if (nchrom == chrom) and (nstep == step) and (npos == pos): 
      if options.verbose: sys.stderr.write("Continued: %s %s\n"%(wpsline.rstrip(),cline.rstrip()))
    else:
      if options.verbose: sys.stderr.write("Both reset: %s %d %d vs %s %d %d\n"%(nchrom,npos,nstep,chrom,pos,step))
      # Make sure that we processed all positions in the last WIG region
      if (len(windowValues) == windowSize):
        posInWindow += 1
        windowMedian = medianS(windowValuesSorted)
        while (posInWindow < windowSize):
          WPSval,Cov = windowValues[posInWindow]
          #ivalue = min(100.0,max(-100.0,(0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100))
          ivalue = (0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100
          ipos = pos-windowSize+posInWindow+1
          sys.stdout.write("%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n"%(chrom,ipos,ivalue,WPSval-windowMedian,WPSval,Cov))
          # Increase position in window towards end of the list
          posInWindow += 1
      # Clean-up for next region
      chrom = nchrom
      pos = npos
      step = nstep
      cstart,cend = None,None
      windowValues = []
      windowValuesSorted.clear()
      posInWindow = 0
    cline = cov.readline()
    wpsline = wps.readline()
  elif wpsline.startswith("fixedStep") and not cline.startswith("fixedStep"):
    #Decide wheather we need to clean-up or whether this block continues 
    fields = wpsline.split()
    nchrom = fields[1].split("=")[1].replace("chr","")
    npos = int(fields[2].split("=")[1])-1
    nstep = int(fields[3].split("=")[1])
    if (nchrom == chrom) and (nstep == step) and (npos == pos): 
      if options.verbose: sys.stderr.write("Continued [WPS]: %s %s\n"%(wpsline.rstrip(),cline.rstrip()))
    else:
      if options.verbose: sys.stderr.write("WPS reset: %s %d %d vs %s %d %d\n"%(nchrom,npos,nstep,chrom,pos,step))
      # Make sure that we processed all positions in the last WIG region
      if (len(windowValues) == windowSize):
        posInWindow += 1
        windowMedian = medianS(windowValuesSorted)
        while (posInWindow < windowSize):
          WPSval,Cov = windowValues[posInWindow]
          #ivalue = min(100.0,max(-100.0,(0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100))
          ivalue = (0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100
          ipos = pos-windowSize+posInWindow+1
          sys.stdout.write("%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n"%(chrom,ipos,ivalue,WPSval-windowMedian,WPSval,Cov))
          # Increase position in window towards end of the list
          posInWindow += 1
      # Clean-up for next region
      chrom = nchrom
      pos = npos
      step = nstep
      cstart,cend = None,None
      windowValues = []
      windowValuesSorted.clear()
      posInWindow = 0
    wpsline = wps.readline()
  elif not wpsline.startswith("fixedStep") and cline.startswith("fixedStep"):
    #Decide wheather we need to clean-up or whether this block continues 
    fields = cline.split()
    nchrom = fields[1].split("=")[1].replace("chr","")
    npos = int(fields[2].split("=")[1])-1
    nstep = int(fields[3].split("=")[1])
    if (nchrom == chrom) and (nstep == step) and (npos == pos): 
      if options.verbose: sys.stderr.write("Continued [Cov]: %s %s\n"%(wpsline.rstrip(),cline.rstrip()))
    else:
      if options.verbose: sys.stderr.write("Coverage reset: %s %d %d vs %s %d %d\n"%(nchrom,npos,nstep,chrom,pos,step))
      # Make sure that we processed all positions in the last WIG region
      if (len(windowValues) == windowSize):
        posInWindow += 1
        windowMedian = medianS(windowValuesSorted)
        while (posInWindow < windowSize):
          WPSval,Cov = windowValues[posInWindow]
          #ivalue = min(100.0,max(-100.0,(0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100))
          ivalue = (0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100
          ipos = pos-windowSize+posInWindow+1
          sys.stdout.write("%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n"%(chrom,ipos,ivalue,WPSval-windowMedian,WPSval,Cov))
          # Increase position in window towards end of the list
          posInWindow += 1
      # Clean-up for next region
      chrom = nchrom
      pos = npos
      step = nstep
      cstart,cend = None,None
      windowValues = []
      windowValuesSorted.clear()
      posInWindow = 0
    cline = cov.readline()
  else:
    pos += step
    value = int(wpsline),int(cline)
    
    # Window buffer has not reached its size yet
    if len(windowValues) < windowSize:
      windowValues.append(value)
      windowValuesSorted.add(value[0])
    # Window buffer has reached its size
    elif len(windowValues) == windowSize:
      # Buffer filled but current position is not yet centered in the window,
      # work with current buffer state
      if (posInWindow < halfWindow):
        windowMedian = medianS(windowValuesSorted)
        while (posInWindow <= halfWindow):
          # Evaluate adjusted value
          WPSval,Cov = windowValues[posInWindow]
          #ivalue = min(100.0,max(-100.0,(0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100))
          ivalue = (0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100
          ipos = pos-windowSize+posInWindow
          sys.stdout.write("%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n"%(chrom,ipos,ivalue,WPSval-windowMedian,WPSval,Cov))
          posInWindow += 1
        # Overstepped by 1, reset to center
        posInWindow -= 1
        
      # We are centered
      windowValuesSorted.discard(windowValues.pop(0)[0])
      windowValues.append(value)
      windowValuesSorted.add(value[0])
      windowMedian = medianS(windowValuesSorted)
      # Evaluate adjusted value
      WPSval,Cov = windowValues[posInWindow]
      #ivalue = min(100.0,max(-100.0,(0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100))
      ivalue = (0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100
      ipos = pos-windowSize+posInWindow+1
      sys.stdout.write("%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n"%(chrom,ipos,ivalue,WPSval-windowMedian,WPSval,Cov))

    cline = cov.readline()
    wpsline = wps.readline()

# Make sure that we processed all positions in the last WIG region
if (len(windowValues) == windowSize):
  posInWindow += 1
  windowMedian = medianS(windowValuesSorted)
  while (posInWindow < windowSize):
    # Evaluate adjusted value
    WPSval,Cov = windowValues[posInWindow]
    #ivalue = min(100.0,max(-100.0,(0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100))
    ivalue = (0.0 if Cov == 0 else (WPSval-windowMedian)/float(Cov))*100
    ipos = pos-windowSize+posInWindow+1
    sys.stdout.write("%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\n"%(chrom,ipos,ivalue,WPSval-windowMedian,WPSval,Cov))
    posInWindow += 1

