#!/usr/bin/env python

"""
:Author: Martin Kircher
:Contact: mkircher@uw.edu
:Date: *13.04.2015
"""

import sys, os
from optparse import OptionParser
import math
from sortedcontainers import SortedList
import numpy as np

parser = OptionParser()
parser.add_option("-s","--smoother", dest="smoother", help="Use smoother on the window protection score (def Off)",default=False,action="store_true")
parser.add_option("-d","--debug", dest="debug", help="Turn debug mode on (produce no calls, just report adjusted values)",default=False,action="store_true")
(options, args) = parser.parse_args()

minlength = 50
maxlength = 150
chromNames = map(str,range(1,23))
chromNames.append("X")
chromNames.append("Y")
chromNames = set(chromNames)

windowSize = 1000
halfWindow = windowSize//2
windowValuesSorted = SortedList(load=int(math.sqrt(windowSize)))
windowValues = []
posInWindow = 0

variCutoff = 5

##############################
# Primitive window smoother
##############################
##v1:
#fweights = [0.1,0.2,0.4,0.2,0.1]
#lweights = [0.15,0.3,0.55]
#rweights = [0.55,0.3,0.15]
##v2:
#fweights = [0.15,0.20,0.3,0.20,0.15]
#lweights = [0.3,0.3,0.4]
#rweights = [0.4,0.3,0.3]
#
#def smoother(values,ctype=None): #True: left, False: right
  #global fweights, lweights, rweights
  #res = 0.0
  #if ctype == None:
    #for ind,weight in enumerate(fweights):
      #res+=weight*values[ind]
  #elif ctype == True:
    #for ind,weight in enumerate(lweights):
      #res+=weight*values[ind]
  #else:
    #for ind,weight in enumerate(rweights):
      #res+=weight*values[ind]
  #return res

######################################
# Savitzky-Golay filter as smoother
######################################
SGF_window_size = 21 # window_size size must be a positive odd number
SGF_order = 2 # max order: window_size - 1
SGF_order_range = range(SGF_order+1)
SGF_half_window = (SGF_window_size -1) // 2
# precompute coefficients
SGF_b = np.mat([[k**i for i in SGF_order_range] for k in range(-SGF_half_window, SGF_half_window+1)])
SGF_m = (np.linalg.pinv(SGF_b).A[0])[::-1]

def smoother(values,ctype=None): #True: fill-in left, False: fill-in right
  values = np.array(values)
  global SGF_m, SGF_half_window, SGF_window_size
 
  # pad the signal at the extremes with
  # values taken from the signal itself 
  # if signal window is too small
  if ctype == True:
    firstvals = values[0] - np.abs( values[1:min(SGF_window_size-len(values)+1,SGF_half_window+1)][::-1] - values[0] )
    y = np.concatenate((firstvals, values))
  elif ctype == False: 
    lastvals = values[-1] + np.abs( values[-min(SGF_window_size-len(values),SGF_half_window)-1:-1][::-1] - values[-1] )
    y = np.concatenate((values, lastvals))
  else:
    y = values
  return np.convolve( SGF_m, y, mode='valid')[0]


def quantiles(values,points):
  helper = list(values)
  helper.sort()
  res = []
  lVal = len(helper)
  for point in points:
    if round(lVal*point) == int(lVal*point):
      res.append(helper[int(lVal*point)])
    else:
      res.append((helper[int(round(lVal*point))]+helper[int(lVal*point)])*0.5)
  return res

def median(values):
  return quantiles(values,[0.5])[0]

def medianS(sortedValues):
  if len(sortedValues) % 2 == 0:
    return 0.5*(sortedValues[len(sortedValues)//2]+sortedValues[len(sortedValues)//2-1])
  else:
    return sortedValues[len(sortedValues)//2]

def quantileS(sortedValues,points):
  res = []
  lVal = len(sortedValues)
  for point in points:
    if round(lVal*point) == int(lVal*point):
      res.append(sortedValues[int(lVal*point)])
    else:
      res.append((sortedValues[int(round(lVal*point))]+sortedValues[int(lVal*point)])*0.5)
  return res

def continousWindows(region):
  res = []
  cstart,cend = None,None
  cmax = None
  #cmax,cmaxPos = None,None
  csum = 0
  for (pos,val) in region:
    if cmax == None:
      cend = pos
      cstart = pos
      cmax = val
      csum = val
      #cmaxPos = [pos]
    else:
      if cend+1 == pos:
        cend = pos
        csum += val
        if cmax < val:
          cmax = val
          #cmaxPos = [pos]
        #if cmax == val:
          #cmaxPos.append(pos)
      else:
        res.append((csum,cstart,cend,cmax)) # int(round(sum(cmaxPos)/float(len(cmaxPos))))
        cstart,cend = None,None
        #cmax,cmaxPos = None,None
        cmax = None
        csum = 0
  if cmax != None:
    res.append((csum,cstart,cend,cmax)) # int(round(sum(cmaxPos)/float(len(cmaxPos)))),
  return res

def evaluateValues(chrom,start,end,values,report = True):
  global minlength, maxlength, variCutoff
  global chromNames
  #if start != None:
    #sys.stdout.write("%s\t%d\t%d\t%d\n"%(chrom,start,end,len(values)))
  if (maxlength >= len(values) >= minlength):
    cMed = median(values)
    res = filter(lambda (x,y): y >= cMed, zip(range(start,end+1),values))
    res = continousWindows(res)
    res.sort()
    score,cstart,cend,cval = res[-1]
    if (chrom in chromNames) and report and (cval > variCutoff):
      cmiddle = cstart+round((cend-cstart)*0.5)
      sys.stdout.write("chr%s\t%d\t%d\t%s:%d-%d\t%d\t+\t%d\t%d\n"%(chrom,cstart-1,cend,chrom,cstart,cend,cval,cmiddle-1,cmiddle))
  
  elif (3*maxlength >= len(values) >= maxlength):
    cMed = median(values)
    res = filter(lambda (x,y): y >= cMed, zip(range(start,end+1),values))
    for (score,cstart,cend,cval) in continousWindows(res):
      if (maxlength >= (cend-cstart+1) >= minlength):
        if (chrom in chromNames) and report and (cval > variCutoff):
          cmiddle = cstart+round((cend-cstart)*0.5)
          sys.stdout.write("chr%s\t%d\t%d\t%s:%d-%d\t%d\t+\t%d\t%d\n"%(chrom,cstart-1,cend,chrom,cstart,cend,cval,cmiddle-1,cmiddle))
        
chrom = "1"
pos = 0
step = 1
cstart,cend = None,None
clist = []
cPass = not options.debug
ivalue = 0
for line in sys.stdin:
  # A new WIG region starts
  if line.startswith("fixedStep"):
    #Decide wheather we need to clean-up or whether this block continues 
    fields = line.split()
    nchrom = fields[1].split("=")[1].replace("chr","")
    npos = int(fields[2].split("=")[1])-1
    nstep = int(fields[3].split("=")[1])
    #if True:
    if (nchrom == chrom) and (nstep == step) and (npos == pos): continue
    else:
      # Make sure that we processed all positions in the last WIG region
      if (len(windowValues) == windowSize):
        posInWindow += 1
        #windowLow,windowMedian,windowHigh = quantileS(windowValuesSorted,[0.25,0.5,0.75])
        windowMedian = medianS(windowValuesSorted)
        while (posInWindow < windowSize):
          # Evaluate adjusted value
          if options.smoother:
            helper = windowValues[max(0,posInWindow-SGF_half_window):posInWindow+SGF_half_window+1]
            if len(helper) == SGF_window_size:
              ivalue = smoother(helper)-windowMedian
            else:
              helper = windowValues[max(0,posInWindow-SGF_half_window):posInWindow+1]
              ivalue = smoother(helper,ctype=False)-windowMedian
          else:
            ivalue = windowValues[posInWindow]-windowMedian
          ipos = pos-windowSize+posInWindow+1
          if options.debug:
            sys.stdout.write("%s\t%d\t%.1f\t%.1f\t(%d)\n"%(chrom,ipos,ivalue,windowValues[posInWindow],posInWindow))
          #print chrom,ipos,ivalue
          if ivalue > 0:
            if (cend != None) and (ipos <= (cend+5*step)):
              while (cend+step < ipos):
                cend+=step
                clist.append(0)
              clist.append(ivalue)
              cend = ipos
              #cPass = cPass and (windowHigh-windowLow >= variCutoff)
            else:
              evaluateValues(chrom,cstart,cend,clist,cPass)
              clist = [ivalue]
              cstart = ipos
              cend = ipos
              cPass = not options.debug #and (windowHigh-windowLow >= variCutoff)
          # Increase position in window towards end of the list
          posInWindow += 1
        evaluateValues(chrom,cstart,cend,clist,cPass)
      # Clean-up for next region
      chrom = nchrom
      pos = npos
      step = nstep
      cstart,cend = None,None
      cPass = not options.debug
      clist = []
      windowValues = []
      windowValuesSorted.clear()
      posInWindow = 0
  else:
    if len(line.strip()) == 0: continue 
    pos += step
    value = int(line)
    
    # Window buffer has not reached its size yet
    if len(windowValues) < windowSize:
      windowValues.append(value)
      windowValuesSorted.add(value)
    # Window buffer has reached its size
    elif len(windowValues) == windowSize:
      # Buffer filled but current position is not yet centered in the window,
      # work with current buffer state
      if (posInWindow < halfWindow):
        #windowLow,windowMedian,windowHigh = quantileS(windowValuesSorted,[0.25,0.5,0.75])
        windowMedian = medianS(windowValuesSorted)
        while (posInWindow <= halfWindow):
          # Evaluate adjusted value
          if options.smoother:
            helper = windowValues[max(0,posInWindow-SGF_half_window):posInWindow+SGF_half_window+1]
            if len(helper) == SGF_window_size:
              ivalue = smoother(helper)-windowMedian
            else:
              helper = windowValues[posInWindow:posInWindow+SGF_half_window+1]
              ivalue = smoother(helper,ctype=True)-windowMedian
          else:
            ivalue = windowValues[posInWindow]-windowMedian
          ipos = pos-windowSize+posInWindow
          if options.debug:
            sys.stdout.write("%s\t%d\t%.1f\t%.1f\t(%d)\n"%(chrom,ipos,ivalue,windowValues[posInWindow],posInWindow))
          #print chrom,ipos,ivalue
          if ivalue > 0:
            if (cend != None) and (ipos <= (cend+5*step)):
              while (cend+step < ipos):
                cend+=step
                clist.append(0)
              clist.append(ivalue)
              cend=ipos
              #cPass = cPass and (windowHigh-windowLow >= variCutoff)
            else:
              evaluateValues(chrom,cstart,cend,clist,cPass)
              clist = [ivalue]
              cstart = ipos
              cend = ipos
              cPass = not options.debug #and (windowHigh-windowLow >= variCutoff)
          # Increase position in window towards center
          posInWindow += 1
        # Overstepped by 1, reset to center
        posInWindow -= 1
        
      # We are centered
      windowValuesSorted.discard(windowValues.pop(0))
      windowValues.append(value)
      windowValuesSorted.add(value)
      #windowLow,windowMedian,windowHigh = quantileS(windowValuesSorted,[0.25,0.5,0.75])
      windowMedian = medianS(windowValuesSorted)
      # Evaluate adjusted value
      if options.smoother:
        helper = windowValues[max(0,posInWindow-SGF_half_window):posInWindow+SGF_half_window+1]
        ivalue = smoother(helper)-windowMedian
      else:
        ivalue = windowValues[posInWindow]-windowMedian
      ipos = pos-windowSize+posInWindow+1
      if options.debug:
        sys.stdout.write("%s\t%d\t%.1f\t%.1f\t(%d)\n"%(chrom,ipos,ivalue,windowValues[posInWindow],posInWindow))
      #print chrom,ipos,ivalue
      if ivalue > 0:
        if (cend != None) and (ipos <= (cend+5*step)):
          while (cend+step < ipos):
            cend+=step
            clist.append(0)
          clist.append(ivalue)
          cend = ipos
          #cPass = cPass and (windowHigh-windowLow >= variCutoff)
        else:
          evaluateValues(chrom,cstart,cend,clist,cPass)
          clist = [ivalue]
          cstart = ipos
          cend = ipos
          cPass = not options.debug #and (windowHigh-windowLow >= variCutoff)

# Make sure that we processed all positions in the last WIG region
if (len(windowValues) == windowSize):
  posInWindow += 1
  #windowLow,windowMedian,windowHigh = quantileS(windowValuesSorted,[0.25,0.5,0.75])
  windowMedian = medianS(windowValuesSorted)
  while (posInWindow < windowSize):
    # Evaluate adjusted value
    if options.smoother:
      helper = windowValues[max(0,posInWindow-SGF_half_window):posInWindow+SGF_half_window+1]
      if len(helper) == SGF_window_size:
        ivalue = smoother(helper)-windowMedian
      else:
        helper = windowValues[max(0,posInWindow-SGF_half_window):posInWindow+1]
        ivalue = smoother(helper,ctype=False)-windowMedian
    else:
      ivalue = windowValues[posInWindow]-windowMedian
    ipos = pos-windowSize+posInWindow+1
    if options.debug:
      sys.stdout.write("%s\t%d\t%d\t%.1f\t(%d)\n"%(chrom,ipos,ivalue,windowValues[posInWindow],posInWindow))
    #print chrom,ipos,ivalue
    if ivalue > 0:
      if (cend != None) and (ipos <= (cend+5*step)):
        while (cend+step < ipos):
          cend+=step
          clist.append(0)
        clist.append(ivalue)
        cend = ipos
        #cPass = cPass and (windowHigh-windowLow >= variCutoff)
      else:
        evaluateValues(chrom,cstart,cend,clist,cPass)
        clist = [ivalue]
        cstart = ipos
        cend = ipos
        cPass = not options.debug #and (windowHigh-windowLow >= variCutoff)
    # Increase position in window towards end of the list
    posInWindow += 1
  evaluateValues(chrom,cstart,cend,clist,cPass)
