#!/usr/bin/env python
# -*- coding: ASCII -*-

"""

Merge/Adapter trim reads stored in BAM

:Author: Martin Kircher
:Contact: mkircher@uw.edu

"""

import sys
import os
import math
import pysam
from optparse import OptionParser,OptionGroup
import string
table = string.maketrans('TGCAMRWSYKVHDBtgcamrwsykvhdb','ACGTKYWSRMBDHVacgtkywsrmbdhv') # COMPLEMENT DNA

quality_offset =  33

parser = OptionParser("%prog [options] BAMfile")
parser.add_option("-p","--PIPE",dest="pipe",help="Read BAM from and write it to PIPE",default=False,action="store_true")
parser.add_option("-o", "--outdir", dest="outdir", help="Create output files in another directory.")
parser.add_option("-c", "--consensus", dest="consensus", help="Report PCR duplicate consensus instead of sequence with highest sum of base qualities.",default=False,action="store_true")
parser.add_option("", "--outprefix", dest="outprefix", help="Prefix for output files (default 'PCRconsensus').",default="PCRconsensus")
parser.add_option("", "--outnewconsensus", dest="outnewconsensus", help="Save reads with new consensus sequence for realignment to separate FastQ output files with this prefix (default OFF).",default=None)
parser.add_option("", "--SAM", dest="SAM", help="Input/Output SAM not BAM.",default=False,action="store_true")
parser.add_option("", "--ignore_RG", dest="ignore_RG", help="Ignore the RG when looking for PCR duplicates. The consensus reads gets the RG of the template used.",default=False,action="store_true")
parser.add_option("", "--fixID", dest="fixID", help="Fix read ID, take only first part up to first / character",default=False,action="store_true")
parser.add_option("--library", dest="library", help="Use library name from RG read header rather than the readgroup ID",default=False,action="store_true")
parser.add_option("--UMI", dest="UMI", help="Use unique molecule identifier (UMI, in second index field) for grouping ",default=False,action="store_true")
parser.add_option("-v", "--verbose", dest="verbose", help="Turn all messages on",default=False,action="store_true")

group = OptionGroup(parser, "Filter options")
group.add_option("--include_qcfail",dest="include_qcfail",help="Consider reads that have the QC fail flag",default=False,action='store_true')
group.add_option("--rescue_qcfail",dest="rescue_qcfail",help="Remove fail quality flag of reads if PCR duplicates are observed",default=False,action='store_true')
group.add_option("--frequency_cutoff",dest="frequency",help="Keep only sequences with at least X PCR duplicates (default X = None)",default=None,type="int")
group.add_option("-f","--5p",dest="fivePrime",help="Cluster reads on five prime coordinate",default=False,action='store_true')
group.add_option("--max_length",dest="MaxLength",help="Longest possible read length stored in the input BAM (only relevant for SR 5' clustering and PE reads spanning multiple contigs, def 800)",default=800,type="int")
group.add_option("-k","--keep",dest="keep",help="Keep unmapped sequences",default=False,action='store_true')
group.add_option("--buffer",dest="rbuffer",help="Lowest number of PE reads buffered before write (def 5000)",default=5000,type="int")
group.add_option("-m","--merged",dest="merged",help="Keep only SR reads that are merged",default=False,action='store_true')
parser.add_option_group(group)
(options, args) = parser.parse_args()


#------------------------------
# PE data
#------------------------------

#The --5p parameter does not trigger any special code in the PE handling. The outer coordinates of PE reads are defined by chromosome, coordinate forward read and coordinate reverse read. 

#(1) chromosome is a string or can be a tuple of strings for PEs mapped across contigs/chromosomes
#(2) Coordinate of a read aligned without reversal: reported five prime position in the BAM
#(3) Coordinate of a read aligned as reverse complement: five prime position + alignment length

#When PE reads are collected, they are first collected incomplete pairs and then as complete pairs. Complete pairs are stored until a buffer limit is reached, incomplete pairs until the end of the script -- where they are essentially forgotten, ehh, removed from the BAM file ;-) . So if the buffer is full, a PE read cluster is processed if:

#(1) The chromosome is a string and we are already on a different chromosome
#(2) or, we are on the same chromosome, but more than molecule length away from the largest 5' coordinate of the cluster
#(3) or, chromosome is a tuple and none of the strings in the tuple matches the current chromosome
#(4) or, one of the strings is the current chromosome (implicating that we have finished the other chromosome since we have complete PE read clusters at hand) and we are more than 5p_max_length bases away from the 5' position of that read.

#------------------------------
#Merged parameter
#------------------------------

#Removes reads that are not flagged "paired in sequencing" and where the read ID does not start in "M_". If your incomplete PE reads do not have the paired in sequencing flag, this would remove them. If they are marked as paired in sequencing, but the second read is just missing from the file that would of course not help -- but those are by default removed (see above).


def extractLengthCigar(cigarlist):
  #Op BAM Description
  #M 0 alignment match (can be a sequence match or mismatch)
  #I 1 insertion to the reference
  #D 2 deletion from the reference
  #N 3 skipped region from the reference
  #S 4 soft clipping (clipped sequences present in SEQ)
  #H 5 hard clipping (clipped sequences NOT present in SEQ)
  #P 6 padding (silent deletion from padded reference)
  #= 7 sequence match
  #X 8 sequence mismatch

  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 1 or operation == 4 or operation == 7 or operation == 8: tlength+=length
  return tlength

def aln_length(cigarlist):
  tlength = 0
  for operation,length in cigarlist:
    if operation == 0 or operation == 2 or operation == 3 or operation >= 6: tlength += length
  return tlength

def calc_consensus(reads):
  if options.consensus:
    # DETERMINE CONSENSUS SEQUENCE FOR THIS VARIANT
    if len(reads) > 1:
      seqstring,qualstring = "",""
      for pos in xrange(len(reads[0].seq)):
        bases = [0,0,0,0]
        count = 0
        base,qualchr = None,None
        for elem in reads:
          base = elem.seq[pos]
          qualchr = elem.qual[pos]
          if base == 'N': continue
          count += 1
          qual = (ord(elem.qual[pos])-quality_offset)/-10.0
          if qual == 0: qual = -0.1
          rev_qual = math.log10(1.0-10**qual)-math.log10(3.0)
          for i,b in enumerate('ACGT'):
            if b == base: bases[i]+=qual
            else: bases[i]+=rev_qual
        if count > 0:
          total_prob = math.log10(max(0.000000001,sum(map(lambda x:10**x,bases))))
          max_base,max_qual,min_val = 'N',chr(quality_offset),0
          for i,b in enumerate('ACGT'):
            cval = bases[i]-total_prob
            if cval < min_val:
              min_val = cval
              max_base,max_qual = b, chr(int(round(min(60,-10.0*(cval))))+quality_offset)
          seqstring+=max_base
          qualstring+=max_qual
        else:
          seqstring+=base
          qualstring+=qualchr
    else:
      seqstring = reads[0].seq
      qualstring = reads[0].qual

    count = len(reads)
    outread = None
    minmapq = 255
    for read in reads: # REUSE READ OBJECT FROM ENTRY WITH IDENTICAL SEQUENCE BUT ALWAYS REPORT ID OF FIRST READ
      ofields = None
      if read.mapq < minmapq: minmapq = read.mapq
      if read.tags != None:
        ofields = []
        for key,value in read.tags: # LOOK FOR PREVIOUS PCR DUPLICATE COUNTS
          if key == "XP": count += value
          else: ofields.append((key,value))
      #sys.stderr.write('%s %s %s\n'%(read.seq==seqstring,read.seq,seqstring))
      if read.seq == seqstring:
        outread = read
        outread.qname = reads[0].qname
        outread.qual = qualstring
        if ofields != None: outread.tags = ofields
    if outread == None: # CONSENSUS SEQUENCE DOES NOT MATCH ONE OF THE ORIGINAL READS: WE DO NOT KNOW THAT READ IS STILL ALIGNED CORRECTLY!
      #if options.verbose:
        #sys.stderr.write('Consensus sequence does not match one of the original %d reads.\n'%(len(reads)))
        #sys.stderr.write('%s   %s\n'%(seqstring,qualstring))
        #for read in reads:
          #sys.stderr.write('%s   %s\n'%(read.seq,read.qual))
      outread = reads[0]
      outread.is_unmapped = True
      outread.mapq = minmapq
      outread.seq = seqstring
      outread.qual = qualstring
  else: # REUSE READ OBJECT WITH HIGHEST SUM OF QUALITIES BUT ALWAYS REPORT ID OF FIRST READ
    count = len(reads)
    outread = None
    maxsumqual = 0
    for read in reads:
      nsum = sum(map(ord,read.qual))
      if nsum  > maxsumqual:
        outread = read
        maxsumqual = nsum
      if read.tags != None:
        ofields = []
        for key,value in read.tags: # LOOK FOR PREVIOUS PCR DUPLICATE COUNTS
          if key == "XP": count += value
          else: ofields.append((key,value))
        outread.tags = ofields
    outread.qname = reads[0].qname

  outread.is_duplicate = False
  if outread.tags == None: outread.tags = [("XP",count)]
  else: outread.tags = outread.tags+[("XP",count)]

  if options.frequency == None or count >= options.frequency:
    return outread
  else:
    return None


def get_consensus(reads):
  by_cigar = {}
  cigar_count = {}
  # DETERMINE MOST FREQUENT CIGAR LINE PAIR
  for (read1,read2) in reads:
    cigars = (tuple(read1.cigar),tuple(read2.cigar))
    if cigars in by_cigar:
      cigar_count[cigars]+=1
      by_cigar[cigars][0].append(read1)
      by_cigar[cigars][1].append(read2)
    else:
      cigar_count[cigars]=1
      by_cigar[cigars]=([read1],[read2])
  to_sort = map(lambda (x,y): (y,-len(str(x)),x),cigar_count.iteritems())
  to_sort.sort()
  selcigar = to_sort[-1][-1]
  reads = by_cigar[selcigar]
  del by_cigar
  del cigar_count
  del to_sort
  forward,reverse = calc_consensus(reads[0]),calc_consensus(reads[1])
  if len(reads) > 1 and (forward.is_qcfail or reverse.is_qcfail) and options.rescue_qcfail:
    new_tags1 = []
    for key,value in forward.tags:
      if key == "ZQ" and value == "Q":
        forward.is_qcfail = False
      else:
        new_tags1.append((key,value))
    new_tags2 = []
    for key,value in reverse.tags:
      if key == "ZQ" and value == "Q":
        reverse.is_qcfail = False
      else:
        new_tags2.append((key,value))
    if not forward.is_qcfail and not reverse.is_qcfail:
      forward.tags = new_tags1
      reverse.tags = new_tags2
      return forward,reverse
    else: return None,None
  elif (forward.is_qcfail or reverse.is_qcfail) and not options.include_qcfail: return None,None
  else: return forward,reverse



def get_consensus_SR(reads):
  global options
  # DETERMINE MOST FREQUENT CIGAR LINE
  by_cigar = {}
  cigar_count = {}
  for read in reads:
    tcigar = tuple(read.cigar)
    if tcigar in by_cigar:
      cigar_count[tcigar]+=1
      by_cigar[tcigar].append(read)
    else:
      cigar_count[tcigar]=1
      by_cigar[tcigar]=[read]
  to_sort = map(lambda (x,y): (y,-len(str(x)),x),cigar_count.iteritems())
  to_sort.sort()
  selcigar = to_sort[-1][-1]
  reads = by_cigar[selcigar]
  del by_cigar
  del cigar_count
  del to_sort
  consensus = calc_consensus(reads)
  if len(reads) > 1 and consensus.is_qcfail and options.rescue_qcfail:
    new_tags = []
    for key,value in consensus.tags:
      if key == "ZQ" and value == "Q":
        consensus.is_qcfail = False
      else:
        new_tags.append((key,value))
    if not consensus.is_qcfail:
      consensus.tags = new_tags
      return consensus
    else: return None
  elif consensus.is_qcfail and not options.include_qcfail: return None
  else: return consensus


if options.outprefix == "":
  sys.stderr.write("Outprefix can not be empty!\n")
  sys.exit()

if options.outdir != None and not os.path.isdir(options.outdir):
  sys.stderr.write("Output folder does not exist!\n")
  sys.exit()
elif options.outdir == None:
  options.outdir = ""
else:
  options.outdir = options.outdir.rstrip('/')+'/'

cfastq_SR, cfastq_PE = 0,0
if options.outnewconsensus != None:
  outfilename = options.outdir+options.outnewconsensus+"_SR.fastq"
  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
  fastq_SR = open(outfilename,'w')
  outfilename = options.outdir+options.outnewconsensus+"_r1.fastq"
  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
  fastq_r1 = open(outfilename,'w')
  outfilename = options.outdir+options.outnewconsensus+"_r2.fastq"
  if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
  fastq_r2 = open(outfilename,'w')

## CREATE OUTPUT FILE(S)/STREAM
fileflags = 'wb'
if options.SAM: fileflags = 'w'

files = args
if options.pipe: files = [None]
if len(files) > 1:
  files=files[:1]
  sys.stderr.write("This script supports only one input file! Limiting processing to first filename.\n")
outfile = None


sys.stderr.write("WARNING: This script will 'cluster' reads based on both outer coordinates and pick a representative sequence of the cluster as the one with the dominant CIGAR string. Other sequences are lost and will not be considered in consensus calling.\n")

for filename in files:
  if filename == None and not options.SAM:
    infile = pysam.Samfile( "-", 'rb' )
  elif filename == None and options.SAM:
    infile = pysam.Samfile( "-", 'r' )
  else:
    infile = pysam.Samfile( filename, 'rb' )

  id2lib = {}
  if options.library and 'RG' in infile.header:
    for rgroup in infile.header['RG']:
      if 'LB' in rgroup and 'ID' in rgroup:
        id2lib[rgroup['ID']] = rgroup['LB']

  if outfile == None:
    if options.verbose: sys.stderr.write("Creating output files/streams...\n")

    if options.pipe:
      outfile = pysam.Samfile( "-", fileflags, template = infile)
      if options.verbose: sys.stderr.write("BAM/SAM output on stdout...\n")
    else:
      outfilename = options.outdir+options.outprefix+".bam"
      if options.verbose: sys.stderr.write("Creating: %s\n"%outfilename)
      outfile = pysam.Samfile( outfilename , fileflags, template = infile)

    if ('HD' in outfile.header) and ('SO' in outfile.header['HD']):
      outfile.header['HD']['SO'] = 'unsorted'
    else:
      outfile.header['HD'] = {'VN': '1.4','SO':'unsorted'}

  lcheck = None,None
  variants = {}
  incomplete_variants = {}
  curpos = None
  curvariants = {}
  total_reads = 0
  out_reads = 0
  out_reads_SR = 0
  out_reads_kept = 0
  for read in infile:
    if read.qual == None: continue
    if options.fixID: read.qname = read.qname.split('/')[0]
    total_reads += 1
    if options.verbose and total_reads % 100000 == 0: 
      sys.stderr.write("Reads in %d / PCR dups out %d PE | %d SR / Unmapped out %d / FastQ realignment %d PE | %d SR ( %.2f%% ; current pos: %s)\n"%(total_reads,out_reads,out_reads_SR,out_reads_kept,cfastq_PE,cfastq_SR,(out_reads*2+out_reads_SR+out_reads_kept+cfastq_SR+cfastq_PE*2)/float(total_reads)*100,str(curpos)))

    if read.is_qcfail and not options.include_qcfail and not options.rescue_qcfail: 
      #if options.verbose: sys.stderr.write("QC FAIL\n")
      continue
    if read.is_unmapped and not options.keep: 
      #if options.verbose: sys.stderr.write("UNMAPPED\n")
      continue
    elif read.is_unmapped and options.keep: 
      if not read.is_qcfail:
        outfile.write(read)
        out_reads_kept += 1
      continue

    if not read.is_paired and options.merged and not read.qname.startswith("M_"): 
      #if options.verbose: sys.stderr.write("MERGED\n")
      continue
    RG = None
    if not options.ignore_RG:
      if read.tags != None:
        for key,value in read.tags:
          if key == "RG":
            if value in id2lib: RG = id2lib[value]
            else: RG = value
            break
    if options.UMI:
      if read.tags != None:
        for key,value in read.tags:
          if key == "XJ":
            RG = value
            break
        
    if RG not in variants: variants[RG] = {}
    if RG not in incomplete_variants: incomplete_variants[RG] = {}
    if RG not in curvariants: curvariants[RG] = {}

    if sum(map(len,variants.itervalues())) > options.rbuffer and ((lcheck[0] != curpos[0]) or (lcheck[1]+options.MaxLength < curpos[1])):
      lcheck = (curpos[0],curpos[1])
      if options.verbose: sys.stderr.write("Full buffer (%d)"%sum(map(len,variants.itervalues()))+str(curpos)+" \n")
      for cRG in variants.keys():
        hvariants = {}
        for (hchr,outpos,outpos_r2),reads in variants[cRG].iteritems():
          if ((type(hchr) != type(()) and # SINGLE CHROM MAPPING PE
                ((hchr != curpos[0]) or 
                ((hchr == curpos[0]) and (max(outpos[1],outpos_r2[1])+options.MaxLength < curpos[1])))) or
              (type(hchr) == type(()) and # CROSS CONTIG MAPPING PE
                (((hchr[0] != curpos[0]) and (hchr[1] != curpos[0])) or
                ((hchr[0] == curpos[0]) and (outpos[1]+options.MaxLength < curpos[1])) or 
                ((hchr[1] == curpos[0]) and (outpos_r2[1]+options.MaxLength < curpos[1]))))):
            forward,reverse = get_consensus(reads)
            if forward == None or reverse == None:
              #if options.verbose: sys.stderr.write("FAILED CONSENSUS\n")
              continue
            elif (forward.is_unmapped or reverse.is_unmapped) and options.outnewconsensus != None:
              cfastq_PE+=1
              seq = forward.seq
              qual = forward.qual
              if forward.is_reverse:
                seq = seq.translate(table)[::-1]
                qual = qual[::-1]
              fastq_r1.write("@%s/1\n%s\n+\n%s\n"%(forward.qname,seq,qual))
              seq = reverse.seq
              qual = reverse.qual
              if reverse.is_reverse:
                seq = seq.translate(table)[::-1]
                qual = qual[::-1]
              fastq_r2.write("@%s/2\n%s\n+\n%s\n"%(reverse.qname,seq,qual))
            else:
              outfile.write(forward)
              outfile.write(reverse)
              out_reads += 1
          else:
            hvariants[(hchr,outpos,outpos_r2)]=reads
        if (len(hvariants) > 0) or (RG == cRG): variants[cRG] = hvariants
        else: del variants[cRG]
      if options.verbose: sys.stderr.write("- Full buffer (%d)"%sum(map(len,variants.itervalues()))+str(curpos)+" \n")

    if read.is_paired: # PE DATA
      if read.mate_is_unmapped and options.keep:
        outfile.write(read)
        out_reads_kept += 1
        continue
      #else:
        #if options.verbose: sys.stderr.write("UNMAPPED\n")

      curpos = (read.tid,read.pos)
      hchr = read.tid
      outpos = (read.tid,read.pos)
      if read.is_reverse: outpos = (read.tid,read.pos+aln_length(read.cigar))

      if  read.is_read1: #FORWARD READ
        if read.qname not in incomplete_variants[RG]:
          incomplete_variants[RG][read.qname] = [read,outpos]
        else:
          read_r2,outpos_r2 = incomplete_variants[RG][read.qname]
          del incomplete_variants[RG][read.qname]
          if outpos_r2[0] != hchr: hchr = hchr,outpos_r2[0]
          if (hchr,outpos,outpos_r2) not in variants[RG]:
            variants[RG][(hchr,outpos,outpos_r2)] = [(read,read_r2)]
          else:
            variants[RG][(hchr,outpos,outpos_r2)].append((read,read_r2))
      elif read.is_read2:  #REVERSE READ
        if read.qname not in incomplete_variants[RG]:
          incomplete_variants[RG][read.qname] = [read,outpos]
        else:
          read_r1,outpos_r1 = incomplete_variants[RG][read.qname]
          del incomplete_variants[RG][read.qname]
          if outpos_r1[0] != hchr: hchr = outpos_r1[0],hchr
          if (hchr,outpos_r1,outpos) not in variants[RG]:
            variants[RG][(hchr,outpos_r1,outpos)] = [(read_r1,read)]
          else:
            variants[RG][(hchr,outpos_r1,outpos)].append((read_r1,read))
      else:
        sys.stderr.write("Should not happen!")
    else: # SR DATA
      if (curpos != None) and ((read.tid,read.pos) != curpos):
        if options.fivePrime and (read.tid == curpos[0]):
          hpos = read.pos-options.MaxLength
          hvariants = {}
          for key,value in curvariants[RG].iteritems():
            if (key[1] < hpos):
              cread = get_consensus_SR(value[0])
              if cread == None: 
                #if options.verbose: sys.stderr.write("FAILED CONSENSUS\n")
                continue
              elif cread.is_unmpapped and options.outnewconsensus != None:
                cfastq_SR+=1
                seq = cread.seq
                qual = cread.qual
                if cread.is_reverse:
                  seq = seq.translate(table)[::-1]
                  qual = qual[::-1]
                fastq_SR.write("@%s\n%s\n+\n%s\n"%(cread.qname,seq,qual))
              else:
                outfile.write(cread)
                out_reads_SR += 1
            else:
              hvariants[key]=value
          curvariants[RG] = hvariants
        else:
          for key,value in curvariants[RG].iteritems():
            cread = get_consensus_SR(value[0])
            if cread == None: 
              #if options.verbose: sys.stderr.write("FAILED CONSENSUS\n")
              continue
            elif cread.is_unmapped and options.outnewconsensus != None:
              cfastq_SR+=1
              seq = cread.seq
              qual = cread.qual
              if cread.is_reverse:
                seq = seq.translate(table)[::-1]
                qual = qual[::-1]
              fastq_SR.write("@%s\n%s\n+\n%s\n"%(cread.qname,seq,qual))
            else:
              outfile.write(cread)
              out_reads_SR += 1
          curvariants[RG] = {}
      curpos = (read.tid,read.pos)
      strand = read.is_reverse
      outpos = curpos[1]
      if strand and options.fivePrime: outpos+=aln_length(read.cigar)

      nkey = (strand,outpos)
      if not options.fivePrime: nkey = (strand,outpos,aln_length(read.cigar))

      if nkey in curvariants[RG]:
        curvariants[RG][nkey][0].append(read)
        curvariants[RG][nkey][1]+=1
      else:
        curvariants[RG][nkey] = [[read],1]

  for RG in curvariants.keys():
    for key,value in curvariants[RG].iteritems():
      read = get_consensus_SR(value[0])
      if read == None: 
        #if options.verbose: sys.stderr.write("FAILED CONSENSUS\n")
        continue
      elif read.is_unmapped and options.outnewconsensus != None:
        cfastq_SR+=1
        seq = read.seq
        qual = read.qual
        if read.is_reverse:
          seq = seq.translate(table)[::-1]
          qual = qual[::-1]
        fastq_SR.write("@%s\n%s\n+\n%s\n"%(read.qname,seq,qual))
      else:
        outfile.write(read)
        out_reads_SR += 1
    del curvariants[RG]

  for RG in variants.keys():
    for key,value in variants[RG].iteritems():
      read1,read2 = get_consensus(value)
      if read1 == None or read2 == None: 
        #if options.verbose: sys.stderr.write("FAILED CONSENSUS\n")
        continue
      elif (read1.is_unmapped or read2.is_unmapped) and options.outnewconsensus != None:
        cfastq_PE+=1
        seq = read1.seq
        qual = read1.qual
        if read1.is_reverse:
          seq = seq.translate(table)[::-1]
          qual = qual[::-1]
        fastq_r1.write("@%s/1\n%s\n+\n%s\n"%(read1.qname,seq,qual))
        seq = read2.seq
        qual = read2.qual
        if read2.is_reverse:
          seq = seq.translate(table)[::-1]
          qual = qual[::-1]
        fastq_r2.write("@%s/2\n%s\n+\n%s\n"%(read2.qname,seq,qual))
      else:
        outfile.write(read1)
        outfile.write(read2)
        out_reads += 1
    del variants[RG]

if options.outnewconsensus != None:
  fastq_SR.close()
  fastq_r1.close()
  fastq_r2.close()
  if cfastq_PE == 0:
    outfilename = options.outdir+options.outnewconsensus+"_r1.fastq"
    if options.verbose: sys.stderr.write("Removing empty file: %s\n"%outfilename)
    os.remove(outfilename)
    outfilename = options.outdir+options.outnewconsensus+"_r2.fastq"
    if options.verbose: sys.stderr.write("Removing empty file: %s\n"%outfilename)
    os.remove(outfilename)
  if cfastq_SR == 0:
    outfilename = options.outdir+options.outnewconsensus+"_SR.fastq"
    if options.verbose: sys.stderr.write("Removing empty file: %s\n"%outfilename)
    os.remove(outfilename)

if options.verbose: 
  sys.stderr.write("Total reads in %d / PCR dups out %d PE | %d SR / Unmapped out %d / FastQ realignment %d PE | %d SR ( %.2f%%)\n"%(total_reads,out_reads,out_reads_SR,out_reads_kept,cfastq_PE,cfastq_SR,(out_reads*2+out_reads_SR+out_reads_kept+cfastq_SR+cfastq_PE*2)/float(total_reads)*100))
