Analysis of 167 bp fragments
---------------------------

We extracted fragments with insert sizes of exactly 167 bp and, for single-stranded library preparation, accounted for asymmetric adapter ligation bias by considering the strand of read 1. We used `BAM2FragmentationPatterns.py` and added 50 bp padding of genomic context on either end of the fragment during dinucleotide frequency calculation.  We performed this analysis for both real and simulated datasets, and for both SSP and DSP libraries.

The results were then visualized in R:

```
realdata <- read.table("nucfreq.SSP_read1.txt", header=F, skip=1,fill=T)
colnames(realdata) <- c("type", "nuc", paste0("U", as.character(seq(50, 1, -1))), paste0("D", as.character(seq(1, 216, 1))))
simdata <- read.table("nucfreq.SSP_sim_read1.txt", header=F, skip=1,fill=T)
colnames(simdata) <- colnames(realdata)

orient <- "ThreePrime"
ATstuff <- c()
GCstuff <- c()
for (mycol in seq(3, ncol(realdata))) {
    TA <- realdata[realdata$nuc == "TA" & realdata$type==orient,mycol]
    TT <- realdata[realdata$nuc == "TT" & realdata$type==orient,mycol]
    AT <- realdata[realdata$nuc == "AT" & realdata$type==orient,mycol]
    AA <- realdata[realdata$nuc == "AA" & realdata$type==orient,mycol]
    mysum <- TA+TT+AT+AA
    ATstuff <- c(ATstuff, mysum)

    CC <- realdata[realdata$nuc == "CC" & realdata$type==orient,mycol]
    CG <- realdata[realdata$nuc == "CG" & realdata$type==orient,mycol]
    GC <- realdata[realdata$nuc == "GC" & realdata$type==orient,mycol]
    GG <- realdata[realdata$nuc == "GG" & realdata$type==orient,mycol]
    mysum <- CC+CG+GC+GG
    GCstuff <- c(GCstuff, mysum)
}

simATstuff <- c()
simGCstuff <- c()
for (mycol in seq(3, ncol(simdata))) {
    TA <- simdata[simdata$nuc == "TA" & simdata$type==orient,mycol]
    TT <- simdata[simdata$nuc == "TT" & simdata$type==orient,mycol]
    AT <- simdata[simdata$nuc == "AT" & simdata$type==orient,mycol]
    AA <- simdata[simdata$nuc == "AA" & simdata$type==orient,mycol]
    mysum <- TA+TT+AT+AA
    simATstuff <- c(simATstuff, mysum)

    CC <- simdata[simdata$nuc == "CC" & simdata$type==orient,mycol]
    CG <- simdata[simdata$nuc == "CG" & simdata$type==orient,mycol]
    GC <- simdata[simdata$nuc == "GC" & simdata$type==orient,mycol]
    GG <- simdata[simdata$nuc == "GG" & simdata$type==orient,mycol]
    mysum <- CC+CG+GC+GG
    simGCstuff <- c(simGCstuff, mysum)
}

myATratio <- log2(ATstuff / simATstuff)
myGCratio <- log2(GCstuff / simGCstuff)
mydata <- data.frame(pos=seq(-50,215,1), ATratio=myATratio, GCratio=myGCratio)
ggplot(mydata, aes(pos, ATratio)) + geom_line(size=1.5) + geom_line(aes(pos, GCratio), size=1.5) +
scale_x_continuous(breaks=seq(-37,216,20), labels=as.character(seq(-37,216,20)-83)) + geom_vline(xintercept=c(0,83,167), linetype="dashed", size=c(1,1.4,1), color="gray50") +
ylim(-0.5, 0.8)
```