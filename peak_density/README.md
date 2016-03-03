Peak density within A/B compartments
------------------------------------

We computed median peak density in all 100kb bins genome-wide with `makePeakDensityBins.py`, as a BED file (available here as `median.peak.density.100kb.CH01.txt`).

We downloaded A/B compartment calls from GSE63525 for GM12878, provided with minimum resolution of 100kb. We split each call into 100kb subwindows:

```
chr1    1700000 1800000 A1
chr1    1800000 1900000 A1
chr1    1900000 2000000 B1
chr1    2000000 2100000 B1
chr1    2100000 2200000 A1
chr1    2200000 2300000 A1
```

We then visualized the results for the autosomes for compartments A1, A2, B1, B2, and B3:
```
pmed <- read.table("median.peak.density.100kb.CH01.txt", header=F)
colnames(pmed) <- c("chrom", "start", "stop", "name", "med")
pmed <- subset(pmed, pmed$chrom != "chrY" & pmed$chrom != "chrX")
comps <- read.table("GSE63525_GM12878_subcompartments.100kb.sorted.bed", header=F)
colnames(comps) <- c("chrom", "start", "stop", "compartment")
pmed$compartment = comps$compartment
ggplot(pmed, aes(med)) + geom_density(..counts.., aes(color=compartment))
```

For visualization of each chromosome separately as an ideogram, we first split the subcompartment bed file by chromosome, then plotted as follows: 

```
peak.meds <- read.table("median.peak.density.100kb.txt", header=F)
colnames(peak.meds) <- c("chrom", "start", "stop", "name", "med")
for (mychrom in paste0("chr", seq(1, 22, 1))) {
    peak.meds.sub <- subset(peak.meds, peak.meds$chrom == mychrom)
    comps.sub <- read.table(paste0("GSE63525_GM12878_subcompartments.100kb.sorted.", mychrom, ".bed"), header=F)
    colnames(comps.sub) <- c("chrom", "start", "stop", "compartment")
    comps.sub$med <- peak.meds.sub$med

    ggplot(peak.meds.sub, aes(start, med)) + 
      geom_point(size=0.7) + xlab("") + ylab("") +
      geom_rect(data=comps.sub, aes(xmax=stop, xmin=start, ymax=170, ymin=155, fill=compartment)) + 
      scale_fill_manual(values=c("#B2182B", "#D6604D", "#2166AC", "#4393C3", "#92C5DE","black")) +
      scale_y_continuous(limits=c(155,195), breaks=c(175, 185, 195)) + scale_x_continuous(limits=c(0,249250621))
   ggsave(paste0("ABplots/", mychrom, ".pdf"), width=mywidth, height=myheight)
}
```

Peak density by GC content
--------------------------

We computed average GC content in 10kb windows genome wide for hg19 using `bedtools nuc`. We calculated the number of peaks in CH01 within each 10kb bin, and plotted the relationship between these two variables, subsetting only on regions of high mappability:  
```
gccontent <- read.table("hg19.10kbBins.GCcontent.bed"(
colnames(gccontent) <- c("chrom", "start", "stop", "interval", "atperc", "gcperc", "acount", "ccount", "gcount", "tcount", "ncount", "othercount", "seqlen")
gccontent$chromcat <- "autosome"
gccontent[gccontent$chrom=="chrY", "chromcat"] <- "chrY"
gccontent[gccontent$chrom=="chrX", "chromcat"] <- "chrX"
combinedpeak <- read.table("CH01.peakdensity.bed", header=F)
gccontent$combined <- combinedpeak$V5
bedmap <- read.table("hg19.10kbBins.mappability.bed", header=F)
gccontent$map <- bedmap$V5
ggplot(subset(gccontent, gccontent$ncount < 1 & gccontent$map >= 0.5), aes(gcperc,combined)) + geom_point(aes(colour=chromcat), alpha=0.5) +
  stat_smooth(method="lm", size=1.3) + scale_x_continuous(breaks=c(0.3, 0.4, 0.5, 0.6, 0.7), labels=c("30%", "40%", "50%", "60%", "70%"))
```
