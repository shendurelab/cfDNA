Spacing of peaks flanking DHS sites
------------------------------------------

DHS peaks for 349 primary tissue and cell line samples were downloaded from http://www.uwencode.org/proj/Science_Maurano_Humbert_et_al/data/all_fdr0.05_hot.tgz. Samples derived from fetal tissues (233 of 349) were removed from the analysis as they behaved inconsistently within tissue type, possibly because of unequal representation of cell types within each tissue sample. 116 DHS callsets from a variety of cell lineages were retained for analysis. 

We first transformed the peak score for each DHS to percentiles within samples.  Many of the lowest scoring DHSs could in principle be filtered before downstream analysis, but in practice this did not seem to substantially affect results. 

```
% for bed in annots/Maurano.dhs/*.bed; echo "python mkPercentiles.py $bed > annots/Maurano.dhs.with.percentile/`basename $bed`" >> percentile.jobs.sh
```

For the midpoint of each DHS peak in a particular set, the nearest upstream and downstream calls in the CH01 callset were identified, and the distance between the centers of those two calls was calculated. 

```
% for peakcalls in nucleosomeCalls*bed.gz; mkdir dhs/`basename $peakcalls .bed.gz | sed s/nucleosomeCalls_//`

% ( for bed in annots/Maurano.dhs.with.percentile/*.bed; for peakcalls in nucleosomeCalls*bed.gz; echo "python DHSdistances.py $bed $peakcalls > dhs.intervals.by.sample/`basename $peakcalls .bed.gz | sed s/nucleosomeCalls_//`/`basename $bed .bed`.intervals.txt" ) > intervals.jobs.sh
```

The distribution of all such distances was visualized for each DHS peak callset using a smoothed density estimate calculated for distances between 0 and 500 bp.

```
everything <- data.frame()
for (dhs.set in dhs.sets) {
  for (samplename in mysamples) {
    tmpdata <- read.table(paste0("dhs.intervals.by.sample/", samplename, "/", dhs.set), header=F)
    colnames(tmpdata) <- c("chrom", "start", "end", "minneg", "minpos", "totaldistance", "perc")
    tmpdata$dhs <- dhs.set
    everything <- rbind(everything, tmpdata)	
  }
}

everything$dhs <- gsub(".hg19.hotspot.twopass.fdr0.05.merge.intervals.txt", "", everything$dhs)
everything$dhs <- gsub("hotspot.twopass.fdr0.05.merge.intervals.txt", "", everything$dhs)
everything$dhs <- gsub("\\.", "", everything$dhs)

dhs.selected <- c("CD3_CordBlood-DS17706",
		 "CD3-DS17198",
		 "CD8-DS17203",
		 "CD14-DS18065",
		 "CD19-DS17281",
		 "CD20-DS17541",
		 "hTH1-DS17592"
)

dhs.all <- unique(everything$dhs)
dhs.not <- dhs.all[dhs.all %in% dhs.selected == FALSE]
everything$dhsfac <- as.factor(everything$dhs)
everything$dhsfac <- factor(everything$dhsfac, levels=c(dhs.not, dhs.selected))

ggplot(everything, aes(totaldistance, colour=dhsfac)) + geom_line(stat="density") +
     scale_color_manual(values=c(rep("gray50", 104), colorRampPalette(brewer.pal(12,"Paired"))(12))) +
     scale_x_continuous(breaks=seq(0,500, 50), limits=c(0,500))
```
