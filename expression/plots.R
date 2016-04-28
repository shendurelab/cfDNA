samples <- c("BH01","IH01","IH02")

##################################

proteinAtlas <- read.table("protein_atlas/RNAtable.tsv",header=T,as.is=T,sep="\t")
rownames(proteinAtlas) <- proteinAtlas$GeneID

ndata <- proteinAtlas[,-1]
ndata[ndata == 0] <- NA
# at least 3 tissues with non-zero values
ndata <- ndata[apply(ndata,1,FUN=function(x) { length(which(is.na(x))) < length(x)-2 } ),] 
ndata[is.na(ndata)] <- 0.04
logndata <- log2(ndata)
dim(logndata)

##################################

tLabels <- read.table("labels.txt",header=T,as.is=T,sep="\t",quote="\"")

fftColumns <- 29:52 # 160-222
selFreq <- c("193","196","199")

library(gplots)

##################################

  pdf("allFreq_correlation_plot.pdf",width=8,height=6)
  for (sample in samples)
  {
    fdata <- read.table(sprintf("fft_summaries/fft_%s_WPS.tsv.gz",sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]

    res <- cor(fdata[,fftColumns],logndata2,use="pairwise.complete.obs")

    matplot(as.numeric(sub("X","",names(fdata[,fftColumns]))),res,type="l",xlab="Frequency",ylab="Correlation",col="darkgrey",lwd=1,lty=1,main=sprintf("%s: Correlation of intensities across tissues",sample))
    lines(as.numeric(sub("X","",names(fdata[,fftColumns]))),cor(fdata[,fftColumns],logndata2[,"NB.4"],method="pearson",use="pairwise.complete.obs"),col="black",lwd=2,type="b",pch=19,cex=0.6)
    legend("topright","NB-4",fill="black")
  }
  dev.off()

##################################

  pdf("Ave193-199bp_correlation.pdf",width=8,height=15)
  for (sample in samples)
  {
    fdata <- read.table(sprintf("fft_summaries/fft_%s_WPS.tsv.gz",sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]

    res <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")
    res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res))
    textplot(res[order(res$correlation),])
    title(sample)
  }
  dev.off()

##################################

  # Replace BH01 with sample you are using as reference in the rank comparison
  fdata <- read.table(sprintf("fft_summaries/fft_%s_WPS.tsv.gz","BH01"),as.is=T,sep="\t",header=T,comment.char="~")
  colnames(fdata) <- sub("X","",colnames(fdata))
  rownames(fdata) <- fdata[,1]
  fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
  logndata2 <- logndata[fdata[,1],]
  refCorrelation <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")

  pdf("Ave193-199bp_correlation_rank.pdf",width=10,height=15)
  for (sample in samples)
  {
    fdata <- read.table(sprintf("fft_summaries/fft_%s_WPS.tsv.gz",sample),as.is=T,sep="\t",header=T,comment.char="~")
    colnames(fdata) <- sub("X","",colnames(fdata))
    rownames(fdata) <- fdata[,1]
    fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
    logndata2 <- logndata[fdata[,1],]

    res <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")
    res <- data.frame(category=tLabels$Category,description=tLabels$Type,tissue=colnames(res),correlation=as.numeric(res),rankDiff=rank(refCorrelation)-rank(res))
    par(mfrow=c(2,1))
    textplot(head(res[order(res$rankDiff,decreasing=T),],15))
    title(sprintf("By correlation rank difference: %s (vs. normal)",sample))
    textplot(head(res[order(res$correlation),],15))
    title(sprintf("By correlation: %s",sample))
  }
  dev.off()

