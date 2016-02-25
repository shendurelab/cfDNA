args <- commandArgs(trailingOnly = TRUE)
# s1 <- "allNormal"
# s2 <- "allNormal_2mer_sim"
s1 <- args[1]
s2 <- args[2]
factors <- c("TSS","TSE","StartCodon","StopCodon","ExonSpliceAcceptor","ExonSpliceDonor")

library(caTools)

selColors <- rev(c("#b2182b","#ef8a62","#fddbc7","#67a9cf","#2166ac"))

pdf(sprintf("plots/%s_%s_Vs_%s.pdf","all",s1,s2),height=22,width=17)
layout(matrix(c(1,2,4,6,8,10,12,1,3,5,7,9,11,13),ncol=2),heights=c(1.2,4,4,4,4,4,4))
par(mar=c(4, 4, 2, 2) + 0.1)
plot.new()
text(0.5,0.5,"Accumulated window protection scores at features of protein-coding genes by NB-4 expression bins",cex=1.8,font=1.5)
for (factor in factors)
{
  bin <- 1
  data <- read.table(sprintf("short/%s_%s_bin%d.tsv",s1,factor,bin),as.is=T,sep="\t",header=T,comment.char="~")
  newData <- data.frame(position=data[,1]-(dim(data)[1]/2))
  newData2 <- data.frame(position=data[,1]-(dim(data)[1]/2))

  for (bin in 1:5)
  {
    data <- read.table(sprintf("short/%s_%s_bin%d.tsv",s1,factor,bin),as.is=T,sep="\t",header=T,comment.char="~")
    data2 <- read.table(sprintf("short/%s_%s_bin%d.tsv",s2,factor,bin),as.is=T,sep="\t",header=T,comment.char="~")
    newData[[sprintf("sites.all-%d",bin)]] <- rowSums(data[,c(2:11)])
    newData2[[sprintf("sites.all-%d",bin)]] <- rowSums(data2[,c(2:11)])
  }
    
  for (i in seq(2,dim(newData)[2],1)) {
    helper <- abs(mean(c(newData[c(1:500),i],newData[c(4500:5000),i])))
    newData[,i] <- newData[,i]/helper
    helper <- abs(mean(c(newData2[c(1:500),i],newData2[c(4500:5000),i])))
    newData2[,i] <- newData2[,i]/helper
  }

  yLims = c(min(newData[,-1]-newData2[,-1]),max(newData[,-1]-newData2[,-1]))
  matplot(newData[,1],newData[,-1]-newData2[,-1],type="l",lty=1,lwd=2,col=selColors,main=sprintf("35-80bp fraction: %s",factor),xlab="Position",ylab="Data - Control",xlim=c(-2000,2000),cex=3.0,ylim=yLims,las=1)
  abline(v=c(-1750,-1500,-1250,-1000,-750,-500,-250,0,250,500,750,1000,1250,1500,1750),lty=2,lwd=c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,1))
  legend("topright",c("Low","20-40%","Medium","60-80%","High"),fill=selColors,bty='n',border=selColors)

  bin <- 1
  data <- read.table(sprintf("long/%s_%s_bin%d.tsv",s1,factor,bin),as.is=T,sep="\t",header=T,comment.char="~")
  newData <- data.frame(position=data[,1]-(dim(data)[1]/2))
  newData2 <- data.frame(position=data[,1]-(dim(data)[1]/2))

  for (bin in 1:5)
  {
    data <- read.table(sprintf("long/%s_%s_bin%d.tsv",s1,factor,bin),as.is=T,sep="\t",header=T,comment.char="~")
    data2 <- read.table(sprintf("long/%s_%s_bin%d.tsv",s2,factor,bin),as.is=T,sep="\t",header=T,comment.char="~")
    newData[[sprintf("sites.all-%d",bin)]] <- rowSums(data[,c(2:11)])
    newData2[[sprintf("sites.all-%d",bin)]] <- rowSums(data2[,c(2:11)])
  }
    
  for (i in seq(2,dim(newData)[2],1)) {
    helper <- abs(mean(c(newData[c(1:500),i],newData[c(4500:5000),i])))
    newData[,i] <- newData[,i]/helper-runmean(newData[,i]/helper,200)
    helper <- abs(mean(c(newData2[c(1:500),i],newData2[c(4500:5000),i])))
    newData2[,i] <- newData2[,i]/helper-runmean(newData2[,i]/helper,200)
  }

  yLims = c(min(newData[,-1]-newData2[,-1]),max(newData[,-1]-newData2[,-1]))
  matplot(newData[,1],newData[,-1]-newData2[,-1],type="l",lty=1,lwd=2,col=selColors,main=sprintf("120-180bp fraction: %s",factor),xlab="Position",ylab="Data - Control",xlim=c(-2000,2000),cex=3.0,ylim=yLims,las=1)
  abline(v=c(-1750,-1500,-1250,-1000,-750,-500,-250,0,250,500,750,1000,1250,1500,1750),lty=2,lwd=c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,1))
  legend("topright",rev(c("Low","20-40%","Medium","60-80%","High")),fill=rev(selColors),bty='n',border=rev(selColors))
}

dev.off()

