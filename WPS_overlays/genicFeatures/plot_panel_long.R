args <- commandArgs(trailingOnly = TRUE)
# s1 <- "bulk-cfDNA"
# s2 <- "bulk-cfDNA_2mer_sim"
s1 <- args[1]
s2 <- args[2]
factors <- c("TSS","TSE","StartCodon","StopCodon","exonSpliceAcceptor","exonSpliceDonor")

library(caTools)

png(sprintf("plots/%s_%s_Vs_%s.png","panel",s1,s2),width=1200,height=1500,type="cairo")
par(mfrow=c(3,2))
for (factor in factors)
{
  data <- read.table(sprintf("long/%s_%s.tsv",s1,factor),as.is=T,sep="\t",header=T,comment.char="~")
  data2 <- read.table(sprintf("long/%s_%s.tsv",s2,factor),as.is=T,sep="\t",header=T,comment.char="~")
  newData <- data.frame(position=data[,1]-(dim(data)[1]/2),sites=rowSums(data[,c(2:11)]))
  newData2 <- data.frame(position=data2[,1]-(dim(data2)[1]/2),sites=rowSums(data2[,c(2:11)]))
  for (i in 1:10) {
    newData[[sprintf("sites.%d",i)]] <- newData$sites-data[,i+1]
    newData2[[sprintf("sites.%d",i)]] <- newData2$sites-data2[,i+1]
  }

  for (i in seq(2,12,1)) {
    helper <- abs(mean(c(newData[c(1:500),i],newData[c(4500:5000),i])))
    newData[,i] <- newData[,i]/helper-runmean(newData[,i]/helper,200)
    helper <- abs(mean(c(newData2[c(1:500),i],newData2[c(4500:5000),i])))
    newData2[,i] <- newData2[,i]/helper-runmean(newData2[,i]/helper,200)
  }

  matplot(newData[,1],newData[,2:12]-newData2[,2:12],type="l",lty=c(1,2),lwd=c(2,rep(1,10)),col="black",main=sprintf("%s: long fraction",factor),xlab="Position",ylab="Data - Control",xlim=c(-2000,2000),cex=3.0)
  abline(v=0,lty=2)
  
dev.off()

