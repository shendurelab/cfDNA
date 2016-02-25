args <- commandArgs(trailingOnly = TRUE)
# s1 <- "allNormal"
# s2 <- "allNormal_2mer_sim"
# factor <- "TSS"
# factor <- "CTCF_Core-2"
s1 <- args[1]
s2 <- args[2]
factor <- args[3]

library(caTools)

selColors <- rev(c("#b2182b","#ef8a62","#fddbc7","#67a9cf","#2166ac"))

png(sprintf("plots/%s_%s_Vs_%s.png",factor,s1,s2),width=2500,height=1000,type="cairo")
par(mfrow=c(2,3))

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
  for (i in 1:10) {
    newData[[sprintf("sites.%d-%d",i,bin)]] <- newData[[sprintf("sites.all-%d",bin)]]-data[,i+1]
    newData2[[sprintf("sites.%d-%d",i,bin)]] <- newData2[[sprintf("sites.all-%d",bin)]]-data2[,i+1]
  }
}
  
for (i in seq(2,dim(newData)[2],1)) {
  helper <- abs(mean(c(newData[c(1:500),i],newData[c(4500:5000),i])))
  newData[,i] <- newData[,i]/helper
  helper <- abs(mean(c(newData2[c(1:500),i],newData2[c(4500:5000),i])))
  newData2[,i] <- newData2[,i]/helper
}

matplot(newData[,1],newData[,-1],type="l",lty=c(1,rep(2,10)),lwd=c(2,rep(1,10)),col=rep(selColors,each=11),main=sprintf("[Data] %s: short",factor),xlab="Position",ylab="Normalized sum of counts",xlim=c(-2000,2000),ylim=c(min(newData[,-1],newData2[,-1]),max(newData[,-1],newData2[,-1])))
abline(v=0,lty=2)

matplot(newData[,1],newData2[,-1],type="l",lty=c(1,rep(2,10)),lwd=c(2,rep(1,10)),col=rep(selColors,each=11),main=sprintf("[Control] %s: short",factor),xlab="Position",ylab="Normalized sum of counts",xlim=c(-2000,2000),ylim=c(min(newData[,-1],newData2[,-1]),max(newData[,-1],newData2[,-1])))
abline(v=0,lty=2)

matplot(newData[,1],(newData[,-1]-newData2[,-1]),type="l",lty=c(1,rep(2,10)),lwd=c(2,rep(1,10)),col=rep(selColors,each=11),main=sprintf("[Diff] %s: short",factor),xlab="Position",ylab="Data - Control",xlim=c(-2000,2000))
abline(v=0,lty=2)
# -apply(newData[,-1]-newData2[,-1],2,FUN=function(x)runmean(x,200))

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
  for (i in 1:10) {
    newData[[sprintf("sites.%d-%d",i,bin)]] <- newData[[sprintf("sites.all-%d",bin)]]-data[,i+1]
    newData2[[sprintf("sites.%d-%d",i,bin)]] <- newData2[[sprintf("sites.all-%d",bin)]]-data2[,i+1]
  }
}
  
for (i in seq(2,dim(newData)[2],1)) {
  helper <- abs(mean(c(newData[c(1:500),i],newData[c(4500:5000),i])))
  newData[,i] <- newData[,i]/helper-runmean(newData[,i]/helper,200)
  helper <- abs(mean(c(newData2[c(1:500),i],newData2[c(4500:5000),i])))
  newData2[,i] <- newData2[,i]/helper-runmean(newData2[,i]/helper,200)
}

matplot(newData[,1],newData[,-1],type="l",lty=c(1,rep(2,10)),lwd=c(2,rep(1,10)),col=rep(selColors,each=11),main=sprintf("[Data] %s: long",factor),xlab="Position",ylab="Normalized sum of counts [runMeanAdj]",xlim=c(-2000,2000),ylim=c(min(newData[,-1],newData2[,-1]),max(newData[,-1],newData2[,-1])))
abline(v=0,lty=2)
matplot(newData[,1],newData2[,-1],type="l",lty=c(1,rep(2,10)),lwd=c(2,rep(1,10)),col=rep(selColors,each=11),main=sprintf("[Control] %s: long",factor),xlab="Position",ylab="Normalized sum of counts [runMeanAdj]",xlim=c(-2000,2000),ylim=c(min(newData[,-1],newData2[,-1]),max(newData[,-1],newData2[,-1])))
abline(v=0,lty=2)
matplot(newData[,1],(newData[,-1]-newData2[,-1]),type="l",lty=c(1,rep(2,10)),lwd=c(2,rep(1,10)),col=rep(selColors,each=11),main=sprintf("[Diff] %s: long",factor),xlab="Position",ylab="Data - Control",xlim=c(-2000,2000))
abline(v=0,lty=2)
# -apply(newData[,-1]-newData2[,-1],2,FUN=function(x)runmean(x,200))
dev.off()
