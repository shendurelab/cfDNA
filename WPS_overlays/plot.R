args <- commandArgs(trailingOnly = TRUE)
# s1 <- "allSamples"
# factor <- "TSS"
s1 <- args[1]
factor <- args[2]

png(sprintf("plots/%s_%s.png",factor,s1),width=2000,height=1000,type="cairo")
par(mfrow=c(2,2))

data <- read.table(sprintf("short/%s_%s.tsv",s1,factor),as.is=T,sep="\t",header=T,comment.char="~")
newData <- data.frame(position=data[,1]-(dim(data)[1]/2),sites=rowSums(data[,c(2:11)]))
for (i in 1:10) {
  newData[[sprintf("sites.%d",i)]] <- newData$sites-data[,i+1]
}

for (i in seq(2,12,1)) {
  helper <- abs(mean(newData[,i]))
  newData[,i] <- newData[,i]/helper
}

# msites <- median(c(newData$sites[c(1:500)],newData$sites[c(length(newData$sites)-500:length(newData$sites))]))
msites <- 0
matplot(newData[,1],newData[,2:12]-msites,type="l",lty=c(1,2),lwd=c(2,rep(1,10)),col="black",main=sprintf("%s: short",factor),xlab="Position",ylab="Normalized sum of counts",xlim=c(-500,500))
matplot(newData[,1],newData[,2:12]-msites,type="l",lty=c(1,2),lwd=c(2,rep(1,10)),col="black",main=sprintf("%s: short",factor),xlab="Position",ylab="Normalized sum of counts",xlim=c(-2000,2000))

data <- read.table(sprintf("long/%s_%s.tsv",s1,factor),as.is=T,sep="\t",header=T,comment.char="~")
newData <- data.frame(position=data[,1]-(dim(data)[1]/2),sites=rowSums(data[,c(2:11)]))
for (i in 1:10) {
  newData[[sprintf("sites.%d",i)]] <- newData$sites-data[,i+1]
}

for (i in seq(2,12,1)) {
  helper <- abs(mean(newData[,i]))
  newData[,i] <- newData[,i]/helper
}

# msites <- median(c(newData$sites[c(1:500)],newData$sites[c(length(newData$sites)-500:length(newData$sites))]))
msites <- 0

matplot(newData[,1],newData[,2:12]-msites,type="l",lty=c(1,2),lwd=c(2,rep(1,10)),col="black",main=sprintf("%s: long",factor),xlab="Position",ylab="Normalized sum of counts",xlim=c(-500,500))
matplot(newData[,1],newData[,2:12]-msites,type="l",lty=c(1,2),lwd=c(2,rep(1,10)),col="black",main=sprintf("%s: long",factor),xlab="Position",ylab="Normalized sum of counts",xlim=c(-2000,2000))

dev.off()

