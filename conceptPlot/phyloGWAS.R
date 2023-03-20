library(ggplot2)
library(cowplot)
library(wesanderson)

# make function to read in files and store pvals and rho values
readFiles<-function(fn,i) {
  data<-read.table(paste("~/projects/ILSSims/ILSsims/phyloGWASsims/pGWAS.",fn,".",as.character(i),".txt",sep=""))
  return(data.frame(b=data$V1,val="null"))
}

getTrueNsubs<-function(fn, i) {
  d<-scan(paste("~/projects/ILSSims/ILSsims/phyloGWASsims/pGWAS.",fn,".",as.character(i),".txt",sep=""),nlines=1,what="character", skip=3)
  return (as.numeric(d[5]))
}
