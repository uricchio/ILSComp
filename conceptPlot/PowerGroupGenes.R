library(wesanderson)
library(ggplot2)
library(cowplot)
# ILSsims/simData/simGroup/highDiscP.nLin500.nG1.txt

readData<-function(type,nLin,nG) {
  f<-paste("~/projects/ILSsims/ILSsims/simData/simGroup/",type,"DiscP.nLin",nLin,".nG",nG,".txt",sep="")
  data<-read.table(f)
  tot =sum(data$V2)
  #print(tot)
  #print(sum(data$V1))
  pow = sum(data$V1)/tot
  pow2 = sum(data$V1 > 0)/length(data$V1)
  return(data.frame(pow=pow,disc=type,nLin=nLin,nG=nG,pow2=pow2))
}

allDat<-data.frame()

for (nLin in c(10,20,50,100,200,500) ) {
  for (nG in c(1,2,5,10,20)) {
   allDat<-rbind(allDat,readData("high",nLin,nG))
   #allDat<-rbind(allDat,readData("mid",nLin,nG))
   #allDat<-rbind(allDat,readData("low",nLin,nG))
  }
}

pal <- wes_palette("Darjeeling1", 6, type = "continuous")

plA<-ggplot(allDat)+geom_point(aes(nG,pow,color=as.factor(nLin)),size=2.2)+
     theme_classic()+scale_x_log10(breaks=c(1,2,5,10,20),name="Number of causal genes")+
    geom_line(data=allDat,aes(nG,pow,color=as.factor(nLin)),size=0.6)+
    scale_color_manual(values=pal,name="Number of lineages")+ylab(expression("Power (" * alpha * "=5e-3)"))+theme(legend.position="NA")

plB<-ggplot(allDat)+geom_point(aes(nG,pow2,color=as.factor(nLin)),size=2.2)+
  theme_classic()+scale_x_log10(breaks=c(1,2,5,10,20),name="Number of causal genes")+
  geom_line(data=allDat,aes(nG,pow2,color=as.factor(nLin)),size=0.6)+
  scale_color_manual(values=pal,name="Number of lineages")+ylab(expression("Power (" * alpha * "=5e-3)"))


plot_grid(plA,plB,rel_widths=c(1,1.3),labels=c("A","B"),ncol=2)

ggsave("~/projects/ILSsims/ILSsims/conceptPlot/groupGenesPow.pdf",height=3,width=9)
