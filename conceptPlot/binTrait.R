library(ggplot2)
library(cowplot)
library(wesanderson)

# make function to read in files and store pvals and rho values

readFiles<-function(fn,i) {
  data<-read.table(paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.",fn,".",as.character(i),".txt",sep=""))
  return(data.frame(p=data$V3,r=data$V2,n=data$V3,val="null"))
}

getTrueP<-function(fn, i) {
  d<-scan(paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.",fn,".",as.character(i),".txt",sep=""),nlines=1,what="character", skip=1)
  return (as.numeric(d[3:length(d)]))
}

getTrueR<-function(fn, i) {
  d<-scan(paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.",fn,".",as.character(i),".txt",sep=""),nlines=1,what="character", skip=2)
  return (as.numeric(d[3:length(d)]))
}

getTrueNsubs<-function(fn, i) {
  d<-scan(paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.",fn,".",as.character(i),".txt",sep=""),nlines=1,what="character", skip=3)
  return (as.numeric(d[5]))
}

lowDiscPNull = data.frame()
lowDiscPTrue = data.frame()

midDiscPNull = data.frame()
midDiscPTrue = data.frame()

highDiscPNull = data.frame()
highDiscPTrue = data.frame()

for (i in seq(1,200)) {
  f = paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.","lowDiscP",".",as.character(i),".txt",sep="")
  #print(f)
  if(file.size(f) != 0) {
    lowDiscPNull = rbind(lowDiscPNull,readFiles("lowDiscP",i))
    lowDiscPTrue = rbind(lowDiscPTrue,data.frame(p=getTrueP("lowDiscP",i), r=getTrueR("lowDiscP",i),n=getTrueNsubs("lowDiscP",i),val="true"))
  }
}

for (i in seq(1,200)) {
  f = paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.","midDiscP",".",as.character(i),".txt",sep="")
  #print(f)
  if(file.size(f) != 0) {
    midDiscPNull = rbind(midDiscPNull,readFiles("midDiscP",i))
    midDiscPTrue = rbind(midDiscPTrue,data.frame(p=getTrueP("midDiscP",i), r=getTrueR("midDiscP",i),n=getTrueNsubs("midDiscP",i),val="true"))
  }
}


for (i in seq(1,200)) {
  f = paste("~/projects/ILSSims/ILSsims/phyloGWASsims/bin.","highDiscP",".",as.character(i),".txt",sep="")
  #print(f)
  if(file.size(f) != 0) {
    highDiscPNull = rbind(highDiscPNull,readFiles("highDiscP",i))
    highDiscPTrue = rbind(highDiscPTrue,data.frame(p=getTrueP("highDiscP",i), r=getTrueR("highDiscP",i),n=getTrueNsubs("highDiscP",i),val="true"))
  }
}

highDiscP = rbind(highDiscPNull,highDiscPTrue)
midDiscP = rbind(midDiscPNull,midDiscPTrue)
lowDiscP = rbind(lowDiscPNull,lowDiscPTrue)

#

plC<-ggplot(data=lowDiscP,aes(p,color=val,fill=val,alpha=val))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",wes_palette("Darjeeling1")[2]),name="")+
  scale_fill_manual(values=c("lightgray",wes_palette("Darjeeling1")[2]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("p-value")+ggtitle("Low")

plB<-ggplot(data=midDiscP,aes(p,color=val,fill=val,alpha=val))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",wes_palette("Darjeeling1")[3]),name="")+
  scale_fill_manual(values=c("lightgray",wes_palette("Darjeeling1")[3]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("p-value") + theme(legend.position = "none")+ggtitle("Mid")

plA<-ggplot(data=highDiscP,aes(p,color=val,fill=val,alpha=val))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",wes_palette("Darjeeling1")[1]),name="")+
  scale_fill_manual(values=c("lightgray",wes_palette("Darjeeling1")[1]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("p-value") + theme(legend.position = "none")+ggtitle("High")

top_row<-plot_grid(plA,plB,plC, labels=c("A","B","C"),ncol=3)

# QQ plot
makeQQ<-function(df) {
  obs<- sort(df$p)
  expec<-seq(1/(length(df$p)),1,1/(length(df$p )))
  obs<-(-log10(obs))
  expec<-(-log10(expec))
  return(data.frame(o=obs,e=expec))
}

hQQ<-data.frame(qq=makeQQ(highDiscP),discordance="high")
mQQ<-data.frame(qq=makeQQ(midDiscP),discordance="mid")
lQQ<-data.frame(qq=makeQQ(lowDiscP),discordance="low")

QQ<-data.frame(rbind(lQQ,mQQ,hQQ))

plD<-ggplot(data=QQ,aes(qq.e,qq.o,color=discordance))+geom_point()+geom_line() +
  theme_classic() + xlab(expression("expected " * -log[10] * " p-value")) +
  ylab(expression("observed " * -log[10] * " p-value")) +
  geom_abline(slope=1, intercept=0)+scale_color_manual(values=wes_palette("Darjeeling1"))

# ROC

ROC <- function(df,n) {
  df<-df[order(df$p),]
  curve = c(0)
  x = c(0)
  for (i in seq(1, length(df$p))) {
    if (df$val[i] == "true") {
      curve =c(curve,curve[i-1]+1)
    } else {
      curve =c(curve,curve[i-1])
    }
    x = c(x,x[i-1]+1)
  } 
  return(data.frame(c=curve/length(df$p[df$val=="true"]),x=x/length(x),discordance=n))
}

hROC<-ROC(highDiscP,"high")
mROC<-ROC(midDiscP,"mid")
lROC<-ROC(lowDiscP,"low")

ROCdf<-rbind(hROC,mROC,lROC)

plE<-ggplot(ROCdf,aes(x,c,color=discordance))+geom_line(lwd=1.3)+xlab("FPR")+ylab("TPR")+theme_classic()+scale_color_manual(values=wes_palette("Darjeeling1"))
plE

bottom_row<-plot_grid(plD,plE,labels=c("D","E"),ncol=2)

plot_grid(top_row,bottom_row,ncol=1)

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/FigS2.pdf",width=10,height=6.4)




