
library(wesanderson)
library(ggplot2)
library(cowplot)

# get the corrected sumStat dist
read.table("~/projects/ILSSims/ILSsims/obsData/corrPvalDist.txt")->data
# get the obs sumStat dist
read.table("~/projects/ILSSims/ILSsims/obsData/realPvalDist.txt")->dataR

data$V1<-replace(data$V1, data$V1 =="Penicillium-mb", "Penicillium-SP-mb")
data$V1<-replace(data$V1, data$V1 =="Penicillium-commune", "Penicillium-bioforme")

dataR$V3<-replace(dataR$V3, dataR$V3 =="Penicillium-mb", "Penicillium-SP-mb")
dataR$V3<-replace(dataR$V3, dataR$V3 =="Penicillium-commune", "Penicillium-bioforme")


myData<-rbind(data.frame(r=dataR$V1,species=dataR$V3,type="Observed"),data.frame(r=data$V2,species=data$V1,type="Null"))


# plot of each species with its null dist and obs dist

myPal<-wes_palette("Zissou1", 8, type = "continuous")

plA<-ggplot(data=myData[myData$species=='Penicillium-brevicompactum',],aes(r,color=type,fill=type,alpha=type))+geom_density(bw=0.02)+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[1]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[1]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-brevicompactum")

plB<-ggplot(data=myData[myData$species=='Penicillium-chrysogenum',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[2]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[2]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-chrysogenum")

plC<-ggplot(data=myData[myData$species=='Penicillium-bioforme',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[3]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[3]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-bioforme")

plD<-ggplot(data=myData[myData$species=='Penicillium-cvjetkovicii',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[4]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[4]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-cvjetkovicii")

plE<-ggplot(data=myData[myData$species=='Penicillium-SP-mb',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[5]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[5]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-SP-mb")

plF<-ggplot(data=myData[myData$species=='Penicillium-polonicum',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[6]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[6]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-polonicum")

plG<-ggplot(data=myData[myData$species=='Penicillium-solitum',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[7]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[7]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-solitum")

plH<-ggplot(data=myData[myData$species=='Penicillium-verrucosum',],aes(r,color=type,fill=type,alpha=type))+geom_density()+theme_classic()+ylab("probability density")+
  scale_color_manual(values=c("lightgray",myPal[8]),name="")+
  scale_fill_manual(values=c("lightgray",myPal[8]),name="") + 
  scale_alpha_manual(values=c(1,0.4),guide="none")+ xlab("Correlation coefficient") +ggtitle("Penicillium-verrucosum")

plot_grid(plA,plB,plC,plD,plE,plF,plG,plH,labels=c("A","B","C","D","E","F","G","H"),nrow = 3)

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/Fig3.pdf",height=7.5,width=14*(7.5/9))


# qqplots
read.table("~/projects/ILSSims/ILSsims/obsData/QQplotData.txt")->QQ

QQ$V1<-replace(QQ$V1, QQ$V1 =="Penicillium-mb", "Penicillium-SP-mb")
QQ$V1<-replace(QQ$V1, QQ$V1 =="Penicillium-commune", "Penicillium-bioforme")

# uncorrected
plA<-ggplot(data=QQ,aes(-log10(V5),-log10(V6),color=V1))+geom_point()+geom_line() +
  theme_classic() + xlab(expression("expected " * -log[10] * " p-value")) +
  ylab(expression("observed " * -log[10] * " p-value")) +
  geom_abline(slope=1, intercept=0)+scale_color_manual(values=myPal,name="Species",guide="none")

# corrected with MSC
plB<-ggplot(data=QQ,aes(-log10(V5),-log10(V4),color=V1))+geom_point()+geom_line() +
  theme_classic() + xlab(expression("expected " * -log[10] * " p-value")) +
  ylab(expression("observed " * -log[10] * " p-value")) +
  geom_abline(slope=1, intercept=0)+scale_color_manual(values=myPal,name="Species")

plot_grid(plA,plB,labels=c("A","B"),rel_widths = c(1,1.7))

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/Fig4.pdf",height=3,width=9)
