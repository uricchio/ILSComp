library(ggplot2)
library(wesanderson)

read.table("~/projects/ILSSims/ILSsims/obsData/corr.5719.txt")->data
myPal<-wes_palette("Zissou1", 8, type = "continuous")


plA<-ggplot(data=data[data$V3=='Penicillium-brevicompactum',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[1]),name="")+theme(legend.position = "none")
plB<-ggplot(data=data[data$V3=='Penicillium-chrysogenum',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[2]),name="")+theme(legend.position = "none")
plC<-ggplot(data=data[data$V3=='Penicillium-commune',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[3]),name="")+theme(legend.position = "none")
plD<-ggplot(data=data[data$V3=='Penicillium-cvjetkovicii',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[4]),name="")+theme(legend.position = "none")
plE<-ggplot(data=data[data$V3=='Penicillium-mb',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[5]),name="")+theme(legend.position = "none")
plF<-ggplot(data=data[data$V3=='Penicillium-polonicum',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[6]),name="")+theme(legend.position = "none")
plG<-ggplot(data=data[data$V3=='Penicillium-solitum',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[7]),name="")+theme(legend.position = "none")
plH<-ggplot(data=data[data$V3=='Penicillium-verrucosum',],aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[8]),name="")+theme(legend.position = "none")

plI<-ggplot(data=data,aes(V1,V2,color=V3))+theme_void()+geom_point()+xlim(-0.1,0)+ylim(-0.1,0)+scale_color_manual(values=myPal,name="Species")+theme(legend.position = c(0.5, 0.5))+ theme(legend.text=element_text(size=11))

plot_grid(plA,plB,plC,plD,plE,plF,plG,plH,plI,labels=c("A","B","C","D","E","F","G","H",""),nrow = 3)

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/Fig5.pdf",width=9,height=7.8)

# sp tree

library(ggplot2)
library(wesanderson)

read.table("~/projects/ILSSims/ILSsims/obsData/corr.5719.txt")->data
myPal<-wes_palette("Zissou1", 8, type = "continuous")


plA<-ggplot(data=data[data$V3=='Penicillium-brevicompactum',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[1]),name="")+theme(legend.position = "none")
plB<-ggplot(data=data[data$V3=='Penicillium-chrysogenum',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[2]),name="")+theme(legend.position = "none")
plC<-ggplot(data=data[data$V3=='Penicillium-commune',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[3]),name="")+theme(legend.position = "none")
plD<-ggplot(data=data[data$V3=='Penicillium-cvjetkovicii',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[4]),name="")+theme(legend.position = "none")
plE<-ggplot(data=data[data$V3=='Penicillium-mb',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[5]),name="")+theme(legend.position = "none")
plF<-ggplot(data=data[data$V3=='Penicillium-polonicum',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[6]),name="")+theme(legend.position = "none")
plG<-ggplot(data=data[data$V3=='Penicillium-solitum',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[7]),name="")+theme(legend.position = "none")
plH<-ggplot(data=data[data$V3=='Penicillium-verrucosum',],aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm')+theme_classic()+xlab("Phylogenetic distance")+ylab("RCC")+scale_color_manual(values=c(myPal[8]),name="")+theme(legend.position = "none")

plI<-ggplot(data=data,aes(V4,V2,color=V3))+theme_void()+geom_point()+xlim(-0.1,0)+ylim(-0.1,0)+scale_color_manual(values=myPal,name="Species")+theme(legend.position = c(0.5, 0.5))+ theme(legend.text=element_text(size=11))

plot_grid(plA,plB,plC,plD,plE,plF,plG,plH,plI,labels=c("A","B","C","D","E","F","G","H",""),nrow = 3)


ggsave("~/projects/ILSSims/ILSsims/conceptPlot/realData.SPtree.pdf",width=9,height=7.8)

cor.test(data$V2,data$V4)
cor.test(data[data$V3=='Penicillium-verrucosum',]$V4,data[data$V3=='Penicillium-verrucosum',]$V2)

