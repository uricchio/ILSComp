library(ggplot2)
library(cowplot)

library(png)
library(grid)
library(wesanderson)

mydata<-rbind(data.frame(dist=dt$V1/max(dt$V1),trait=dt$V3,tree="Species Tree"),data.frame(dist=dt$V2/max(dt$V2),trait=dt$V3,tree = "Gene Tree"))

plB<-ggplot(mydata, aes(dist,trait,col=tree,linetype=tree))+geom_smooth(method="lm",alpha=0.1,se=FALSE)+geom_point() + theme_classic() +scale_color_manual(name="",values=c("skyblue","darkgray"))
plB<-plB+xlab("Phylogenetic distance")+scale_linetype_manual(values=c(1,2),guide="none")+theme(text = element_text(size=14))

#final<-plot_grid(plA,plB,labels=c("A","B"),rel_widths = c(1.2,1))

plB

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/conceptPlotPD.pdf",width=8,height=3.9)


