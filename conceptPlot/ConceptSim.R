library(ggplot2)
library(cowplot)

library(png)
library(grid)
library(wesanderson)

img <- readPNG("/Users/uricchio/projects/ILSSims/ILSsims/conceptPlot/UglyTree.png")
g <- rasterGrob(img, interpolate=TRUE)

plA<-qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  geom_blank()+theme_void()

read.table("projects/ILSSims/ILSsims/conceptPlot/conceptData.txt")->dt

mydata<-rbind(data.frame(dist=dt$V1/max(dt$V1),trait=dt$V3,tree="Species Tree"),data.frame(dist=dt$V2/max(dt$V2),trait=dt$V3,tree = "Gene Tree"))

plB<-ggplot(mydata, aes(dist,trait,col=tree,linetype=tree))+geom_smooth(method="lm",alpha=0.1)+geom_point() + theme_classic() +scale_color_manual(values=c("skyblue","darkgray"))
plB<-plB+xlab("Phylogenetic distance")

final<-plot_grid(plA,plB,labels=c("A","B"),rel_widths = c(1.2,1))

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/conceptPlot.pdf",width=11,height=3.6)


