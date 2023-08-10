library(ggplot2)
library(cowplot)
library(ggtree)

read.table("~/projects/ILSSims/ILSsims/obsData/pValsNull/sumDataDiffs.txt")->data

myTree<-read.tree("~/projects/ILSSims/ILSsims/obsData/coalUnitTreePlot.txt")

plA<-ggplot(data,aes(V1,V2))+geom_point()+geom_smooth(se=FALSE)+theme_classic()+xlab(expression(italic(N[e])))+ylab(expression(italic(D)))
plA<-plA+geom_point(data=data.frame(V1=c(0.0699),V2=c(0.0242)),color='red')+geom_hline(yintercept=0.0242, linetype="dashed")+scale_x_continuous(limits=c(0.01,0.13),breaks=c(0.025,0.05,0.075,0.1,0.125))

plB<-ggtree(myTree)+ theme_tree2()+geom_tiplab(size=8,as_ylab=TRUE) 
plot_grid(plA,plB,labels=c("A","B"),rel_widths = c(1,1.5))

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/InferPopSize.pdf",width=9,height=2.5)
