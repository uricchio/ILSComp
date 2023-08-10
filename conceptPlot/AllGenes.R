library(ggplot2)
library(cowplot)
library(wesanderson)
library(ape)
library(ggtree)

read.table("~/projects/ILSSims/ILSsims/obsData/corr.5719.txt")->data
myPal<-wes_palette("Zissou1", 8, type = "continuous")


myTree<-read.tree("~/projects/ILSSims/ILSsims/obsData/specTree.txt")
sNames<-c("Penicillium-mb", "Penicillium-polonicum", "Penicillium-verrucosum",  "Penicillium-commune", "Penicillium-solitum", "Penicillium-chrysogenum",  "Penicillium-brevicompactum",  "Penicillium-cvjetkovicii" )

#names(myPal)<-sNames
#myTree$edge.length <- rep(1,length(myTree$edge.length))
plTree<-ggtree(myTree,branch.length = 'none')+geom_tiplab(as_ylab=TRUE)


plAll<-ggplot(data=data,aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm',color='black')+theme_classic()+xlab("Phylogenetic distance")+ylab("competition estimate")+scale_color_manual(values = myPal,name='')

plot_grid(plTree,plAll,labels=c("A","B"))

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/AllSpec_SpecTree.pdf",height=2.3*1.2,width=11*1.2)

read.table("~/projects/ILSSims/ILSsims/obsData/corr.5719.txt")->data


myTree<-read.tree("~/projects/ILSSims/ILSsims/obsData/gene5719.tree")
#myTree$edge.length <- rep(1,length(myTree$edge.length))
plTree<-ggtree(myTree,branch.length = 'none')+geom_tiplab(as_ylab=TRUE)

plAll<-ggplot(data=data,aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm',color='black')+theme_classic()+xlab("Phylogenetic distance")+ylab("competition estimate")+scale_color_manual(values = myPal,name='')

plot_grid(plTree,plAll,labels=c("A","B"))

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/AllSpec_5719Tree.pdf",height=2.3*1.2,width=11*1.2)
