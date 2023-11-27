library(ggplot2)
library(cowplot)
library(wesanderson)
library(ape)
library(ggtree)


# this script can be used to recreate Figure 5 of Louw et al 2023, https://www.biorxiv.org/content/10.1101/2023.05.11.540388v3
# Figure 5 shows a cladogram of the phylogenetic tree structure for the locus that demonstrated the strongest phylogenetic signal in our dataseet, and
# plots of the correlations between trait values and pairwise distance for each focal species of Penecillium in our dataset

# first we read in the trait data and replace species names with current nomenclature for two species
read.table("~/projects/ILSSims/ILSsims/obsData/corr.5719.txt")->data

data$V3<-replace(data$V3, data$V3 =="Penicillium-mb", "Penicillium-SP-mb")
data$V3<-replace(data$V3, data$V3 =="Penicillium-commune", "Penicillium-bioforme")

# same for phylogenetic data related to the tree
myTree<-read.tree("~/projects/ILSSims/ILSsims/obsData/gene5719.tree")
myTree$tip.label[8]<-"Penicillium-SP-mb"
myTree$tip.label[5]<-"Penicillium-bioforme"


# now we use ggtree to plot a cladogram of the tree
plTree<-ggtree(myTree,branch.length = 'none')+geom_tiplab(as_ylab=TRUE)


# now we plot the correlations
plAll<-ggplot(data=data,aes(V1,V2,color=V3))+geom_point()+geom_smooth(method='lm',color='black',se=FALSE)+theme_classic()+xlab("Phylogenetic distance")+ylab("competition estimate (RCC)")+scale_color_manual(values = myPal,name='')
plAll<-plAll+geom_smooth(data=data,aes(V1,V2,color=V3), method='lm',se=FALSE,fullrange=TRUE,lty=2,aplha=0.5,size=0.4)
plot_grid(plTree,plAll,labels=c("A","B"),rel_widths=c(1,1.3))

# gene tree stats for table
cor.test(data$V2[data$V3=="Penicillium-SP-mb"],data$V1[data$V3=="Penicillium-SP-mb"])
cor.test(data$V2[data$V3=="Penicillium-bioforme"],data$V1[data$V3=="Penicillium-bioforme"])
cor.test(data$V2[data$V3=="Penicillium-polonicum"],data$V1[data$V3=="Penicillium-polonicum"])
cor.test(data$V2[data$V3=="Penicillium-verrucosum"],data$V1[data$V3=="Penicillium-verrucosum"])
cor.test(data$V2[data$V3=="Penicillium-solitum"],data$V1[data$V3=="Penicillium-solitum"])
cor.test(data$V2[data$V3=="Penicillium-chrysogenum"],data$V1[data$V3=="Penicillium-chrysogenum"])
cor.test(data$V2[data$V3=="Penicillium-brevicompactum"],data$V1[data$V3=="Penicillium-brevicompactum"])
cor.test(data$V2[data$V3=="Penicillium-cvjetkovicii"],data$V1[data$V3=="Penicillium-cvjetkovicii"])

ggsave("~/projects/ILSsims/ILSsims/conceptPlot/Fig5Updated.pdf",width=12,height=3)

# the following code performs the same analysis as above, but now using the species tree instead of the locus of interest

myTree<-read.tree("~/projects/ILSSims/ILSsims/simData/spec.txt")
myTree$tip.label[8]<-"Penicillium-SP-mb"
myTree$tip.label[5]<-"Penicillium-bioforme"

plTree<-ggtree(myTree,branch.length = 'none')+geom_tiplab(as_ylab=TRUE)

plAll<-ggplot(data=data,aes(V4,V2,color=V3))+geom_point()+geom_smooth(method='lm',color='black',se=FALSE)+theme_classic()+xlab("Phylogenetic distance")+ylab("competition estimate (RCC)")+scale_color_manual(values = myPal,name='')
plAll<-plAll+geom_smooth(data=data,aes(V4,V2,color=V3), method='lm',se=FALSE,fullrange=TRUE,lty=2,aplha=0.5,size=0.4)+ylim(c(1,1.6))
plot_grid(plTree,plAll,labels=c("A","B"),rel_widths=c(1,1.3))

# species tree stats
cor.test(data$V2[data$V3=="Penicillium-SP-mb"],data$V4[data$V3=="Penicillium-SP-mb"])
cor.test(data$V2[data$V3=="Penicillium-bioforme"],data$V4[data$V3=="Penicillium-bioforme"])
cor.test(data$V2[data$V3=="Penicillium-polonicum"],data$V4[data$V3=="Penicillium-polonicum"])
cor.test(data$V2[data$V3=="Penicillium-verrucosum"],data$V4[data$V3=="Penicillium-verrucosum"])
cor.test(data$V2[data$V3=="Penicillium-solitum"],data$V4[data$V3=="Penicillium-solitum"])
cor.test(data$V2[data$V3=="Penicillium-chrysogenum"],data$V4[data$V3=="Penicillium-chrysogenum"])
cor.test(data$V2[data$V3=="Penicillium-brevicompactum"],data$V4[data$V3=="Penicillium-brevicompactum"])
cor.test(data$V2[data$V3=="Penicillium-cvjetkovicii"],data$V4[data$V3=="Penicillium-cvjetkovicii"])

ggsave("~/projects/ILSsims/ILSsims/conceptPlot/FigS4Updated.pdf",width=12,height=3)

# get the corrected sumStat dist
#read.table("~/projects/ILSSims/ILSsims/obsData/corrPvalDist.txt")->data
# get the obs sumStat dist
#read.table("~/projects/ILSSims/ILSsims/obsData/realPvalDist.txt")->dataR

#p_mb<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-mb",]$V1,data,"Penicillium-mb")
#p_ch<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-chrysogenum",]$V1,data,"Penicillium-chrysogenum")
#p_co<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-commune",]$V1,data,"Penicillium-commune")
#p_cv<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-cvjetkovicii",]$V1,data,"Penicillium-cvjetkovicii")
#p_po<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-polonicum",]$V1,data,"Penicillium-polonicum")
#p_ve<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-verrucosum",]$V1,data,"Penicillium-verrucosum")
#p_so<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-solitum",]$V1,data,"Penicillium-solitum")
#p_br<-get_p(dataR[dataR$V4=="5719" & dataR$V3=="Penicillium-brevicompactum",]$V1,data,"Penicillium-brevicompactum")

