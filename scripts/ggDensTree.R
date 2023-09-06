library(dplyr)
library(ggtree)

trees<-read.tree("~/projects/ILSSims/ILSsims/obsData/newickTreesForDensiplot.txt")

for (i in seq(1,length(trees))) {
  j = 1
  for (spec in trees[[i]]$tip.label) {
    if(spec == "\'Penicillium-mb\'") {
      trees[[i]]$tip.label[j] <- "\'Penicillium-SP-mb\'"
    } 
    
    if (spec == "\'Penicillium-commune\'") {
      trees[[i]]$tip.label[j] <- "\'Penicillium-bioforme\'"
    }
    
    j <- j+ 1
  }
}

trees.fort <- list(trees[[1]] %>% fortify %>% mutate(tree="gene"))
for (i in seq(1,250)) {
  trees.fort<-c(trees.fort, list(trees[[i]] %>% fortify %>% mutate(tree="gene")));
}

trees.fort<-c(trees.fort, list(trees[[251]] %>% fortify %>% mutate(tree="species")));




plA<-ggdensitree(trees.fort,aes(col=tree,size=as.factor(tree)),alpha=0.5,jitter=0.1)+ geom_tiplab(size=3,alpha=1,col='black') +scale_color_manual(values=c("steelblue","black"))+scale_size_manual(values=c(0.1,1))+ggplot2::xlim(-7,3)
plA<-plA+theme(legend.position="none")

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/DiscTree.pdf",height=3,width=8)
