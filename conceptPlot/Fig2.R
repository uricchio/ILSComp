library(ape)
library(geiger)
library(nlme)
library(phytools)

read.table("~/projects/ILSSims/ILSsims/simData/lowDiscP.txt")->vals
read.tree("~/projects/ILSSims/ILSsims/simData/spec.txt")->spec

pv<-c()
tPV<-c()

for (mySpec in spec$tip.label) { 

  valsTemp<-vals[vals$V6==0 & vals$V4 == mySpec,]
  modTrue<-gls(V3  ~ V2, correlation = corBrownian(phy = spec, form = ~V5),data =valsTemp, method = "ML",na.action=na.omit)



    for (i in seq(1,max(vals$V6))) {
      valsTemp<-vals[vals$V6==i & vals$V4 == mySpec,]
      mod<-gls(V3  ~ V2, correlation = corBrownian(phy = spec, form = ~V5),data =valsTemp, method = "ML",na.action=na.omit)
      pv<-c(pv,coef(summary(mod))[,4][2])
    }
 
  tPV<-c(tPV,coef(summary(modTrue))[,4][2])
  pv<-c(t(t(pv)))

}



