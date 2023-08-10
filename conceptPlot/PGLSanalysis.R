library(ape)
library(geiger)
library(nlme)
library(phytools)

# function to read in PGLS data
readAndDoPGLS<-function(n,i) {
  
    f = paste('~/projects/ILSSims/ILSsims/PGLSsims/pgls.',n, 'DiscP.',as.character(i),'.txt',sep='')
    #print(f)
    read.table(f,row.names=2)->vals # nrows set so that we only get one focal species, cvjetkovicii
    vals<-cbind(vals, data.frame(V2=row.names(vals)))
    
    read.tree("~/projects/ILSSims/ILSsims/simData/spec.txt")->spec
    
    nullmods = c()
    truemods = c()
    for(j in seq(2,2)) { # switch to seq(0,7) for all species -- 2 is P. commune
      
        locVals <- vals[vals$V3 == j,]
        
        # scale phenotype vals
        locVals$V5 <- scale(locVals$V5)
        
        newspec <- spec
        for (i in seq(length(newspec$tip.label))) {
            newspec$tip.label[i]= paste(newspec$tip.label[i],'-',as.character(j),sep='')
        }
    
        obj<-name.check(newspec,locVals)
        dropSpec<-drop.tip(newspec,obj$tree_not_data)
        
        #print(dropSpec$tip.label)
        
        # get 'true' p-value
        bm <- corBrownian(1, dropSpec, form = ~ V2)
        mod<-try(gls(V4 ~ V5, data=locVals, correlation=bm))
        myCoef <- try(coef(summary(mod))[,4][2])
        if (is.numeric(myCoef)) {
          truemods <- c(truemods,myCoef)
         }
        
        
        # brownian term
        
        # get null dis of p-vals
        for (c in seq(6,length(colnames(locVals))-1)) {
          regVals =data.frame(tr=locVals$V4,gd=locVals[,c],sp=locVals$V2)
          row.names(regVals)<-regVals$sp
          obj<-name.check(newspec,regVals)
          dropSpec<-drop.tip(newspec,obj$tree_not_data)
          
          bm <- corBrownian(1, dropSpec, form = ~sp)

          mod<-try(gls( tr ~ gd, data=regVals, correlation=bm))
          myCoef <- try(coef(summary(mod))[,4][2])
          if (is.numeric(myCoef)) {
              nullmods <- c(nullmods,myCoef)
          }
        }
        
    }
    tm <- data.frame(p=truemods,val="true")
    nm <- data.frame(p=nullmods,val="null")
    
    ret<-rbind(tm,nm)
    
    return(ret)
}

highDiscP <- data.frame()
midDiscP <- data.frame()
lowDiscP <- data.frame()

for (i in seq(1,100)) {
  print(i)
  #try(res<-readAndDoPGLS("high",i))
  highDiscP<-rbind(highDiscP,readAndDoPGLS("high",i) )
  lowDiscP<-rbind(lowDiscP,readAndDoPGLS("low",i) )
  midDiscP<-rbind(midDiscP,readAndDoPGLS("mid",i) )
}


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

plD<-ggplot(data=QQ,aes(qq.e,qq.o,color=discordance))+geom_line(lwd=1.4) +
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

bottom_row<-plot_grid(plD,plE,labels=c("D","E"),ncol=2,rel_heights = c(0.7,1))

plot_grid(top_row,bottom_row,ncol=1)

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/FigS1.pdf",width=10,height=5)




