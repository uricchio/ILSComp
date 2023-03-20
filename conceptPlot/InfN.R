library(ggplot2)

# get the corrected sumStat dist
read.table("~/projects/ILSSims/ILSsims/simData/obsInfN.txt")->data

ggplot(data,aes(V2,V1))+geom_abline(slope=1,lwd=0.8)+geom_point(col='red',size=2.1)+geom_abline(slope=1,lwd=0.8)+scale_x_log10()+scale_y_log10()+theme_classic()+xlab(expression("True population size (" * italic(N) * ")"))+ylab(expression("Inferred population size (" * italic(N) * ")"))

ggsave("~/projects/ILSSims/ILSsims/conceptPlot/InfPopSize.pdf",width=6,height=4)
