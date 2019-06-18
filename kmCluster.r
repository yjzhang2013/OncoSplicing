setwd("/home/u1357/")

library(ConsensusClusterPlus)
data <- read.csv("/home/u1357/kmCluster.csv", header=TRUE,row.names=1)
title = "/home/u1357/spc/"

d <- data.matrix(data)

results = ConsensusClusterPlus(d,maxK=15,reps=100,pItem=0.8,pFeature=1,
    title=title,clusterAlg="km",distance="euclidean",seed=1262118388,plot="png")

x3<-results[[3]][["consensusClass"]]
x4<-results[[4]][["consensusClass"]]
x5<-results[[5]][["consensusClass"]]
x6<-results[[6]][["consensusClass"]]
x7<-results[[7]][["consensusClass"]]
x8<-results[[8]][["consensusClass"]]
x9<-results[[9]][["consensusClass"]]
x10<-results[[10]][["consensusClass"]]
x11<-results[[11]][["consensusClass"]]
x12<-results[[12]][["consensusClass"]]
x13<-results[[13]][["consensusClass"]]
x14<-results[[14]][["consensusClass"]]
x15<-results[[15]][["consensusClass"]]

z<-cbind(x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15)
write.csv(z,file="/home/u1357/resultKMcluster.csv")





