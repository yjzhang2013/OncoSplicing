setwd("/home/u1357/")
library(survival)

mydata<- read.csv("/home/u1357/Cancer_Type.csv")
colnames(mydata)[2:3]<-c("spell","event")
attach(mydata)
time <- spell
event <- event

N<-length(mydata)
result<-list()
HR<-list()

for(i in 4:N)
{
  result[[i]]<-anova(coxph(Surv(time, event)~mydata[,i],data=mydata, ties="breslow"))$Pr[2]
  HR[[i]]<-summary(coxph(Surv(time, event)~mydata[,i],data=mydata, ties="breslow"))$conf.int 
}
data<-t(mydata[1,c(4:N)])
pval<-unlist(result)
gene<-rownames(data)
HRR<-do.call(rbind, lapply(HR, `[`, c(1:4)))
HRR<-data.frame(HRR)
names(HRR) <- c("exp(coef)","exp(-coef)","lower.95","upper.95")
resultsata<-data.frame(gene,pval,HRR)


sigData<-sigData1[sigData1$pval<0.05,]

write.csv(resultsata,file = "cox_cancer.csv",row.names = F)
write.csv(sigData,file = "cox_cancer_0.05p.csv",row.names = F)

