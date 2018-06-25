library(ggplot2)
library(data.table)

runs = fread("filtertest15.csv",header=F)
names(runs) = c("model","dimensions","repNum",
                "step","numFinkel","numFranken",
                "numPart","nIter","mhAccepts",
                "kl","covdiv","meandiv","entropydiv")
runmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv)),
                by=list(model,dimensions,step,numFinkel,nIter)]

pd = position_dodge(2)

(ggplot(data=runmeans[model=="finkel",],mapping=
         aes(x=nIter,y=meancov,color=as.factor(numFinkel),
             shape=as.factor(step),group=as.factor(step))) +
  geom_line(position=pd) + geom_point(position=pd) + 
    geom_errorbar(aes(ymin=meancov-sdcov,ymax=meancov+sdcov),position=pd,width=2))
dim(runs[model=="finkel"&covdiv<100,])
