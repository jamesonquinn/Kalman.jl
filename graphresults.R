library(ggplot2)
library(data.table)

pd = position_dodge(20)

runs = fread("filtertest14.csv",header=F)
names(runs) = c("model","dimensions","repNum",
                "step","numFinkel","numFranken",
                "numPart","nIter","mhAccepts",
                "kl","covdiv","meandiv","entropydiv")

runmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv)),
                by=list(model,dimensions,step,numFinkel,nIter)]

rf = runmeans[model=="finkel",]

table(rf[,list(meankl<10,samp)])
table(rf[,list(meankl<10,histPerLoc)])
table(rf[,list(meankl<10,step)])
table(rf[,list(meankl<10,numFinkel)])
table(runmeans[,list(meankl<10,model)])

(ggplot(data=runmeans[model=="finkel"&meancov<100&step<6,],mapping=
         aes(x=nIter,y=meancov,shape=as.factor(numFinkel),
             color=as.factor(step),group=as.factor(paste0(step,numFinkel)))) +
  geom_line(position=pd) + geom_point(position=pd) + 
    geom_errorbar(aes(ymin=meancov-sdcov,ymax=meancov+sdcov),position=pd,width=40.5))
dim(runs[model=="finkel"&covdiv<100,])





runs1 = fread("filtertest18.csv", fill=TRUE)
runs2 = fread("filtertest19.csv", fill=TRUE)
runs = rbind(runs, runs2, fill=TRUE)
runs[,samp:="log7.5_20"]

#runs = rbind(runs,runs3)
runmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv)),
                by=list(model,dimensions,step,numFinkel,nIter,histPerLoc,samp)]

rf = runmeans[model=="finkel",]



rf = runmeans[model=="finkel",]

table(rf[,list(meankl<10,samp)])
table(rf[,list(meankl<10,histPerLoc)])
table(rf[,list(meankl<10,step)])
table(rf[,list(meankl<10,numFinkel)])
table(runmeans[,list(meankl<10,model)])
