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





runs1 = fread("filtertest21.csv", fill=TRUE)
runs2 = fread("filtertest22.csv", fill=TRUE)
runs3 = fread("filtertest23.csv", fill=TRUE)
runs4 = fread("filtertest24.csv", fill=TRUE)
runs = rbind(runs1, runs2,runs3,runs4, fill=TRUE)
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]

#runs = rbind(runs,runs3)
runmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),meantime=mean(runtime)),
                by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType)]
runmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),meantime=mean(runtime)),
                by=list(model,particles,nIter,histPerLoc,sampType)] #no steps

rf = runmeans[model=="finkel",]
rnf = runmeans2[model!="finkel",]


table(rf[,list(meankl<10,sampType)])
table(rf[,list(meankl<10,histPerLoc)])
table(rf[,list(meankl<10,steps)])
table(rf[,list(meankl<10,nfp)])
table(runmeans[,list(meankl<10,model)])
