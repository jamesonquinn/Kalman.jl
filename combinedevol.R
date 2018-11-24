setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(data.table)
library(gridExtra)
library(ggsci)
library(plyr)
source("theme_publication.R")
pd = position_dodge(20)


runs = fread("lorenz.csv")
if (F) { #temporarily commenting out
      runs = rbind(runs,fread("fixedest10.csv"))
      runs = rbind(runs,fread("fixedest11.csv"))
      runs = rbind(runs,fread("fixedest7.2.csv"))
      #runs = rbind(runs,fread("fixedest8.csv"))
      runs = rbind(runs,fread("fixedest5.2.csv"))
        runs = rbind(runs,fread("fixedest7.2.csv"))
        runs = rbind(runs,fread("fixedest8.csv"))
}
if (F) { #peek at old runs
  runs = fread("fixedest4.csv")
}
#runs = fread("newtest1.csv")
#runs = rbind(runs,fread("newtest2.csv"))
#runs = rbind(runs,fread("newtest3.csv"))
#runs = rbind(runs,fread("newtest4.csv"))
#runs = rbind(runs,fread("newtest5.csv"))
#runs = rbind(runs,fread("newtest6.csv"))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[,Algorithm:=relevel(factor(revalue(model, c("ideal"="Kalman filter",
                                  "truth"="True value",
                                  "observed"="Observed value",
                                  "particle"="Bootstrap PF",
                                  "franken"="Block PF",
                                  "finkel"="Finkelstein PF"))),"Kalman filter")]

runs[,mnfp:=min(nfp,na.rm=T),by=worldid]
#runs=runs[mnfp>350,]
runs[,mhpl:=max(histPerLoc,na.rm=T),by=worldid]
runs[,list(mhpl=max(histPerLoc,na.rm=T)),by=worldid]
#runs = runs[mhpl>40,]
besties = runs[,list(mit=max(nIter,na.rm=T),mf=max(useForward,na.rm=T),md=max(dimension,na.rm=T)),by=worldid]



runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
# tracks over time of everything
wid = 1#27
runs[,step:=steps-1]
#runs[model=="observed",step:=steps]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
runs[model=="truth",mvlvar:=NA]

examplerun = runs[worldid==wid&(!nIter %in% c(0)),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)


(pp1 = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
     )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
     ) + xlab("Time step") + ylab("Sum of loci 3-5") 
     #+ geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2)) 
  + theme_bw() + scale_color_uchicago()
) 

# tracks over time of all 3 algorithms
(pp2 = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
)) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
) + xlab("Time step") + ylab("sum of loci 3:5") 
  #+ geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2)) 
 + theme_bw()+ scale_color_uchicago()
)

grid.arrange(pp1, pp2, ncol=1)
#ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken")&nIter!=30,],aes(y=mvlmean,x=steps,color=paste(model,sampType,histPerLoc,nIter,useForward,mhType))) + geom_line() 
