runs = fread("fixedest4.csv")
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 4
runs[,step:=steps-1]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                    by=list(model,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward)] #no steps
rf = runmeans[model=="finkel"|model=="ideal",]
rnf = runmeans2[model!="finkel",]
ggplot(data=runmeans[steps>4,],aes(x=sqrt(meankl),y=sqrt(meansq),shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],not),aes(x=meankl,y=meancov,shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
nice(d=arunmeans[model %in% c("finkel","franken"),],not)
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],not),aes(x=meankl,y=meancov,shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
ggplot(data=nice(d=arunmeans,not),aes(x=sqrt(meankl),y=sqrt(meansq),shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
ggplot(data=runmeans[steps>4,],aes(x=sqrt(meankl),y=sqrt(meansq),shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=runmeans[steps>4&model %in% c("finkel","franken"),],aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=nice(d=runmeans[steps>4&model %in% c("finkel","franken"),],not),aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],it),aes(x=sqrt(meankl),y=sqrt(meansq),shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
arunmeans[model %in% c("finkel","franken"),],it)
nice(d=arunmeans[model %in% c("finkel","franken"),],it)
nice(d=allrunmeans2[model %in% c("finkel","franken"),],it)
nice(d=allrunmeans2[model %in% c("finkel","franken"),])
ggplot(data=nice(d=runmeans[steps>4&model %in% c("finkel","franken"),],not),aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=nice(d=arunmeans[steps>4&model %in% c("finkel","franken"),],it),aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
nice(d=arunmeans[steps>4&model %in% c("finkel","franken"),],it)
ggplot(data=arunmeans[steps>4&model %in% c("finkel","franken")&particles>40,],aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
runs = fread("fixedest5.csv")
runs = fread("fixedest5.csv")
runs = fread("fixedest5.2.csv")
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
?fread
runs = rbind(runs,fread("fixedest8.csv",skip=1))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 4
head(runs[dimension==60,]
)
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 19
runs[,step:=steps-1]
runs[step==0,mvlmean:=0]
runs[model=="observed"&step==0,mvlmean:=NA]
runs[model!="finkel",nIter:=NA]
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "observed", "truth")
                  &nIter %in% c(NA,80),],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
unique(runs[dimension==60,worldid])
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken", "particle")
                  &nIter %in% c(NA,80),],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken", "particle")
                  &nIter %in% c(NA,40),],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
unique(runs[dimension==90,worldid])
wid = 27
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0,mvlmean:=0]
runs[model=="observed"&step==0,mvlmean:=NA]
runs[model!="finkel",nIter:=NA]
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "observed", "truth")
                  &nIter %in% c(NA,80),],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken", "particle")
                  &nIter %in% c(NA,40),],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken", "particle")
                  ,],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken", "particle")
                  &nIter %in% c(NA,120),],aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
                  )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
                  ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[model!="finkel",nIter:=NA]
runs[model %in% c("ideal")&step==0,mvlvar:=NA]
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[model %in% c("ideal")&step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs[step==0,mvlvar:=NA]
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
examplerun = rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs[model=="truth",mvlvar:=NA]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(labels = c("Kalman", "Finkelstein PF", "block PF", "bootstrap PF"),
                       values = c("ideal", "finkel", "franken", "particle"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(labels = c("Kalman", "Finkelstein PF", "block PF", "bootstrap PF"),
                       breaks = c("ideal", "finkel", "franken", "particle"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(labels = c("Kalman", "Finkelstein PF", "block PF", "bootstrap PF"),
                       breaks = c("ideal", "finkel", "franken", "particle"),
                       values = c("#000","#E33","#4F8","#411"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=model,group=model,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(labels = c("Kalman", "Finkelstein PF", "block PF", "bootstrap PF"),
                       breaks = c("ideal", "finkel", "franken", "particle"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
levels(runs[,model])
unique(runs[,model])
runs[,Algorithm:=revalue(model, c("truth"="True value",
                                  "observed"="Observed value",
                                  "ideal"="Kalman filter",
                                  "particle"="Bootstrap PF",
                                  "franken"="Block PF",
                                  "finkel"="Finkelstein PF"))]
library(plyr)
runs[,Algorithm:=revalue(model, c("truth"="True value",
                                  "observed"="Observed value",
                                  "ideal"="Kalman filter",
                                  "particle"="Bootstrap PF",
                                  "franken"="Block PF",
                                  "finkel"="Finkelstein PF"))]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman Filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=model,
             )) + geom_line() + ggtitle("Comparison of"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs[,Algorithm:=revalue(model, c("truth"="True value",
                                  "observed"="Observed value",
                                  "ideal"="Kalman filter",
                                  "particle"="Bootstrap PF",
                                  "franken"="Block PF",
                                  "finkel"="Finkelstein PF"))]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
unique(examplerun[,Algorithm])
head(runs)
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[,Algorithm:=revalue(model, c("truth"="True value",
                                  "observed"="Observed value",
                                  "ideal"="Kalman filter",
                                  "particle"="Bootstrap PF",
                                  "franken"="Block PF",
                                  "finkel"="Finkelstein PF"))]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
runs[model=="truth",mvlvar:=NA]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(model="finkel",step=0,mvlmean=0),fill=T)
unique(examplerun[,Algorithm])
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Evolution of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=model,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman Filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=model,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + scale_linetype_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                          values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("#000000","#E03030","#40F080","#401010"))
  + scale_linetype_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF")
  )
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("Kalman filter"="#000000","Finkelstein PF"="#E03030","Block PF"="#40F080","Bootstrap PF"="#401010"))
  + scale_linetype_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                          values = c("#000000","#E03030","#40F080","#401010"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + scale_color_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                       values = c("Kalman filter"="#000000","Finkelstein PF"="#E03030","Block PF"="#40F080","Bootstrap PF"="#401010"))
  + scale_linetype_manual(breaks = c("Kalman filter", "Finkelstein PF", "Block PF", "Bootstrap PF"),
                          values = c("Kalman filter"="solid","Finkelstein PF"="#dotted","Block PF"="#dashed","Bootstrap PF"="#dotdash"))
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Comparison of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series of filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("Sum of loci 3-5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs[,Algorithm:=relevel(revalue(model, c("truth"="True value",
                                          "observed"="Observed value",
                                          "ideal"="Kalman filter",
                                          "particle"="Bootstrap PF",
                                          "franken"="Block PF",
                                          "finkel"="Finkelstein PF")), "Kalman filter")]
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
runs[,i:=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]
runs[,Algorithm:=relevel(revalue(model, c("truth"="True value",
                                          "observed"="Observed value",
                                          "ideal"="Kalman filter",
                                          "particle"="Bootstrap PF",
                                          "franken"="Block PF",
                                          "finkel"="Finkelstein PF")), "Kalman filter")]
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
runs[model=="truth",mvlvar:=NA]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("Sum of loci 3-5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
runs[,Algorithm]
runs[,Algorithm:=relevel(revalue(model, c("truth"="True value",
                                          "observed"="Observed value",
                                          "ideal"="Kalman filter",
                                          "particle"="Bootstrap PF",
                                          "franken"="Block PF",
                                          "finkel"="Finkelstein PF")), "Kalman filter")]
runs[,Algorithm:=revalue(model, c("ideal"="Kalman filter",
                                  "truth"="True value",
                                  "observed"="Observed value",
                                  "particle"="Bootstrap PF",
                                  "franken"="Block PF",
                                  "finkel"="Finkelstein PF"))]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("Sum of loci 3-5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
?revalue
runs[,Algorithm:=relevel(factor(revalue(model, c("ideal"="Kalman filter",
                                                 "truth"="True value",
                                                 "observed"="Observed value",
                                                 "particle"="Bootstrap PF",
                                                 "franken"="Block PF",
                                                 "finkel"="Finkelstein PF"))),"Kalman filter")]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("Sum of loci 3-5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
             ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
) + theme_bw()
(pp = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
             aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
             )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
             ) + xlab("Time step") + ylab("Sum of loci 3-5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
  + theme_bw()
)
(pp2 = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
              aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
              )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
              ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
  + theme_bw()
)
grid.arrange(pp1, pp2, ncol=1)
(pp1 = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
              aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
              )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
              ) + xlab("Time step") + ylab("Sum of loci 3-5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
  + theme_bw()
)
(pp2 = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
              aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
              )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
              ) + xlab("Time step") + ylab("sum of loci 3:5")
  + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
  + theme_bw()
)
grid.arrange(pp1, pp2, ncol=1)
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                    by=list(model,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward)] #no steps
rf = runmeans[model=="finkel"|model=="ideal",]
rnf = runmeans2[model!="finkel",]
nicepar = function(d=allrunmeans2) {d[sampType %in% c("none","bkf.SampleLog") & d[,mhType] %in% c("none","bkf.MhCompromise"),]}
nicehpl = function(d=allrunmeans2) {d[histPerLoc %in% c(0,15),]}
niceit = function(d=allrunmeans2) {d[nIter %in% c(NA,40),]}
nicen2 = function(d=allrunmeans2) {d[particles %in% c(200,40000,8000),]}
nicen = function(d=allrunmeans2) {d[particles %in% c(400,160000,32000),]}
nicefw = function(d=allrunmeans2) {d[useForward %in% c(NA,1.),]}
nicedim = function(d=allrunmeans2) {d[dimension==30,]}
niceval = function(d=allrunmeans2) {d[meankl<50&meansq<50,]}
nice = function(...,d=allrunmeans2) {
  myargs = sapply(substitute(list(...)),as.character) #includes a spurious "list" but who cares.
  ma = myargs
  if ("not" %in% myargs) {
    myargs = setdiff(c("par", "hpl", "it", "n", "fw", "dim", "val"),myargs)
  }
  if ("par" %in% myargs) {
    d = nicepar(d)
  }
  if ("hpl" %in% myargs) {
    d = nicehpl(d)
  }
  if ("it" %in% myargs) {
    d = niceit(d)
  }
  if ("n" %in% myargs) {
    if ("alt" %in% ma) {
      d = nicen2(d)
    } else {
      d = nicen(d)
    }
  }
  if ("fw" %in% myargs) {
    d = nicefw(d)
  }
  if ("dim" %in% myargs) {
    d = nicedim(d)
  }
  if ("val" %in% myargs) {
    d = niceval(d)
  }
  d
}
ggplot(data=arunmeans[steps>4&model %in% c("finkel","franken")&particles>40&nIter!=0,],aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=arunmeans[steps>4&model %in% c("finkel","franken")&particles>40&!(nIter %in% 0,],aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],it),aes(x=sqrt(meankl),y=sqrt(meansq),shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
ggplot(data=runmeans[model=="finkel"&particles==400&nIter>10,],aes(x=sqrt(meankl),y=sqrt(meansq),group=paste(params,particles,nIter,histPerLoc))
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params,shape=factor(histPerLoc))
)
ggplot(data=runmeans[model=="ideal"|(particles==400&nIter>10),],aes(x=sqrt(meankl),y=sqrt(meansq),group=paste(model,params,particles,nIter,histPerLoc,dimension))
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params,shape=factor(histPerLoc))
)
ggplot(data=runmeans[model=="ideal"|(particles>300&nIter>10&params=="bkf.SampleLog bkf.MhCompromise"),],aes(x=sqrt(meankl),y=sqrt(meansq),group=paste(model,params,particles,nIter,histPerLoc,dimension))
) + geom_path(aes(linetype=factor(nIter),color=factor(nIter))
) + geom_point(aes(shape=factor(histPerLoc))
)
#+ scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,8))
ggplot(data=runmeans[model!="finkel"|(nIter==70&particles==800),],aes(x=sqrt(meankl),y=sqrt(meansq),shape=paste(mhType,sampType),group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=runmeans[model=="franken",],aes(x=sqrt(meankl),y=sqrt(meansq),shape=model,group=model)
) + geom_path(aes(linetype=model,color=factor(particles))
) + geom_point(
)
+ scale_x_continuous(limits=c(0,20)) + scale_y_continuous(limits=c(1,10))
qplot(sqrt(meankl),meansq,data=runmeans[,],color=model,main="zoom 30 colors")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans2[,],shape=paste(model,nIter),color=log(particles),main="finkel by particles")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(particles),main="finkel by particles")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(particles),main="finkel by particles zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(steps),main="finkel by steps zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(nIter),main="finkel by nIter zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(sampType),main="finkel by sampType zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(mhType),main="finkel by mhType zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(histPerLoc),main="finkel by histPerLoc zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(paste(sampType,mhType)),main="finkel by samp/mh zoom")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[,],shape=model,color=as.factor(paste(sampType,mhType)),main="finkel by samp/mh zoom")
qplot(sqrt(meankl),sqrt(meansq),data=rf[sampType=="bkf.SampleLog"&mhType=="bkf.MhCompromise",],shape=as.factor(nIter),color=as.factor(particles),size=steps,main="log/compromise by nIter/particles")
qplot(sqrt(meankl),sqrt(meansq),data=runmeans[(model!="finkel"|(sampType=="bkf.SampleLog"&mhType=="bkf.MhCompromise")),],shape=model,color=as.factor(pmin(particles,1001)),main="nonfinkel plus log/compromise by nIter/particles")
ggplot(data=arunmeans[steps>4&model %in% c("finkel","franken")&particles>40&!(nIter %in% 0),],aes(x=meankl,y=meancov,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
)
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],it),aes(x=sqrt(meankl),y=sqrt(meansq),shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
runs[,params:=paste(sampType,mhType)]
runs[,params:=paste(sampType,mhType)]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                    by=list(model,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],it),aes(x=sqrt(meankl),y=sqrt(meansq),shape=params,group=params)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
rf = runmeans[model=="finkel"|model=="ideal",]
rnf = runmeans2[model!="finkel",]
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],it),aes(x=sqrt(meankl),y=sqrt(meansq),shape=model,group=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
ggplot(data=nice(d=arunmeans[model %in% c("finkel","franken"),],it),aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
somnice = nice(d=arunmeans[model %in% c("finkel","franken"),],it)
somnice = nice(d=arunmeans[model %in% c("finkel","franken"),],it)
unique(somenice[,particles])
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it)
unique(somenice[,particles])
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n)
unique(somenice[,particles])
ggplot(data=nice(d=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
qplot(meankl,meansq,data=nice(val),shape=model,color=model,
      main="all from run 30")
qplot(meankl,meansq,data=nice(not,n),shape=model,color=as.factor(particles),
      main="finkel by particles")
#qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(steps),main="finkel by steps zoom")
#finkel by steps zoom
qplot(meankl,meansq,data=nice(not,it),shape=model,color=paste(nIter),
      main="finkel by nIter")
#qplot(meankl,meansq,data=allrunmeans2[meankl<5&meankl>0&meansq<50,],shape=paste(useForward),color=paste(nIter),main=
#        "by nIter useForward ")
#qplot(meankl,meansq,data=runmeans[model=="finkel"&meankl<5&meankl>0&meansq<50,],shape=paste(useForward),color=paste(nIter),main=
#       "finkel by nIter useForward ")
qplot(meankl,meancov,data=nice(not,it,fw),shape=paste(useForward),color=paste(nIter),main=
        "meancov by nIter useForward ")
qplot(meandiff,meankl-meancov-meandiff,data=nice(not,it,fw),shape=paste(useForward),color=paste(nIter),main=
        "meancov by nIter useForward ")
qplot(meankl,meansq,data=nice(n),shape=paste(useForward),color=paste(nIter),main=
        "sqerr by nIter useForward ")
qplot(meankl,meansq,data=nice(alt,n),shape=paste(useForward),color=paste(nIter),main=
        "sqerr by nIter useForward fewer")
qplot(meandiff,meancov,data=nice(),shape=paste(useForward),color=paste(nIter),main=
        "finkel by nIter useForward ")
ggplot(data=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
somenice[steps==2,]
sn = somenice[steps==2,]
View(sn)
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl)
unique(somenice[,particles])
dim(sn)
View(sn)
View(sn)
dim(nicehpl(sn))
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl)
nice = function(...,d=allrunmeans2) {
  myargs = sapply(substitute(list(...)),as.character) #includes a spurious "list" but who cares.
  ma = myargs
  print(myargs)
  if ("not" %in% myargs) {
    myargs = setdiff(c("par", "hpl", "it", "n", "fw", "dim", "val"),myargs)
  }
  if ("par" %in% myargs) {
    d = nicepar(d)
  }
  if ("hpl" %in% myargs) {
    d = nicehpl(d)
  }
  if ("it" %in% myargs) {
    d = niceit(d)
  }
  if ("n" %in% myargs) {
    if ("alt" %in% ma) {
      d = nicen2(d)
    } else {
      d = nicen(d)
    }
  }
  if ("fw" %in% myargs) {
    d = nicefw(d)
  }
  if ("dim" %in% myargs) {
    d = nicedim(d)
  }
  if ("val" %in% myargs) {
    d = niceval(d)
  }
  d
}
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl)
nice = function(...,d=allrunmeans2) {
  myargs = sapply(substitute(list(...)),as.character) #includes a spurious "list" but who cares.
  ma = myargs
  print(myargs)
  if ("not" %in% myargs) {
    myargs = setdiff(c("par", "hpl", "it", "n", "fw", "dim", "val"),myargs)
  }
  if ("par" %in% myargs) {
    d = nicepar(d)
  }
  if ("hpl" %in% myargs) {
    print("hpl!!")
    d = nicehpl(d)
  }
  if ("it" %in% myargs) {
    d = niceit(d)
  }
  if ("n" %in% myargs) {
    if ("alt" %in% ma) {
      d = nicen2(d)
    } else {
      d = nicen(d)
    }
  }
  if ("fw" %in% myargs) {
    d = nicefw(d)
  }
  if ("dim" %in% myargs) {
    d = nicedim(d)
  }
  if ("val" %in% myargs) {
    d = niceval(d)
  }
  d
}
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl)
unique(somenice[,hpl])
unique(somenice[,histPerLoc])
sn = somenice[steps==2,]
dim(sn)
sn[model=="franken",]
somenice[,mean(meansq),by=list(model,dimension)]
somenice[steps>4,mean(meansq),by=list(model,dimension)]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                                             cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                       cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                                              cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                          cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl)
unique(somenice[,particles])
sn = somenice[steps==2,]
dim(sn)
View(sn)
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl,fw)
unique(somenice[,particles])
sn = somenice[steps==2,]
dim(sn)
sn[,list(model,dimension,params,cases)]
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],it,n,hpl,fw,dim)
unique(somenice[,particles])
sn = somenice[steps==2,]
dim(sn)
ggplot(data=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
arunmeans[steps==2,list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params,cases)]
caselist = arunmeans[steps==2,list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params,cases)]
caselist[useForward==1,]
caselist[useForward==1,]
caselist = arunmeans[steps==2,list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[useForward==1,]
caselist[useForward==1,]
substr(c("asdtasdt","q"),5)
substr(c("asdtasdt","q"),5,99)
runs[,mhType:=substr(mhType,5,99)]
runs[,sampType:=substr(sampType,5,99)]
runs[,params:=paste(sampType,mhType)]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                                             cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                       cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                                              cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                          cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
caselist = arunmeans[steps==2,list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[useForward==1,]
runs[,mhType:=substr(mhType,2,99)]
runs[,sampType:=substr(sampType,7,99)]
runs[,params:=paste(sampType,mhType)]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                                             cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                       cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                                              cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                          cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
caselist = arunmeans[steps==2,list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[useForward==1,]
runs[,mhType:=substr(mhType,2,99)]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runs[,params:=paste(sampType,mhType)]
runmeans = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                             meancov=mean(covdiv),sdcov=sd(covdiv),
                                             meandiff=mean(meandiv),sddiff=sd(meandiv),
                                             meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                                             cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                       cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1),list(meankl=mean(kl),sdkl=sd(kl),
                                              meancov=mean(covdiv),sdcov=sd(covdiv),
                                              meandiff=mean(meandiv),sddiff=sd(meandiv),
                                              meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                                              cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                          meancov=mean(covdiv),sdcov=sd(covdiv),
                          meandiff=mean(meandiv),sddiff=sd(meandiv),
                          meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                          cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
caselist = arunmeans[steps==2,list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[useForward==1,]
caselist = arunmeans[steps==2,list(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[useForward==1,]
caselist[useForward==1&particles>40,]
caselist[useForward %in% c(1,NA)&particles>40,]
nice = function(...,d=allrunmeans2) {
  myargs = sapply(substitute(list(...)),as.character) #includes a spurious "list" but who cares.
  ma = myargs
  #print(myargs)
  if ("not" %in% myargs) {
    myargs = setdiff(c("par", "hpl", "it", "n", "fw", "dim", "val"),myargs)
  }
  if ("par" %in% myargs) {
    d = nicepar(d)
  }
  if ("hpl" %in% myargs) {
    d = nicehpl(d)
  }
  if ("it" %in% myargs) {
    d = niceit(d)
  }
  if ("n" %in% myargs) {
    if ("alt" %in% ma) {
      d = nicen2(d)
    } else {
      d = nicen(d)
    }
  }
  if ("fw" %in% myargs) {
    d = nicefw(d)
  }
  if ("dim" %in% myargs) {
    d = nicedim(d)
  }
  if ("val" %in% myargs) {
    d = niceval(d)
  }
  if ("hit" %in% myargs) {
    d = d[nIter %in% c(NA,160),]
  }
  d
}
somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],hit,n,hpl,fw,dim)
unique(somenice[,particles])
sn = somenice[steps==2,]
dim(sn)
sn
ggplot(data=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
caselist[useForward %in% c(1,NA)&particles>40,]
caselist[order(model),][useForward %in% c(1,NA)&particles>40,]
?order
caselist[order(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
somenice = arunmeans[model %in% c("finkel","franken")&
                       dimensions==60 &
                       nIter %in% c(NA,120)&
                       particles %in% c(600,72000),]
somenice = arunmeans[model %in% c("finkel","franken")&
                       dimension==60 &
                       nIter %in% c(NA,120)&
                       particles %in% c(600,72000),]
sn = somenice[steps==2,]
dim(sn)
ggplot(data=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
somenice = somenice[steps>2,]
somenice = somenice[steps>2,]
ggplot(data=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
names(runs)
runs[,kl2:=kl]
runs[model=="franken",kl2:=klideal]
mean(runs[model=="franken",kl2-kl])
mean(runs[model=="franken",kl2-kl],na.rm=T)
runmeans = runs[useForward %in% c(NA,1)
                & steps>2,
                list(meankl=mean(kl),sdkl=sd(kl),
                     meancov=mean(covdiv),sdcov=sd(covdiv),
                     meandiff=mean(meandiv),sddiff=sd(meandiv),
                     meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                     cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
caselist = runmeans[steps==2,list(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[order(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
size(runmeans)
dim(runmeans)
caselist = runmeans[,list(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[order(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
runmeans2 = runs[useForward %in% c(NA,1)
                 & steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[steps>2,
                    list(meankl=mean(kl),sdkl=sd(kl),
                         meancov=mean(covdiv),sdcov=sd(covdiv),
                         meandiff=mean(meandiv),sddiff=sd(meandiv),
                         meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                         cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
caselist = runmeans2[,list(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[order(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
caselist = runmeans2[,list(model,dimension,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[order(model,dimension,steps,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
caselist[order(model,dimension,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
F
runs = fread("fixedest4.csv")
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
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
runs[model=="truth",mvlvar:=NA]
runs[,mhType:=substr(mhType,7,99)]
runs[,sampType:=substr(sampType,12,99)]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runs[,params:=paste(sampType,mhType)]
runs[,kl2:=kl]
runs[model=="franken",kl2:=klideal]
mean(runs[model=="franken",kl2-kl],na.rm=T)
runmeans = runs[useForward %in% c(NA,1)
                & steps>2,
                list(meankl=mean(kl),sdkl=sd(kl),
                     meancov=mean(covdiv),sdcov=sd(covdiv),
                     meandiff=mean(meandiv),sddiff=sd(meandiv),
                     meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                     cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1)
                 & steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[steps>2,
                    list(meankl=mean(kl),sdkl=sd(kl),
                         meancov=mean(covdiv),sdcov=sd(covdiv),
                         meandiff=mean(meandiv),sddiff=sd(meandiv),
                         meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                         cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
caselist = runmeans2[,list(model,dimension,particles,nIter,histPerLoc,useForward,params,cases)]
caselist[order(model,dimension,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
somenice = arunmeans[model %in% c("finkel","franken")&
                       dimension==30 &
                       nIter %in% c(NA,160)&
                       particles %in% c(400,32000),]
somenice = somenice[steps>2,]
ggplot(data=somenice,aes(x=sqrt(meankl),y=sqrt(meansq),shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
oldcl = caselist[order(model,dimension,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
if (F) { #peek at old runs
  runs = fread("fixedest4.csv")
}
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
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
if (F) { #draw "evol" graphs
  (pp1 = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
                aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
                )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
                ) + xlab("Time step") + ylab("Sum of loci 3-5")
   + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
   + theme_bw()
  )
  # tracks over time of all 3 algorithms
  (pp2 = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
                aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
                )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
                ) + xlab("Time step") + ylab("sum of loci 3:5")
    + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
    + theme_bw()
  )
  grid.arrange(pp1, pp2, ncol=1)
  #ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken")&nIter!=30,],aes(y=mvlmean,x=steps,color=paste(model,sampType,histPerLoc,nIter,useForward,mhType))) + geom_line()
  # calculate differences
  #getidealmean = Vectorize(function(astep,wid) {runs[worldid==wid&model=="ideal"&steps==astep,mvlmean]})
  runs[,idealmean:=mvlmean[match(T,(model=="ideal"))],by=list(worldid,steps)]
  runs[,err:=mvlmean-idealmean]
  #ggplot(runs[worldid==wid&model %in% c("finkel", "franken")],aes(y=err,x=steps,color=paste(model,sampType,histPerLoc,nIter))) + geom_line()
  runs[,sum(err^2),by=paste(model,sampType,histPerLoc,nIter)]
}
runs = runs[model!="",]
runs[,mhType:=substr(mhType,7,99)]
runs[,sampType:=substr(sampType,11,99)]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runs[,params:=paste(sampType,mhType)]
runs[,kl2:=kl]
runs[model=="franken",kl2:=klideal]
mean(runs[model=="franken",kl2-kl],na.rm=T)
runmeans = runs[useForward %in% c(NA,1)
                & steps>2,
                list(meankl=mean(kl),sdkl=sd(kl),
                     meancov=mean(covdiv),sdcov=sd(covdiv),
                     meandiff=mean(meandiv),sddiff=sd(meandiv),
                     meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                     cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1)
                 & steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[steps>2,
                    list(meankl=mean(kl),sdkl=sd(kl),
                         meancov=mean(covdiv),sdcov=sd(covdiv),
                         meandiff=mean(meandiv),sddiff=sd(meandiv),
                         meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                         cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
caselist = runmeans2[,list(model,dimension,particles,nIter,histPerLoc,useForward,params,cases)]
newcl = caselist[order(model,dimension,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
newcl
View(oldcl)
View(newcl)
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       nIter %in% c(NA,120)&
                       histPerLoc %in% c(NA,30)&
                       particles %in% c(400,32000),]
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
oldcl[order(histPerLoc),]
newcl[order(histPerLoc),]
oldcl[order(histPerLoc),]
max(oldcl[,histPerLoc])
max(newcl[,histPerLoc])
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       nIter %in% c(NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
arunmeans[model=="franken",nIter]
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       nIter %in% c(0,NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
sn = somenice[steps==3,]
dim(sn)
sn
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       useForward %in% c(NA,1)
                     nIter %in% c(0,NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       useForward %in% c(NA,1),
                     nIter %in% c(0,NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
somenice = somenice[steps>2,]
unique(somenice[,particles])
sn = somenice[steps==3,]
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       useForward %in% c(NA,1),
                     nIter %in% c(0,NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
somenice = somenice[steps>2,]
unique(somenice[,particles])
dim(somenice)
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       #useForward %in% c(NA,1),
                       nIter %in% c(0,NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
somenice = somenice[steps>2,]
dim(somenice)
sn = somenice[steps==3,]
dim(sn)
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       useForward %in% c(NA,1)&
                       nIter %in% c(0,NA,120)&
                       histPerLoc %in% c(0,NA,30)&
                       particles %in% c(400,32000),]
somenice = somenice[steps>2,]
dim(somenice)
sn = somenice[steps==3,]
dim(sn)
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
)
runmeans2[order(meankl),]
runs = fread("fixedest4.csv")
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
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
runs[model=="truth",mvlvar:=NA]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)
if (F) { #draw "evol" graphs
  (pp1 = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
                aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
                )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
                ) + xlab("Time step") + ylab("Sum of loci 3-5")
   + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
   + theme_bw()
  )
  # tracks over time of all 3 algorithms
  (pp2 = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
                aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
                )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
                ) + xlab("Time step") + ylab("sum of loci 3:5")
    + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
    + theme_bw()
  )
  grid.arrange(pp1, pp2, ncol=1)
  #ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken")&nIter!=30,],aes(y=mvlmean,x=steps,color=paste(model,sampType,histPerLoc,nIter,useForward,mhType))) + geom_line()
  # calculate differences
  #getidealmean = Vectorize(function(astep,wid) {runs[worldid==wid&model=="ideal"&steps==astep,mvlmean]})
  runs[,idealmean:=mvlmean[match(T,(model=="ideal"))],by=list(worldid,steps)]
  runs[,err:=mvlmean-idealmean]
  #ggplot(runs[worldid==wid&model %in% c("finkel", "franken")],aes(y=err,x=steps,color=paste(model,sampType,histPerLoc,nIter))) + geom_line()
  runs[,sum(err^2),by=paste(model,sampType,histPerLoc,nIter)]
}
runs = runs[model!="",]
runs[,mhType:=substr(mhType,7,99)]
runs[,sampType:=substr(sampType,11,99)]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runs[,params:=paste(sampType,mhType)]
runs[,kl2:=kl]
runs[model=="franken",kl2:=klideal]
mean(runs[model=="franken",kl2-kl],na.rm=T)
runmeans = runs[useForward %in% c(NA,1)
                & steps>2,
                list(meankl=mean(kl),sdkl=sd(kl),
                     meancov=mean(covdiv),sdcov=sd(covdiv),
                     meandiff=mean(meandiv),sddiff=sd(meandiv),
                     meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                     cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1)
                 & steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[steps>2,
                    list(meankl=mean(kl),sdkl=sd(kl),
                         meancov=mean(covdiv),sdcov=sd(covdiv),
                         meandiff=mean(meandiv),sddiff=sd(meandiv),
                         meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                         cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
runmeans2[order(meankl),]
runmeans2[order(meansq),]
head(runmeans2[order(meansq),],15)
head(runmeans2[order(meancov),],15)
runs = fread("fixedest5.2.csv")
runs = rbind(runs,fread("fixedest7.2.csv"))
runs = rbind(runs,fread("fixedest8.csv"))
if (F) { #peek at old runs
  runs = fread("fixedest4.csv")
}
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
runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
wid = 27
runs[,step:=steps-1]
runs[step==0@model!="truth",mvlmean:=NA]
runs[model %in% c("ideal")&step==0,mvlmean:=0]
runs[step==0,mvlvar:=NA]
runs[model!="finkel",nIter:=NA]
runs[model=="truth",mvlvar:=NA]
examplerun = runs[worldid==wid&nIter %in% c(NA,120),]
examplerun = rbind(examplerun,data.table(Algorithm="Finkelstein PF",step=0,mvlmean=0),fill=T)
if (F) { #draw "evol" graphs
  (pp1 = ggplot(examplerun[model %in% c("ideal", "observed", "truth"),],
                aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
                )) + geom_line() + ggtitle("Time series of truth, observation, and ideal Kalman filter"
                ) + xlab("Time step") + ylab("Sum of loci 3-5")
   + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
   + theme_bw()
  )
  # tracks over time of all 3 algorithms
  (pp2 = ggplot(examplerun[model %in% c("ideal", "finkel", "franken", "particle"),],
                aes(y=mvlmean,x=factor(step),color=Algorithm,group=Algorithm,linetype=Algorithm,
                )) + geom_line() + ggtitle("Time series using various filtering estimation algorithms"
                ) + xlab("Time step") + ylab("sum of loci 3:5")
    + geom_errorbar(aes(ymin=mvlmean-1.96*sqrt(mvlvar), ymax=mvlmean+1.96*sqrt(mvlvar)), width=.5, position=position_dodge(.2))
    + theme_bw()
  )
  grid.arrange(pp1, pp2, ncol=1)
  #ggplot(runs[worldid==wid&model %in% c("ideal", "finkel", "franken")&nIter!=30,],aes(y=mvlmean,x=steps,color=paste(model,sampType,histPerLoc,nIter,useForward,mhType))) + geom_line()
  # calculate differences
  #getidealmean = Vectorize(function(astep,wid) {runs[worldid==wid&model=="ideal"&steps==astep,mvlmean]})
  runs[,idealmean:=mvlmean[match(T,(model=="ideal"))],by=list(worldid,steps)]
  runs[,err:=mvlmean-idealmean]
  #ggplot(runs[worldid==wid&model %in% c("finkel", "franken")],aes(y=err,x=steps,color=paste(model,sampType,histPerLoc,nIter))) + geom_line()
  runs[,sum(err^2),by=paste(model,sampType,histPerLoc,nIter)]
}
runs = runs[model!="",]
runs[,mhType:=substr(mhType,7,99)]
runs[,sampType:=substr(sampType,11,99)]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]
runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]
runs[,params:=paste(sampType,mhType)]
runs[,kl2:=kl]
runs[model=="franken",kl2:=klideal]
mean(runs[model=="franken",kl2-kl],na.rm=T)
runmeans = runs[useForward %in% c(NA,1)
                & steps>2,
                list(meankl=mean(kl),sdkl=sd(kl),
                     meancov=mean(covdiv),sdcov=sd(covdiv),
                     meandiff=mean(meandiv),sddiff=sd(meandiv),
                     meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                     cases=.N),
                by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
arunmeans = runs[steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,useForward,params)]
runmeans2 = runs[useForward %in% c(NA,1)
                 & steps>2,
                 list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                      cases=.N),
                 by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
allrunmeans2 = runs[steps>2,
                    list(meankl=mean(kl),sdkl=sd(kl),
                         meancov=mean(covdiv),sdcov=sd(covdiv),
                         meandiff=mean(meandiv),sddiff=sd(meandiv),
                         meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T),
                         cases=.N),
                    by=list(model,Algorithm,dimension,particles,nIter,histPerLoc,sampType,mhType,useForward,params)] #no steps
runmeans2[order(meankl),]
runmeans2[order(meansq),]
head(runmeans2[order(meankl),],15)
head(runmeans2[order(meankl),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
head(runmeans2[order(meankl),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
head(runmeans2[order(meansq),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
head(runmeans2[order(meancov),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
head(runmeans2[order(meancov),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],30)
head(runmeans2[order(meankl-meandiff),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],30)
head(runmeans2[order(meandiff),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],30)
head(runmeans2[order(meansq),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
?savehistory
savehistory()
