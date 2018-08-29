#starting line 5858
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