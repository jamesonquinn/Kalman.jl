#line 3427 is good? No, way too early.
#search for + geom_path(aes(linetype=params,color=params)
#Line 5924 looks promising?
#search backwards for fread
#5859-5924??

runs = fread("fixedest4.csv")
=.I]
starts = runs[model=="truth"&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
=worldnum(i)]
=0]
unique(runs[,worldid])
wid = 4
=steps-1]
=`(sampType="none",mhType="none")]
=npf]
=nfapf]
=nfp]
=1]
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