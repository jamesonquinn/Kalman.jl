library(ggplot2)
library(data.table)

pd = position_dodge(20)

runs = fread("filtertest60incomplete.csv")
runs[,i:=.I]
starts = runs[model=="truth"&rep==1&steps==1,i]
worldnum = Vectorize(function(i) {sum(i>=starts)})
runs[,worldid:=worldnum(i)]



runs[model != "finkel",nIter:=0]
# tracks over time of everything
ggplot(runs[worldid==1,],aes(y=mvlmean,x=steps,color=paste(model,sampType,histPerLoc))) + geom_line() 

# tracks over time of all 3 algorithms
ggplot(runs[model %in% c("ideal", "finkel", "franken")],aes(y=mvlmean,x=steps,color=paste(model,sampType,histPerLoc,nIter))) + geom_line() 

# calculate differences
getidealmean = Vectorize(function(astep) {runs[model=="ideal"&steps==astep,mvlmean]})
runs[,idealmean:=getidealmean(steps)]
runs[,err:=mvlmean-idealmean]
ggplot(runs[model %in% c("finkel", "franken")],aes(y=err,x=steps,color=paste(model,sampType,histPerLoc,nIter))) + geom_line() 
runs[,sum(err^2),by=paste(model,sampType,histPerLoc,nIter)]





runs[sampType=="log7.5_20",`:=`(sampType="bkf.SampleLog",mhType="bkf.MhSampled")]
runs[sampType=="compromise2.1..log7.5_20",`:=`(sampType="bkf.SampleLog",mhType="bkf.MhCompromise")]
runs[sampType=="compromise2.1..uniform",`:=`(sampType="bkf.SampleUniform",mhType="bkf.MhCompromise")]
runs[sampType=="sampled..uniform",`:=`(sampType="bkf.SampleUniform",mhType="bkf.MhSampled")]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]




runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]

#runs = rbind(runs,runs3)
runmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType)]
runmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                by=list(model,particles,nIter,histPerLoc,sampType,mhType)] #no steps




rf = runmeans[model=="finkel"|model=="ideal",]
rnf = runmeans2[model!="finkel",]



qplot(meankl,meansq,data=runmeans[meankl<50&meankl>0&meansq<50,],shape=model,color=steps,main="all from run 30")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],color=model,main="zoom 30 colors")


qplot(meankl,meansq,data=runmeans[meankl<50&meankl>0&meansq<50,],shape=model,color=as.factor(particles),main="finkel by particles")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(particles),main="finkel by particles zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(steps),main="finkel by steps zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(nIter),main="finkel by nIter zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(sampType),main="finkel by sampType zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(mhType),main="finkel by mhType zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(histPerLoc),main="finkel by histPerLoc zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(paste(sampType,mhType)),main="finkel by samp/mh zoom")
qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(paste(sampType,mhType)),main="finkel by samp/mh zoom")
qplot(meankl,meansq,data=rf[meankl<5&sampType=="bkf.SampleLog"&mhType=="bkf.MhCompromise",],shape=as.factor(nIter),color=as.factor(particles),size=steps,main="log/compromise by nIter/particles")
qplot(meankl,meansq,data=runmeans[meankl<50&(model!="finkel"|(sampType=="bkf.SampleLog"&mhType=="bkf.MhCompromise")),],shape=model,color=as.factor(pmin(particles,1001)),main="nonfinkel plus log/compromise by nIter/particles")





runs[is.na(sampType),sampType:="NA"]
runs[model !="finkel",sampType:="none"]

runs[is.na(histPerLoc),histPerLoc:=0]
runs[is.na(nIter),nIter:=0]

times= lm(runtime ~ model + nIter + sampType + particles + particles*model + I(particles**2*(model=="finkel")), runs)
times= lm(runtime ~ model + nIter + histPerLoc + I(histPerLoc*nIter) + sampType + particles + particles*model + I(particles**2*(model=="finkel")) + I(histPerLoc*nIter*(sampType=="compromise2.1..log7.5_20")), runs)
summary(times)

runs[,invpart:=1/sqrt(particles)]
runs[,invIter1:=exp(-nIter/50)]
runs[,invIter2:=exp(-nIter/500)]
runs[,invIter1:=exp(-nIter/50)]
runs[,invIter2:=exp(-nIter/500)]
runs[,lkl:=log(kl)]
errors = lm(lkl ~ model + invpart + invpart*model + histPerLoc + sampType + invIter1 + invIter2,runs)
summary(errors)

nona = copy(runs)
nona[is.na(histPerLoc),histPerLoc:=0]

table(rf[,list(meankl<10,sampType)])
table(rf[,list(meankl<10,histPerLoc)])
table(rf[,list(meankl<10,steps)])
table(rf[,list(meankl<10,nfp)])
table(runmeans[,list(meankl<10,model)])



























#for 33 and on

runs = fread("filtertest37.csv", fill=TRUE)
runs = runs[model=="finkel",]
runs = rbind(runs,fread("filtertest36.csv", fill=TRUE))
runs = rbind(runs,fread("filtertest38.csv", fill=TRUE))
runs = rbind(runs,fread("filtertest39.csv", fill=TRUE))
#runs = rbind(runs,fread("filtertest40.csv", fill=TRUE))
runs = rbind(runs,fread("filtertest41.csv", fill=TRUE))
runs = rbind(runs,fread("filtertesty.csv", fill=TRUE))


runs2 = fread("filtertest40.csv", fill=TRUE)
runs = rbind(runs,runs2)

#above is largely focused on 400/.85. Now, some 600/.7:
runs = fread("filtertesty.csv", fill=TRUE)




runs[sampType=="log7.5_20",`:=`(sampType="bkf.SampleLog",mhType="bkf.MhSampled")]
runs[sampType=="compromise2.1..log7.5_20",`:=`(sampType="bkf.SampleLog",mhType="bkf.MhCompromise")]
runs[sampType=="compromise2.1..uniform",`:=`(sampType="bkf.SampleUniform",mhType="bkf.MhCompromise")]
runs[sampType=="sampled..uniform",`:=`(sampType="bkf.SampleUniform",mhType="bkf.MhSampled")]
runs[is.na(sampType),`:=`(sampType="none",mhType="none")]

runs[,params:=paste(sampType,mhType)]




runs[model=="particle",particles:=npf]
runs[model=="franken",particles:=nfapf]
runs[model=="finkel",particles:=nfp]
runs[model=="ideal",particles:=1]

#runs = rbind(runs,runs3)
runmeans = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                      meancov=mean(covdiv),sdcov=sd(covdiv),
                      meandiff=mean(meandiv),sddiff=sd(meandiv),
                      meantime=mean(runtime),meansq=mean(sqerr,na.rm=T)),
                by=list(model,dimension,steps,particles,nIter,histPerLoc,sampType,mhType,params)]
runmeans2 = runs[,list(meankl=mean(kl),sdkl=sd(kl),
                       meancov=mean(covdiv),sdcov=sd(covdiv),
                       meandiff=mean(meandiv),sddiff=sd(meandiv),
                       meantime=mean(runtime,na.rm=T),meansq=mean(sqerr,na.rm=T)),
                 by=list(model,particles,nIter,histPerLoc,sampType,mhType,params)] #no steps




rf = runmeans[model=="finkel"|model=="ideal",]
rnf = runmeans2[model!="finkel",]

equilib = runmeans[steps>7,list(kl=mean(meankl),sq=mean(meansq)),by=list(model,params,sampType,mhType)]

# x=time, y=meansq
ggplot(data=runmeans[order(model,steps),],aes(x=steps,y=meansq,shape=model,
                                              group=paste(model,particles,sampType,mhType,nIter,histPerLoc))
) + geom_path(aes(linetype=model)) + geom_point(aes(color=steps)) 


ggplot(data=equilib,aes(x=sqrt(kl),y=sqrt(sq),color=paste(model,sampType),shape=mhType)) + geom_point()

ggplot(data=runmeans[order(model,steps),],aes(x=meankl,y=meansq,shape=model,
                                              group=paste(model,particles,sampType,mhType,nIter,histPerLoc))
) + geom_path(aes(linetype=model)) + geom_point(aes(color=steps))

ggplot(data=runmeans[,],aes(x=sqrt(meankl),y=sqrt(meansq),shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
) 

ggplot(data=runmeans[model=="finkel"&nIter>10,],aes(x=sqrt(meankl),y=sqrt(meansq),shape=params,group=params)
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



