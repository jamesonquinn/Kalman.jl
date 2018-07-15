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
runs5 = fread("filtertest25.csv", fill=TRUE)
runs6 = fread("filtertest26.csv", fill=TRUE)

runs7 = fread("filtertest27.csv", fill=TRUE)
runs8 = fread("filtertest28.csv", fill=TRUE)

runs = rbind(runs1, runs2,runs3,runs4,runs5,runs6,runs7,runs8, fill=TRUE)
runs = runs[model!="finkel",] #code was broken, throw away results

runs9 = fread("filtertest31.csv", fill=TRUE)
runs10 = fread("filtertest32.csv", fill=TRUE)

runs = rbind(runs,runs9,runs10,fill=TRUE)


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



ggplot(data=runmeans[order(model,steps),],aes(x=meankl,y=meansq,shape=model,
                                              group=paste(model,particles,sampType,mhType,nIter,histPerLoc))
) + geom_path(aes(linetype=model)) + geom_point(aes(color=steps))

ggplot(data=runmeans[,],aes(x=meankl,y=meansq,shape=model,group=model)
) + geom_path(aes(linetype=model,color=model)
) + geom_point(
) + scale_x_continuous(limits=c(0,10)) + scale_y_continuous(limits=c(0,8))

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



