library(ggplot2)
library(data.table)
library(gridExtra)
library(ggsci)
library(plyr)
source("theme_publication.R")
pd = position_dodge(20)

runs = fread("fixedest9.csv")
runs = rbind(runs,fread("fixedest10.csv"))
runs = rbind(runs,fread("fixedest11.csv"))
if (F) { #temporarily commenting out
        runs = fread("fixedest5.2.csv")
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



runs[model != "finkel",nIter:=0]
unique(runs[,worldid])
# tracks over time of everything
wid = 27
runs[,step:=steps-1]
#runs[model=="observed",step:=steps]
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

#runs = rbind(runs,runs3)
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

head(runmeans2[order(meankl),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
head(runmeans2[order(meansq),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],15)
head(runmeans2[order(meancov),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],30)
head(runmeans2[order(meankl-meandiff),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],30)
head(runmeans2[order(meandiff),list(model,dimension,particles,nIter,histPerLoc,params,meankl,meancov,meansq)],30)

runmeans2[order(meankl),]
caselist = runmeans2[,list(model,dimension,particles,nIter,histPerLoc,useForward,params,cases)]
#oldcl = caselist[order(model,dimension,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]
newcl = caselist[order(model,dimension,particles,nIter,histPerLoc,useForward,params,cases),][useForward %in% c(1,NA)&particles>40,]


rf = runmeans[model=="finkel"|model=="ideal",]
rnf = runmeans2[model!="finkel",]

nicepar = function(d=allrunmeans2) {d[sampType %in% c("none","bkf.SampleLog") & d[,mhType] %in% c("none","bkf.MhCompromise"),]}
#nicedim = function(dim) {dim == 48}
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

#somenice = nice(d=arunmeans[model %in% c("finkel","franken"),],hit,n,hpl,fw,dim)
somenice = arunmeans[model %in% c("finkel","franken")&
                       #dimension==30 &
                       useForward %in% c(NA,1)&
                       nIter %in% c(NA,160)
                       &abs(meankl)<100
                       #&histPerLoc %in% c(0,NA,30)
                       #&particles %in% c(400,32000)
                     ,]
somenice2 = runs[model %in% c("finkel","franken")&
                       #dimension==30 &
                       useForward %in% c(NA,1)&
                       nIter %in% c(NA,160)&steps>5
                    &abs(kl)<100
                     #&histPerLoc %in% c(0,NA,30)
                     #&particles %in% c(400,32000)
                     ,]
somenice3 = somenice2[,list(meankl=mean(kl),sdkl=sd(kl),
                            meancov=mean(covdiv),sdcov=sd(covdiv),
                            meandiff=mean(meandiv),sddiff=sd(meandiv),
                            meantime=mean(runtime),meansq=mean(sqerr,na.rm=T),
                            cases=.N),by=list(model,params,worldid)]
if (F) {
somenice = arunmeans[model %in% c("finkel","franken")&
                       dimension==60 &
                       nIter %in% c(NA,120)&
                       particles %in% c(600,72000),]
}

somenice = somenice[steps>5,]
dim(somenice)
sn = somenice[steps==8,]
dim(sn)
ggplot(data=somenice,aes(x=meankl,y=meansq,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
) 

ggplot(data=somenice,aes(x=meankl-meandiff,y=meandiff,shape=model)
) + geom_path(aes(linetype=params,color=params)
) + geom_point(aes(color=params)
) 

ggplot(data=somenice3,aes(x=meankl,y=meansq,shape=model)
) + geom_point(aes(color=params)
) 

ggplot(data=somenice2,aes(x=kl,y=sqerr,shape=model,group=paste(model,worldid))
) + geom_path(aes(linetype=factor(worldid),color=params)
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

qplot(meansq,meancov,data=nice(not,it,fw,dim,d=runmeans),shape=paste(useForward),size=1,linetype=paste(useForward,nIter),group=paste(useForward,nIter),main=
        "finkel by nIter useForward allsteps") + geom_path(size=1)

qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(sampType),
      main="finkel by sampType zoom")

qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(mhType),
      main="finkel by mhType zoom")
qplot(meankl,meansq,data=runmeans2[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(mhType),
      main="finkel by mhType zoom")

qplot(meankl,meansq,data=nice(not,hpl),shape=model,color=paste(histPerLoc,nIter),
      main="finkel by histPerLoc zoom")
qplot(meankl,meansq,data=runmeans2[meankl<5&meankl>0&meansq<50,],shape=model,color=paste(histPerLoc,nIter),
      main="finkel by histPerLoc zoom")

qplot(meankl,meansq,data=runmeans[meankl<5&meankl>0&meansq<50,],shape=model,color=as.factor(paste(sampType,mhType)),
      main="finkel by samp/mh zoom")
qplot(meankl,meansq,data=rf[meankl<5&sampType=="bkf.SampleLog"&mhType=="bkf.MhCompromise",],shape=as.factor(nIter),color=as.factor(particles),size=steps,
      main="log/compromise by nIter/particles")
qplot(meankl,meandiff,data=runmeans[meankl<50&(model!="finkel"|(sampType=="bkf.SampleLog"&mhType=="bkf.MhCompromise")),],shape=model,color=paste(as.factor(pmin(particles,1001)),nIter),
      main="nonfinkel plus log/compromise by nIter/particles")





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
ggplot(data=nice(not,it,d=arunmeans),aes(x=steps,y=meansq,shape=factor(nIter),
                                              group=paste(model,particles,sampType,mhType,nIter,histPerLoc))
) + geom_path(aes(linetype=model)) + geom_point(aes(color=steps)) 


ggplot(data=equilib,aes(x=sqrt(kl),y=sqrt(sq),color=paste(model,sampType),shape=mhType)) + geom_point()

ggplot(data=runmeans[order(model,steps),],aes(x=meankl,y=meansq,shape=model,
                                              group=paste(model,particles,sampType,mhType,nIter,histPerLoc))
) + geom_path(aes(linetype=model)) + geom_point(aes(color=steps))

ggplot(data=arunmeans[steps>4&model %in% c("finkel","franken")&particles>40&!(nIter %in% 0),],aes(x=meankl,y=meancov,shape=model,group=model)
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



