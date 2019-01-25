#Analyze nonlinear output conditional on fixed truth and observations

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(data.table)
library(gridExtra)
library(ggsci)
library(plyr)
library(tidyr)
library(stringr)
shifter = function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}
#read data from different formats and merge to common format

d = 20 #dimension

cond = fread("1outcomes_big_400.csv")
for (i in 2:6) {
  cond2 = fread(paste0(i,"outcomes_big_400.csv"))
  cond2[,rep:=rep + (100*i)]
  cond = rbind(cond, cond2)
}

cond = cond[model=="Finkelstein PF"]

for (i in 1:2) {
  cond2 = fread(paste0(i,"outcomes_wideblock_400.csv"))
  cond2[,rep:=rep + 1000+(100*i)]
  cond = rbind(cond, cond2, fill=T)
}

for (i in 2:3) {
  cond2 = fread(paste0("1outcomes_wideblock",i,"_400.csv"))
  cond2[,rep:=rep + 1000*i]
  cond = rbind(cond, cond2, fill=T)
}

condi = cond[nIter %in% c(NA,160),]
condi[,mean(kl),by=model ]
condi[,sqrt(mean(b11^2)),by=list(model,step) ]
tcol = match("t01",names(condi))
condi = condi[,1:(tcol-1),with=F]

blocksize = 5
nblocks = floor(d/blocksize)
measurevars = rep(c(0.36,0.09,1.,4.,   1.,.09,4.,0.36),10)[1:d]
measuresds = sqrt(measurevars)
peersds = shifter(measuresds,-1)
quadneighborsds = shifter(measuresds,-2) *  shifter(measuresds,1)



truth = fread("myLorenzWorldVecObsVar.csv")
for (i in 1:d) {
  cat(i,
      mean((truth[,i,with=F]-truth[,i+d,with=F])[[1]]),
      var(truth[,i,with=F]-truth[,i+d,with=F]),
      measuresds[i],
      "\n")
}

allcols= fread("1outcomes_wideblock3_400.csv")[1:10,]
names(allcols)
tcol = match("t01",names(allcols))
ncol = dim(allcols)[2]
condi = cbind(condi,allcols[condi[,step],tcol:ncol,with=F])
names(condi)
for (i in 1:d) {
  dname = paste0("d",str_pad(i, 2, pad = "0"))
  doname = paste0("do",str_pad(i, 2, pad = "0"))
  bname = paste0("b",str_pad(i, 2, pad = "0"))
  kname = paste0("k",str_pad(i, 2, pad = "0"))
  tname = paste0("t",str_pad(i, 2, pad = "0"))
  oname = paste0("o",str_pad(i, 2, pad = "0"))
  dexpr <- bquote(.(as.name(dname)):=
                   .(as.name(bname))+.(as.name(kname))-.(as.name(tname)))
  doexpr <- bquote(.(as.name(doname)):=
                   .(as.name(bname))+.(as.name(kname))-.(as.name(oname)))
  condi[,eval(dexpr)]
  condi[,eval(doexpr)]
}



stepcol = match("step",names(condi))
pre_d = match("do01",names(condi)) - 1 #number of the last col before the "b's" μ
pre_v = match("v0101",names(condi)) - 1 #number of the last col before the "v's" Σ


biases = data.table()
for (l in 1:nblocks) {
  for (j in 1:blocksize) {
    oneterm = condi[,c(1,stepcol,pre_d+j+blocksize*(l-1)),with=F] #model, step
    onevar = condi[,c(1,stepcol,pre_v+j+blocksize*(l-1)),with=F]
    #oneterm[,loc:=l]
    names(oneterm)[3] = "error"
    names(onevar)[3] = "theestvar"
    oneterm[,estvar:=onevar[,theestvar]]
    oneterm[,mysd:=measuresds[j+blocksize*(l-1)]]
    oneterm[,peersd:=peersds[j+blocksize*(l-1)]]
    oneterm[,quadsd:=quadneighborsds[j+blocksize*(l-1)]]
    biases = rbind(biases,oneterm[,list(bias=mean(error),variance=var(error),
                                               estvar=mean(estvar),
                                        block=l,position=j,
                                        border=(j %in% c(1,2,blocksize)),
                                        doubleborder=(j == 1)),
                                  by=list(model,step)],
                    fill=T)
  }
}
condi[,list(b01=mean(b01),k01=mean(k01),t01=mean(t01),o01=mean(o01),do=mean(do01)),by=list(model,step)]
biases[,list(msqbias=sqrt(mean(bias^2,na.rm=T)),
             mvar=sqrt(mean(variance,na.rm=T)),
             mvar=sqrt(mean(estvar,na.rm=T))
               ),by=list(model,position,block) ]
biases[,list(msqbias=sqrt(mean(bias^2,na.rm=T)),
             mvar=sqrt(mean(variance,na.rm=T)),
             mvar=sqrt(mean(estvar,na.rm=T))
),by=list(model) ]
oneterm[model=="Block PF"&step==1,error]












cond2 = fread("outcomes_big_400.csv")
condi2 = cond2[nIter %in% c(NA,160),]
condi2[,mean(kl),by=model ]
condi2[,mean(b11^2),by=list(model,step) ]
pre_b = 22
pre_v = 22+66
names(condi2)[pre_v+1]
biases2 = data.table()
for (l in 1:22) {
  oneterm_left = condi2[,c(1,14,pre_b-2+3*l),with=F]
  onevar = condi2[,c(1,14,pre_v-2+3*l),with=F]
  #oneterm[,loc:=l]
  names(oneterm_left)[3] = "error"
  names(onevar)[3] = "theestvar"
  oneterm_left[,estvar:=onevar[,theestvar]]
  biases2 = rbind(biases2,oneterm_left[,list(bias=mean(error),variance=var(error),
                                             estvar=estvar,neighborhood=l,loc=-2+3*l,border=T,position=-1),by=list(model,step)],
                 fill=T)
  oneterm_center = condi2[,c(1,14,pre_b-1+3*l),with=F]
  onevar = condi2[,c(1,14,pre_v-1+3*l),with=F]
  #oneterm[,loc:=l]
  names(oneterm_center)[3] = "error"
  names(onevar)[3] = "theestvar"
  oneterm_center[,estvar:=onevar[,theestvar]]
  biases2 = rbind(biases2,oneterm_center[,list(bias=mean(error),variance=var(error),
                                               estvar=estvar,neighborhood=l,loc=-1+3*l,border=F,position=0),by=list(model,step)],
                 fill=T)
  oneterm_right = condi2[,c(1,14,pre_b+3*l),with=F]
  onevar = condi2[,c(1,14,pre_v+3*l),with=F]
  #oneterm[,loc:=l]
  names(oneterm_right)[3] = "error"
  names(onevar)[3] = "theestvar"
  oneterm_right[,estvar:=onevar[,theestvar]]
  biases2 = rbind(biases2,oneterm_right[,list(bias=mean(error),variance=var(error),
                                              estvar=estvar,neighborhood=l,loc=3*l,border=F,position=0),by=list(model,step)],
                 fill=T)
}
biases2[,list(bias=mean(bias^2),variance=mean(variance),sqerr=mean(bias^2)+mean(variance),estvar=mean(estvar)),by=list(model,border) ][order(model),]

#biases2[,mean(bias^2),by=list(model,border) ]
biases[,size:=30]
biases2[,size:=66]
bcomb = rbind(biases,biases2)

bcomb[,list(bias=mean(bias^2),variance=mean(variance),sqerr=mean(bias^2)+mean(variance),estvar=mean(estvar)),by=list(model,border,size) ][order(model),]
toPlotWide = bcomb[,list(bias=mean(bias^2),variance=mean(variance),sqerr=mean(bias^2)+mean(variance),estvar=mean(estvar)),by=list(model,border) ][order(model),]
#biases[,mean(bias^2),by=list(model,border,neighborhood) ][order(model),]
#biases[,mean(bias^2),by=list(loc,model)][order(model),]
toPlot3wide= bcomb[,list(bias=mean(bias^2),variance=mean(variance),sqerr=mean(bias^2)+mean(variance),estvar=mean(estvar)),by=list(model) ][order(model),]
toPlot3wide[,border:="All loci"]
toPlotWide = rbind(toPlotWide, toPlot3wide)
toPlot = data.table(gather(toPlotWide,model:border,sqerr,bias:estvar,factor_key=T))
names(toPlot)[3] = "component"
toPlot[,component:=relevel(component,"variance")]
toPlot[,position:=factor(border,levels=c("All loci","FALSE","TRUE"),labels=c("All loci","Central loci","Peripheral loci"))]
toPlot2 = toPlot[component %in% c("bias","variance"),]
toPlot2[,component:=factor(component,levels=c("variance","bias"),labels=c("variance","bias²"))]

ggplot(data=toPlot2,
       aes(y=sqerr,x=model,fill=component)) +
  geom_col()+
  facet_wrap(~position)+
  labs(x="Algorithm",y="Average squared error per locus\nof estimated filtering distribution mean")+
  theme_bw() + scale_fill_uchicago()
