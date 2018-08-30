library(ggplot2)
library(data.table)
library(gridExtra)
library(ggsci)
library(plyr)
library(tidyr)


cond = fread("outcomes_for_400.csv")
condi = cond[nIter %in% c(NA,160),]
condi[,mean(kl),by=model ]
condi[,mean(b11^2),by=list(model,step) ]
pre_b = 22
pre_v = 52
biases = data.table()
for (l in 1:10) {
  oneterm_left = condi[,c(1,14,pre_b-2+3*l),with=F]
  onevar = condi[,c(1,14,pre_v-2+3*l),with=F]
  #oneterm[,loc:=l]
  names(oneterm_left)[3] = "error"
  names(onevar)[3] = "theestvar"
  oneterm_left[,estvar:=onevar[,theestvar]]
  biases = rbind(biases,oneterm_left[,list(bias=mean(error),variance=var(error),
                                             estvar=estvar,neighborhood=l,loc=-2+3*l,border=T,position=-1),by=list(model,step)],
                  fill=T)
  oneterm_center = condi[,c(1,14,pre_b-1+3*l),with=F]
  onevar = condi[,c(1,14,pre_v-1+3*l),with=F]
  #oneterm[,loc:=l]
  names(oneterm_center)[3] = "error"
  names(onevar)[3] = "theestvar"
  oneterm_center[,estvar:=onevar[,theestvar]]
  biases = rbind(biases,oneterm_center[,list(bias=mean(error),variance=var(error),
                                               estvar=estvar,neighborhood=l,loc=-1+3*l,border=F,position=0),by=list(model,step)],
                  fill=T)
  oneterm_right = condi[,c(1,14,pre_b+3*l),with=F]
  onevar = condi[,c(1,14,pre_v+3*l),with=F]
  #oneterm[,loc:=l]
  names(oneterm_right)[3] = "error"
  names(onevar)[3] = "theestvar"
  oneterm_right[,estvar:=onevar[,theestvar]]
  biases = rbind(biases,oneterm_right[,list(bias=mean(error),variance=var(error),
                                              estvar=estvar,neighborhood=l,loc=3*l,border=F,position=0),by=list(model,step)],
                  fill=T)
}
biases[,mean(bias^2),by=list(model,border) ]




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
toPlot = data.table(gather(toPlotWide,model:border,sqerr,bias:estvar,factor_key=T))
names(toPlot)[3] = "component"
toPlot[,component:=relevel(component,"variance")]
toPlot[,position:=factor(border,levels=c(F,T),labels=c("central loci","peripheral loci"))]
toPlot2 = toPlot[component %in% c("bias","variance"),]
toPlot2[,component:=factor(component,levels=c("bias","variance"),labels=c("bias²","variance"))]
ggplot(data=toPlot2,
       aes(y=sqerr,x=model,fill=component)) +
  geom_col()+
  facet_wrap(~position)+
  labs(x="Algorithm",y="Squared error\nof distribution mean")+
  
