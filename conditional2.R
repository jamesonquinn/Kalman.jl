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

#equal measurement variance edition

d = 15 #dimension

cond = fread("myworld3__260.csv")
cond = rbind(cond,fread("myworld3__261.csv"))
condi = cond[,]#nIter %in% c(NA,60),]
condi[,mean(kl),by=model ]
condi[,sqrt(mean(d11^2, na.rm=T)),by=list(model,step) ]

blocksize = 5
nblocks = floor(d/blocksize)


if (F) { #look at raw values
      truth = fread("nondestructiveWorld.csv")
      for (i in 1:d) {
        cat(i,
            mean((truth[,i,with=F]-truth[,i+d,with=F])[[1]]),
            var(truth[,i,with=F]-truth[,i+d,with=F]),
            measuresds[i],
            "\n")
      }
}


savecols = c(match("step",names(condi)),
             match("nIter",names(condi)),
             match("sampType",names(condi)),
             match("hpl",names(condi)))
pre_d = match("d01",names(condi)) - 1 #number of the last col before the "b's" μ
pre_v = match("v0101",names(condi)) - 1 #number of the last col before the "v's" Σ


biases = data.table()
for (l in 1:nblocks) {
  for (j in 1:blocksize) {
    oneterm = condi[,c(1,savecols,pre_d+j+blocksize*(l-1)),with=F] #model, step
    onevar = condi[,c(1,savecols,pre_v+j+blocksize*(l-1)),with=F]
    #print(names(oneterm))
    #oneterm[,loc:=l]
    names(oneterm)[length(names(oneterm))] = "error"
    names(onevar)[length(names(oneterm))] = "theestvar"
    oneterm[,estvar:=onevar[,theestvar]]
    #print(names(oneterm))
    biases = rbind(biases,oneterm[,list(bias=mean(error),variance=var(error),
                                               estvar=mean(estvar),
                                        block=l,position=j,
                                        border=(j %in% c(1,2,blocksize)),
                                        doubleborder=(j == 1)
                                        ),
                                  by=list(model,nIter,hpl,step,sampType)],
                    fill=T)
  }
}
#condi[,list(b01=mean(b01),k01=mean(k01),t01=mean(t01),o01=mean(o01),do=mean(do01)),by=list(model,step)]
allbiases = biases[,list(msqbias=sqrt(mean(bias^2,na.rm=T)),
             mvar=sqrt(mean(variance,na.rm=T)),
             mvar=sqrt(mean(estvar,na.rm=T)),
             border=border[1],
             doubleborder=doubleborder[1],
             anyborder=border[1] | doubleborder[1]
               ),by=list(model,position,block,step,nIter,hpl) ]
biases[,list(msqbias=mean(bias^2,na.rm=T),
             mvar=mean(variance,na.rm=T),
             mvar=mean(estvar,na.rm=T)
),by=list(model,nIter,hpl,sampType) ]
biases[,list(msqbias=sqrt(mean(bias^2,na.rm=T)),
             mvar=sqrt(mean(variance,na.rm=T)),
             mvar=sqrt(mean(estvar,na.rm=T))
),by=list(model,position) ]
oneterm[model=="Block PF"&step==1,error]


summary(lm(msqbias ~ border + doubleborder, data=allbiases[model=="Finkelstein PF",]))
summary(lm(msqbias ~ border + doubleborder, data=allbiases[model=="Block PF",]))









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
