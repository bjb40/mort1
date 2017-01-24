#dev R 3.2.3; Rstan 2.8.2; RStudio 0.99.0467; Stan 2.8.0
#figures and tables for changepoint model
#Bryce Bartlett
#2/2016


rm(list=ls())

#@@@@@@
#Universals
#@@@@@@

#functions
source('funs.R')

#helper function for tables
eff = function(s,c){
  # calculates mean effect from posterior sample and 95%CI
  #
  # Args:
  #   s: a series of posterior samples
  #   c: a number betwen 0 and 1 for confidence; .95 is default
  #
  # Returns:
  #   a named vector including mean, lower ci, upper ci

  e = mean(s)
  names(e)='mean'
  #calculate probs
  if(missing(c)){c=.95}
  l=(1-c)/2
  u=1-((1-c)/2)

  interval = quantile(s,prob=c(l,u))
  return(c(e,interval))
}


#helper function for printing
printeff = function(s,c){
  #returns string of the form "mean posterior [lower post,upper post]"
  #Input: s = a series of posterior samples,c = confidence pval (0-1)
  #see function (eff) for more info
  #output is string
  e = rnd(eff(s,c))
  return(cat(e[1],'<br>[',paste(e[2],e[3],sep=','),']',sep=''))
}


testlower = function(s){
  #returns bayesian pvalue to test hypotheses that ICD10 RR are less than ICD9 RR
  #Input: s = posterior samples, with column1 from ICD9 and column 2 from ICD10
  #output is a single number
  #not strictly correct; need a yrep or maybe \tilde(rep)

  return(sum(apply(s,1,FUN=function(x) x[2]<x[1]))/nrow(s))

}

#calculate CI of difference in estimated effects
delta=function(s,c,p){
  #returns distribution of the difference between ICD9 and ICD10 effects
  #Input: s = posterior samples, with column1 from ICD9 and column 2 from ICD10
  #input: c = confidence value (see function eff for information)
  #input: p is True or false for printing
  #output is a vecotr of 3 values, mean, and CI
  #not strictly correct; need a yrep or maybe \tilde(rep)

  if(missing(p)){p=FALSE}
  #print(length(p))
  
  d = apply(s,1,FUN=function(x) x[2]-x[1])

  if(p==FALSE){return(eff(d,c))}
  if(p==TRUE){return(printeff(d,c))}
}

#scaled student t for robust functions
rscaled_t = function(n,df,s){
  # returns a vector of n length from a scaled student t with degrees of freedem df, and
  # based on a mixture of normals, see BDA 3 pp. 294,576
  #
  # Args:
  #   n: number of samples to draw
  #   df: degrees of freedom
  #   s: scale (equivalent to sd for normal draw)
  #
  # Returns:
  #   a named vector n length constituting a random draw from the scaled t distribution
  alpha=df/2
  beta=alpha*s
  v_i = 1/rgamma(n,shape=alpha,rate=beta)
  return(rnorm(1000,mean=0,sd=v_i))
}


rawdir = "H:/projects/mort1/dat~/"
outdir = "H:/projects/mort1/output/"
imdir = "H:/projects/mort1/img~/"

#@@@@
#Load and clean data
#@@@@

dat = read.csv(paste0(outdir,'stata-series.csv'))

model=list()
#load previously sved datasets from time_series_analysis2.R
for(m in 1:4){
  load(file=paste0(outdir,paste0('m',m,'samp.gz')),verbose=T)
  model[[m]] = samp; rm(samp)
}


#@@@@@@@@@@@@@@@@@@@
#Put together 2 tables of results
#@@@@@@@@@@@@@@@@@@@

gammanames = c(as.character(seq(40,85,by=5)),'Female','Complex','Home','Ltcare','Oplace','Black')

#prepare yrrac model (models 1 and 2)

sink(paste0(outdir,'yrrac_regtable.txt'))

cat('\n\nTable __. Bayesian regression estimates on logged Relative Rate
    of Acute to Chronic all-cause deaths 1994-2003.\n\n')


#need to check with scott on the "pvalue" test, and read is simulation articles...
#need to integrate uncertainty in measure out of this, probably to get a a \tilde(y)

cat('|      |  Model 1 | Model 2-ICD9 | Model 2-ICD10 | M2: ICD10-ICD9 |\n')
cat('|:-----|---------:|-------------:|--------------:|--------------------------------:|\n')

#year effect - Beta
cat('| Year | ')
cat(printeff(model[[1]]$beta), '|')
cat(printeff(model[[2]]$beta[,1]), '|')
cat(printeff(model[[2]]$beta[,2]), '|')

#bayesian pvalue ICD10 effect < ICD9 effect
#cat(delta(model[[2]]$beta,p=TRUE),'|\n')

#conditional means of effects - gammas
for(e in 1:16){
  cat('|',gammanames[e],'|')
  cat(printeff(model[[1]]$gamma[,e]), '|')
  #model 2 is three dimensional
  cat(printeff(model[[2]]$gamma[,1,e]), '|')
  cat(printeff(model[[2]]$gamma[,2,e]), '|')
  cat(delta(model[[2]]$gamma[,,e],p=TRUE),'|\n')

}

#variances

cat('| Level 2 Cell Var ($\\Sigma$) |')
cat(printeff(model[[1]]$zi),'|')
cat(printeff(model[[2]]$zi[,1]), '|')
cat(printeff(model[[2]]$zi[,2]), '|')
cat(delta(model[[2]]$zi,p=TRUE),'|\n')
cat(delta)

cat('| Level 2 Cell Var ($\\delta^2$) |')
cat(printeff(model[[1]]$delta),'|')
cat(printeff(model[[2]]$delta), '|')
cat(' |')
cat(' |\n')


cat('| Correlation of Level 2 Var (Corr($\\alpha_1,\\alpha_2$)) |   -   |    |')
#calculated from L_Omega which is the lower cholesky of the correlation matrix
# so it is L_Omega %*% t(L_Omega) on the upper or lower diagonal

corr = apply(model[[2]]$L_Omega,1,FUN = function(x) (matrix(x,2,2) %*% t(matrix(x,2,2)))[1,2])
cat(printeff(corr), '| -  |\n')

cat('| Level 1 Var |')
cat(printeff(model[[1]]$sig),'|')
cat(printeff(model[[2]]$sig), '|   |  - |')

cat('\n\nNote: Mean estimates with 95% CI; fit can be input from model summaries')

sink()




#prepare yrrdc model (models 3 and 4)

sink(paste0(outdir,'yrrdc_regtable.txt'))

cat('\n\nTable __. Bayesian regression estimates on logged Relative Rate
    of Acute to Chronic underlying-cause deaths 1994-2003.\n\n')


#need to check with scott on the "pvalue" test, and read is simulation articles...
#need to integrate uncertainty in measure out of this, probably to get a a \tilde(y)

cat('|      |  Model 3 | Model 4-ICD9 | Model 4-ICD10 | M4: ICD10-ICD9 |\n')
cat('|:-----|---------:|-------------:|--------------:|--------------------------------:|\n')

#year effect - Beta
cat('| Year | ')
cat(printeff(model[[3]]$beta), '|')
cat(printeff(model[[4]]$beta[,1]), '|')
cat(printeff(model[[4]]$beta[,2]), '|')

#bayesian pvalue ICD10 effect < ICD9 effect
cat(delta(model[[4]]$beta,p=TRUE),'|\n')

#conditional means of effects - gammas
for(e in 1:16){
  cat('|',gammanames[e],'|')
  cat(printeff(model[[3]]$gamma[,e]), '|')
  #model 2 is three dimensional
  cat(printeff(model[[4]]$gamma[,1,e]), '|')
  cat(printeff(model[[4]]$gamma[,2,e]), '|')
  cat(delta(model[[4]]$gamma[,,e],p=TRUE),'|\n')

}

#variances

cat('| Level 2 Var |')
cat(printeff(model[[3]]$zi),'|')
cat(printeff(model[[4]]$zi[,1]), '|')
cat(printeff(model[[4]]$zi[,2]), '|')
cat(delta(model[[4]]$zi,p=TRUE),'|\n')

cat('| Level 2 Cell Var ($\\delta^2$) |')
cat(printeff(model[[1]]$delta),'|')
cat(printeff(model[[2]]$delta), '|')
cat(' |')
cat(' |\n')


cat('| Correlation of Level 2 Var |   -   |    |')
#calculated from L_Omega which is the lower cholesky of the correlation matrix
# so it is L_Omega %*% t(L_Omega) on the upper or lower diagonal

corr = apply(model[[4]]$L_Omega,1,FUN = function(x) (matrix(x,2,2) %*% t(matrix(x,2,2)))[1,2])
cat(printeff(corr), '| -  |\n')

cat('| Level 1 Var |')
cat(printeff(model[[3]]$sig),'|')
cat(printeff(model[[4]]$sig), '|   |  - |')

cat('\n\nNote: Mean estimates with 95% CI; fit can be input from model summaries')

sink()



#@@@@@@@@@@@@@@@@@@@
#Age-group ppd - male and female
#would need reweighted to look right...
#@@@@@@@@@@@@@@@@@@@

#slices
datuc = dat[is.finite(dat$yrrdc),]
slices=unique(dat[,c('Female','Black')])
slices.names = c('White Male','White Female','Black Male','Black Female')
#visual inspection - need to triple check
plts = array(0,dim=c(3,2,2,4))
dimnames(plts)=list(est=c('mean','lcl','ucl'),
                  dv=c('all-cause','underlying-cause'),
                  icd=c('9','10'),dem=slices.names)

#helper function for limiting
makelim = function(dt){
  #input is dataframe
  #returns fector of true/false for limits
  return(list(Black=dt[,'Black'] == 1,
              Female=dt[,'Female'] == 1,
              icd10=dt$Years>=0,
              tdeaths=dt$tdeaths))
  
}

ppd.ac = cbind(as.data.frame(makelim(dat)),t(model[[2]]$ppd))
ppd.uc = cbind(as.data.frame(makelim(datuc)),t(model[[4]]$ppd))
ppd=list(ac=ppd.ac,uc=ppd.uc); rm(ppd.ac,ppd.uc)

deltas=data.frame(dv=integer(),
                  mean=double(),
                  upper=double(),
                  lower=double(),
                  Female=integer(),
                  Black=integer())

#exponentiate to put on ratio scale
#ppd=exp(ppd)

#par(mfrow=c(1,4))
for(r in 1:nrow(slices)){
  for(dv in 1:2){
      l9 = ppd[[dv]]$Black == slices[r,'Black'] &
              ppd[[dv]]$Female == slices[r,'Female'] &
              ! ppd[[dv]]$icd10
      
      l10 = ppd[[dv]]$Black == slices[r,'Black'] &
            ppd[[dv]]$Female == slices[r,'Female'] &
            ppd[[dv]]$icd10
      
      wt.9 = ppd[[dv]]$tdeaths[l9]/sum(ppd[[dv]]$tdeaths[l9])
      wt.10 = ppd[[dv]]$tdeaths[l10]/sum(ppd[[dv]]$tdeaths[l10])
        
      tmp9 = apply(ppd[[dv]][l9,paste0(1:5000)],2,FUN=function(x) sum(x*wt.9))  
      tmp10 = apply(ppd[[dv]][l10,paste0(1:5000)],2,FUN=function(x) sum(x*wt.10))

      tmpdeltas=eff(exp(tmp10)-exp(tmp9))
      names(tmpdeltas) = c('mean','lower','upper')
      #print(tmpdeltas)
      tmpdeltas$dv = dv
      #print(names(ppd)[dv])
      tmpdeltas$Female = slices[r,'Female']
      tmpdeltas$Black = slices[r,'Black']
      
      deltas = rbind(deltas,tmpdeltas)
      
      
      
          plts[,dv,2,r] = eff(tmp10)
          plts[,dv,1,r] = eff(tmp9) 
  }#end dv loop
}#end dem (r) loop

#print(deltas)
deltas$dv = factor(deltas$dv,labels=c('All Cause','Underlying Cause'))
deltas$Black = factor(deltas$Black,labels=c('White','Black'))
deltas$Female = factor(deltas$Female,labels=c('Male','Female'))


#prepare plots
library(ggplot2)

#deltas plot (1/24/2017)
#####THIS ONE
p = ggplot(deltas,aes(x=Female,y=mean))
  p + geom_point() + 
      geom_errorbar(aes(ymin=lower,ymax=upper),width=0.1) +
      facet_grid(.~Black+dv)

library(dplyr)
  
#omnibus -- calculate PROPORITON CHANGE BETWEEN ICD9 AND ICD 10 AS % DECLINE
  
wt.ac = ppd[['ac']] %>% 
    select(-Black, -Female) %>%  
    group_by(icd10) %>% 
    mutate(wt=tdeaths/sum(tdeaths)) %>%
    mutate_each(funs(.*wt),-tdeaths,-icd10,-wt) %>%
    ungroup


print(sum(wt.ac$wt)) #should equal 2

sum.ac = wt.ac %>% 
         group_by(icd10) %>% 
         summarise_each(funs(sum)) %>%
         select(-tdeaths, -wt) %>%
         ungroup

#tst=as.data.frame(t(sum.ac[,2:5001]))
#ggplot(exp(tst),aes(x=V1,y=V2)) + geom_point()


#bayesian pval
#sum(exp(sum.ac[2,2:5001])-exp(sum.ac[1,2:5001])>0)/5000


wt.uc = ppd[['uc']] %>% 
  select(-Black, -Female) %>%  
  group_by(icd10) %>% 
  mutate(wt=tdeaths/sum(tdeaths)) %>%
  mutate_each(funs(.*wt),-tdeaths,-icd10,-wt) %>%
  ungroup

print(sum(wt.uc$wt)) #should equal 2

sum.uc = wt.uc %>% 
  group_by(icd10) %>% 
  summarise_each(funs(sum)) %>%
  select(-tdeaths, -wt) %>%
  ungroup


delta.ac =as.numeric(exp(sum.ac[2,2:5001])-exp(sum.ac[1,2:5001])) 
scaled.ac=delta.ac/as.numeric(exp(sum.ac[2,2:5001]))
eff(scaled.ac) #proportion change ac

delta.uc=as.numeric(exp(sum.uc[2,2:5001])-exp(sum.uc[1,2:5001]))
scaled.uc=delta.uc/as.numeric(exp(sum.uc[2,2:5001]))
eff(scaled.uc)

#pval 
(sum(scaled.ac>scaled.uc)/5000)*2
all.scaled=data.frame(ac=scaled.ac,uc=scaled.uc)
all.delta=data.frame(ac=delta.ac,uc=delta.uc)
actual.ac=t(exp(sum.ac[,2:5001])); colnames(actual.ac) = c('ic9','icd10')
actual.uc=t(exp(sum.uc[,2:5001])); colnames(actual.uc) = c('ic9','icd10')

library(reshape2)

#View(melt(all.scaled))

#mention the point changes
ggplot(melt(all.scaled)) + geom_density(aes(x=value,fill=variable),alpha=.35)

#THIS ONE                      
ggplot(melt(all.delta)) + geom_density(aes(x=value,fill=variable),alpha=.35)


  #tmp9 = apply(ppd[[dv]][l9,paste0(1:5000)],2,FUN=function(x) sum(x*wt.9))  
  #tmp10 = apply(ppd[[dv]][l10,paste0(1:5000)],2,FUN=function(x) sum(x*wt.10))
  
  
  
par(mfrow=c(2,1))
#max=apply(plts,2,max)
#print(max)
plts=exp(plts)
rg=apply(plts,2,range)
print(rg);
#pltslog = plts

#plot saved manually (laptop)

par(mfrow=c(2,1),mar=c(3,7,1,2),cex=.9)
mp.1 = barplot(plts[1,1,,],beside=TRUE, horiz=TRUE
               , xlim=c(0,max(rg))
               ,xaxt='n'
               #, legend.text=c('ICD-9','ICD-10')
               , las=1,main="All Cause (RRA)")
  arrows(plts[2,1,,],mp.1,plts[3,1,,],mp.1, code=3, angle=90,length=.01)
  
mp.2 = barplot(plts[1,2,,],beside=TRUE, horiz=TRUE
              ,xlim=c(0,max(rg))
              ,legend.text=c('ICD-9','ICD-10')
              ,args.legend=c(x='bottomright',bty='n')
              ,las=1,main="Underlying Cause (RRD)")
  arrows( plts[2,2,,],mp.2,plts[3,2,,],mp.2, code=3, angle=90,length=.01)

#@@@@@@@@@@@@@@@@@@@
#Age Specific Plot
#@@@@@@@@@@@@@@@@@@@

#add complexity and age for predicted effect of complexity and age
ageplot = function(mod){
  #input approppiate model number; must be 2 or 4
    #pull effects and 84% intervals; include predicted effects fro complex = 1
    complex = model[[mod]]$gamma[,,12]
    #complex = matrix(0,2,2)
    yrrac9 = apply(model[[mod]]$gamma[,1,1:10],2,FUN=function(x) eff(exp(x), c=.84))
    yrrac10 = apply(model[[mod]]$gamma[,2,1:10],2,FUN=function(x) eff(exp(x), c=.84))

    #yl=range(c(yrrac9,yrrac10))
    #print(yl)
    yl=c(range(c(yrrac9,yrrac10)))

    plot(1,type='n',ylim=yl,xlim=c(1,10),log='y')
    o = 0 #offset for visibility

    #polygon(c((1:10)-o,rev((1:10)-0)),c(yrrac9[2,],rev(yrrac9[3,])),
    #        border=NA, col=gray(0.75,alpha=.25))
    #polygon(c((1:10)-o,rev((1:10)-0)),c(yrrac10[2,],rev(yrrac10[3,])),
    #        border=NA, col=gray(0.75,alpha=.25))

    lines((1:10)-o,yrrac9['mean',], type='p',pch=10)
    #lines(yrrac9[1,])
    arrows((1:10)-o,yrrac9[2,],(1:10)-o,yrrac9[3,],angle=90,code=3,length=.1)

    lines((1:10)+o,yrrac10['mean',], type='p',pch=16)
    #lines(yrrac10[1,],lty=2)
    arrows((1:10)+o,yrrac10[2,],(1:10)+o,yrrac10[3,],angle=90,code=3,length=.1)

}



#png(paste0(imdir,'ageplots.png'))

par(mfrow=c(2,1), oma=c(3,3,1,1), mar=c(0.5,1.5,0,0), font.main=1)
ageplot(2)
ageplot(4)

#dev.off()

d.all=apply(model[[2]]$gamma[,,1:10],3,delta)
d.under = apply(model[[4]]$gamma[,,1:10],3,delta)

xl = range(c(d.all,d.under))

#need lables for age
#png(paste0(imdir,'age_diff.png'))

par(mfrow=c(1,1))
plot(1,type='n',ylim=xl,xlim=c(1,10))

o=0.1 #offset for display

lines((1:10)-o,d.all[1,], type='p',pch=15)
  segments((1:10)-o,d.all[2,],(1:10)-o,d.all[3,])
  segments((1:10)+o,d.under[2,],(1:10)+o,d.under[3,])
lines((1:10)+o,d.under[1,], type='p',pch=16)
abline(h=0)
tick=axTicks(1) #get axis ticks from bottom side
#plot logged rate
axis(3,at=tick,labels=rnd(exp(tick),2))


#dev.off()

#@@@@@@@@@@@
#time series plot with ppd
#@@@@@@@@@@@

#calculate weight matrix by year

tdeaths = aggregate(dat$tdeaths,by=list(dat$Years),sum)

dat$wt=as.numeric(NA)
for(y in unique(dat$Years)){
  print(y)
  dat$wt[dat$Years==y] = dat$tdeaths[dat$Years==y]/tdeaths$x[tdeaths$Group.1 == y]
  }

#summarize posterior draws
#can take off exponentials, maybe

ppd.yrrac.exp = t(exp(model[[2]]$ppd))
ppd.yrrac = apply(ppd.yrrac.exp,2,FUN=function(x) x*dat$wt)
#ppd.yrrac = cbind(t(exp(model[[2]]$ppd)),dat$wt); wtcol = ncol(tst.yrrac);
#ppd.yrrac = apply(t(model[[2]]$ppd),2,FUN=function(x) x*dat$wt)
  tst.m.yrrac = aggregate(ppd.yrrac.exp,by=list(dat$Years),mean)
  plt.m.yrrac = apply(tst.m.yrrac[,2:ncol(ppd.yrrac)],1,eff)

  tst.yrrac = aggregate(list(ppd.yrrac,dat$wt),by=list(dat$Years),sum)
  plt.yrrac = apply(tst.yrrac[,2:ncol(ppd.yrrac)],1,eff)

lim = is.finite(dat$yrrdc)

ppd.yrrdc.exp = t(exp(model[[4]]$ppd))
ppd.yrrdc = apply(ppd.yrrdc.exp,2,FUN=function(x) x*dat$wt[lim])
#ppd.yrrdc = t(model[[4]]$ppd)
  tst.m.yrrdc = aggregate(ppd.yrrdc.exp,by=list(dat$Years[lim]),mean)
  plt.m.yrrdc = apply(tst.m.yrrdc[,2:ncol(tst.m.yrrdc)],1,eff)

  tst.yrrdc = aggregate(list(ppd.yrrdc,dat$wt[lim]),by=list(dat$Years[lim]),sum)
  plt.yrrdc = apply(tst.yrrdc[,2:ncol(ppd.yrrdc)],1,eff)


ob.m.yrrac=aggregate(exp(dat$yrrac),by=list(dat$Years),mean)
ob.m.yrrdc=aggregate(exp(dat$yrrdc[lim]),by=list(dat$Years[lim]),mean)

ob.yrrac=aggregate(exp(dat$yrrac)*dat$wt,by=list(dat$Years),sum)
ob.yrrdc=aggregate(exp(dat$yrrdc[lim])*dat$wt[lim],by=list(dat$Years[lim]),sum)
#ob.yrrac=aggregate(dat$yrrac,by=list(dat$Years),mean)
#ob.yrrdc=aggregate(dat$yrrdc[is.finite(dat$yrrdc)],by=list(dat$Years[is.finite(dat$yrrdc)]),mean)


#@@@@@@@@@@@@@@@@@@@@@@@
#unweighted ppd plot
#@@@@@@@@@@@@@@@@@@@@@@


png(paste0(imdir,'logscale-series1.png'))
par(mfrow=c(1,1),mar=c(3,3,3,3))


yl=range(c(ob.m.yrrac$x,plt.m.yrrac,ob.m.yrrdc$x,plt.m.yrrdc))
yl=range(c(ob.yrrac$x,plt.m.yrrac,ob.yrrdc$x,plt.yrrdc))

#yl=range(c(ob.yrrac$x,plt.yrrac))
plot(1,type='n',ylim=yl,xlim=c(1,10),xaxt='n',log="y")

polygon(c(1:5,rev(1:5)),
        c(plt.m.yrrac[2,1:5],rev(plt.m.yrrac[3,1:5])),
        border=NA, col=gray(0.9)
)


#lines(1:5,plt.yrrac[1,1:5],type="l")
#lines(1:5,ob.yrrac$x[1:5],type="p",pch=10)
lines(1:5,plt.m.yrrac[1,1:5],type="l", lty=2)
lines(1:5,ob.m.yrrac$x[1:5],type="p",pch=15)


#icd10
polygon(c(6:10,rev(6:10)),
        c(plt.m.yrrac[2,6:10],rev(plt.m.yrrac[3,6:10])),
        border=NA, col=gray(0.75)
)
lines(6:10,plt.m.yrrac[1,6:10],type="l",lty=2)
lines(6:10,ob.m.yrrac$x[6:10],type="p",pch=15)

abline(v=5.5)

polygon(c(1:5,rev(1:5)),
        c(plt.m.yrrdc[2,1:5],rev(plt.m.yrrdc[3,1:5])),
        border=NA, col=gray(0.75)
)
lines(1:5,plt.m.yrrdc[1,1:5],type="l",lty=2)
lines(1:5,ob.m.yrrdc$x[1:5],type="p",pch=16)

#icd10
polygon(c(6:10,rev(6:10)),
        c(plt.m.yrrdc[2,6:10],rev(plt.m.yrrdc[3,6:10])),
        border=NA, col=gray(0.9)
)
lines(6:10,plt.m.yrrdc[1,6:10],type="l",lty=2)
lines(6:10,ob.m.yrrdc$x[6:10],type="p",pch=16)

axis(1,at=1:10,labels=1994:2003)
T=axTicks(2) #get axis ticks from left side
#plot logged rate
axis(4,at=T,labels=rnd(log(T),2))
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@
#weighted ppd plot
#@@@@@@@@@@@@@@@@@@@@@@


png(paste0(imdir,'logscale-series1-wt.png'))
par(mfrow=c(1,1),mar=c(3,3,3,3))

yl=range(c(ob.yrrac$x,plt.m.yrrac,ob.yrrdc$x,plt.m.yrrdc))

#yl=range(c(ob.yrrac$x,plt.yrrac))
plot(1,type='n',ylim=yl,xlim=c(1,10),xaxt='n',log="y")

polygon(c(1:5,rev(1:5)),
        c(plt.m.yrrac[2,1:5],rev(plt.m.yrrac[3,1:5])),
        border=NA, col=gray(0.9)
)


#lines(1:5,plt.yrrac[1,1:5],type="l")
#lines(1:5,ob.yrrac$x[1:5],type="p",pch=10)
lines(1:5,plt.m.yrrac[1,1:5],type="l", lty=2)
lines(1:5,ob.yrrac$x[1:5],type="p",pch=15)


#icd10
polygon(c(6:10,rev(6:10)),
        c(plt.yrrac[2,6:10],rev(plt.yrrac[3,6:10])),
        border=NA, col=gray(0.75)
)
lines(6:10,plt.yrrac[1,6:10],type="l",lty=2)
lines(6:10,ob.yrrac$x[6:10],type="p",pch=15)

abline(v=5.5)

polygon(c(1:5,rev(1:5)),
        c(plt.yrrdc[2,1:5],rev(plt.yrrdc[3,1:5])),
        border=NA, col=gray(0.75)
)
lines(1:5,plt.yrrdc[1,1:5],type="l",lty=2)
lines(1:5,ob.yrrdc$x[1:5],type="p",pch=16)

#icd10
polygon(c(6:10,rev(6:10)),
        c(plt.yrrdc[2,6:10],rev(plt.yrrdc[3,6:10])),
        border=NA, col=gray(0.75)
)
lines(6:10,plt.yrrdc[1,6:10],type="l",lty=2)
lines(6:10,ob.yrrdc$x[6:10],type="p",pch=16)

axis(1,at=1:10,labels=1994:2003)
T=axTicks(2) #get axis ticks from left side
#plot logged rate
axis(4,at=T,labels=rnd(log(T),2))

dev.off()


#distrubution of cell effects
#mu_i = model[[2]]$mu_i









###############
###############
#delta ppd draw (across two levels) - BDA III, p. 118
###############
###############

#load data as saved to
load(paste0(outdir,'yrracdat.RData'))
load(paste0(outdir,'yrrdcdat.RData'))

#reassign "t" to dat$year -- it gets lost for some reason (probably dynamic assignment in code)
#yrracdat$t = dat$Years
#yrrdcdat$t = dat$Years[is.finite(dat$yrrdc)]

#copy structure -- need to draw a random year variable,
#and extend through the entire set of periods (based on a1 and a2)

makeppt=function(m,dat){
  #m is model number - 2 or 4
  # dat is a design matrix of data in a list in the same form as the stan model
  # dat is saved from the analysis set up
  #returns a list of ppts including implied icd9 and icd10

  #make container the size of ppt without periods
  #(to cover both icd9 and icd10)
  ppt = list(model[[m]]$ppd,model[[m]]$ppd)
  names(ppt) = c('icd9','icd10')

  #iterate through posterior
  for(iter in 1:nrow(ppt[[1]])){

    if(iter%%500==0){cat('Coding',iter,'of',nrow(ppt[[1]]),'\n')}

    #collect previously drawn mu_i's
    mu_i.9=cbind(1:dat$IDS,model[[m]]$mu_i[iter,1,])
    mu_i.10=cbind(1:dat$IDS,model[[m]]$mu_i[iter,2,])
    ieff.9=ieff.10=cbind(dat$id,NA)
    for(i in 1:nrow(mu_i.9)){
      ieff.9[ieff.9[,1]==i,2] = mu_i.9[i,2]
      ieff.10[ieff.10[,1]==i,2] = mu_i.10[i,2]
    }

    #calculate E(y) (using only gammas and alpha(drawn from mu))
    yhat.9 = dat$z %*% model[[m]]$gamma[iter,1,] + ieff.9[,2];
    yhat.10 = dat$z %*% model[[m]]$gamma[iter,2,] + ieff.10[,2];

    if(iter%%1000==0){
      cat('yhat (mean,min,max)\n')
      print(c(mean(yhat.9),range(yhat.9)))
      print(c(mean(yhat.10),range(yhat.10)))
    }

    ##draw tilde{delta}, tilde{y}
    yr=range(dat$t)
    mu_t.9 = cbind(yr[1]:yr[2],rscaled_t(n=10,df=9,s=model[[m]]$delta[iter]))
    mu_t.10 = cbind(yr[1]:yr[2],rscaled_t(n=10,df=9,s=model[[m]]$delta[iter]))
    teff.9 = teff.10 = cbind(dat$t,NA)
    for(t in 1:nrow(mu_t.9)){
      teff.9[teff.9[,1]==mu_t.9[t,1],2] = mu_t.9[t,2]
      teff.10[teff.10[,1]==mu_t.10[t,1],2] = mu_t.10[t,2]
    }

    newhat.9 = yhat.9+model[[m]]$beta[iter,1]*dat$t+teff.9[,2]*dat$t
    newhat.10 = yhat.10+model[[m]]$beta[iter,2]*dat$t+teff.10[,2]*dat$t

    if(iter%%1000==0){
      cat('new yhat (mean,min,max)\n')
      print(c(mean(newhat.9),range(newhat.9)))
      print(c(mean(newhat.10),range(newhat.10)))
    }

    n=ncol(ppt[[1]])
    ppt[[1]][iter,] = rnorm(n,mean=newhat.9,sd=model[[m]]$sig[iter,1])
    ppt[[2]][iter,] = rnorm(n,mean=newhat.10,sd=model[[m]]$sig[iter,2])

    if(iter%%1000==0){
      cat('tilde{y} (mean,min,max)\n')
      print(c(mean(ppt[[1]]),range(ppt[[1]])))
      print(c(mean(ppt[[2]]),range(ppt[[2]])))
    }


  }

  return(ppt)

}

ppt.yrrac = makeppt(m=2,dat=yrracdat)
ppt.yrrdc = makeppt(m=4,dat=yrrdcdat)

#generate exponentiated ppdt

yrrac.exp = lapply(ppt.yrrac,FUN=function(x) t(exp(x)))
yrrac.m = lapply(yrrac.exp, FUN= function(x)
          aggregate(x,by=list(yrracdat$t),mean))
yrrac.plt = lapply(yrrac.m,FUN=function(x)
          apply(x[,2:ncol(x)],1,eff))

yrrdc.exp = lapply(ppt.yrrdc,FUN=function(x) t(exp(x)))
yrrdc.m = lapply(yrrdc.exp, FUN= function(x)
  aggregate(x,by=list(yrrdcdat$t),mean))
yrrdc.plt = lapply(yrrdc.m,FUN=function(x)
  apply(x[,2:ncol(x)],1,eff))
yrrdc.med=lapply(yrrdc.m,FUN=function(x) 
  apply(x[,2:ncol(x)],1,quantile, prob=.5))

#@@@@@@@@@@@@@@@@@@@@@@@
#unweighted ppd plot
#@@@@@@@@@@@@@@@@@@@@@@

png(paste0(imdir,'mean-rel-rates.png'),height=4.5,width=6.5,res=500,units='in')

#pdf(paste0(imdir,'series2.pdf'))

gr = grey.colors(2,start=0.4,end=0.6,alpha=0.3)
grdot = grey.colors(1,start=0.6,end=0.8)

par(mfrow=c(1,1),mar=c(2,1,1,2),cex=.75)


yl=range(c(ob.m.yrrac$x,yrrac.plt,yrrdc.plt,ob.m.yrrdc$x))
#yl=range(c(yrrdc.plt,ob.m.yrrdc$x))

plot(1,type='n',ylim=yl,xlim=c(0.5,10),xaxt='n',axes=FALSE)

polygon(c(1:10,rev(1:10)),
        c(yrrac.plt$icd9[2,1:10],rev(yrrac.plt$icd9[3,1:10])),
        border=gr[2] ,col=gr[1]
)

polygon(c(1:10,rev(1:10)),
        c(yrrac.plt$icd10[2,1:10],rev(yrrac.plt$icd10[3,1:10])),
         border=gr[2],col=gr[2]
)


polygon(c(1:10,rev(1:10)),
        c(yrrdc.plt$icd9[2,1:10],rev(yrrdc.plt$icd9[3,1:10])),
        border=gr[1], col=gr[1]
)


polygon(c(1:10,rev(1:10)),
        c(yrrdc.plt$icd10[2,1:10],rev(yrrdc.plt$icd10[3,1:10])),
        border=gr[2], col=gr[2]
)


#lines(1:10,plt.m.yrrac[1,1:10],type="l", lty=1)
lines(1:5,yrrac.plt$icd9[1,1:5],type='l',lty=1)
lines(5:10,yrrac.plt$icd9[1,5:10],type='l',lty=2)

lines(6:10,yrrac.plt$icd10[1,6:10],type='l',lty=1)
lines(1:6,yrrac.plt$icd10[1,1:6],type='l',lty=2)


lines(1:5,yrrdc.med$icd9[1:5],type='l',lty=1)
lines(5:10,yrrdc.med$icd9[5:10],type='l',lty=2)

lines(6:10,yrrdc.med$icd10[6:10],type='l',lty=1)
lines(1:6,yrrdc.med$icd10[1:6],type='l',lty=2)

#need to switch-up colors??
lines(1:10,ob.m.yrrac$x[1:10],type="p",pch=16,col=grdot)
lines(1:10,ob.m.yrrdc$x[1:10],type='p',pch=16,col=grdot)

abline(v=5.5)
axis(4)
axis(1,at=1:10,labels=1994:2003) 

#tick=axTicks(2) #get axis ticks from bottom side
#plot logged rate
#axis(4,at=tick,labels=rnd(log(tick),2))

mtext('Logged Relative Rate',side=4,padj=3,cex=.7)

text(0.5,.3,'All-Cause (RRA)',srt=90)
text(0.75,yrrac.plt$icd9[1,1],'ICD-9', srt=90)
text(0.75,yrrac.plt$icd10[1,1],'ICD-10', srt=90)

text(0.5,.12,'Underlying-Cause (RRD)',srt=90)
text(0.75,yrrdc.plt$icd9[1,1],'ICD-9', srt=90)
text(0.75,yrrdc.plt$icd10[1,1],'ICD-10', srt=90)


text(5.5,.36,'  ICD-10\n  Implemented',adj=0)

legend('topright',c('Observed','Estimated','Projected'),
       pch=c(16,NA,NA),
       lty=c(NA,1,2),
       col=c(grdot,'black','black'),
       bty='n')

dev.off() 


###########
###########
# PPD age-group plots
###########
###########

png(paste0(imdir,'ageplot-99-ppd-84ci.png'),height=5.5,width=5.5,res=500,units='in')

yrrac.agelist=rep(0,nrow(yrracdat$z))
yrrdc.agelist=rep(0,nrow(yrrdcdat$z))
for(a in 1:10){
  yrrac.agelist[yrracdat$z[,a]==1] = a
  yrrdc.agelist[yrrdcdat$z[,a]==1] = a
  }

#collect mean rates by age groups for 1999
###YRRAC
#par(mfrow=c(5,2))
#for(year in unique(yrracdat$t)){
par(mfrow=c(2,1), oma=c(3,2,1,1), mar=c(0.5,1.5,3,0), font.main=1)
year = 0
yrrac.agep= lapply(yrrac.exp,FUN=function(x)
         aggregate(x[yrracdat$t==year,],by=list(yrrac.agelist[yrracdat$t==year]),mean))

yaplt = lapply(yrrac.agep,FUN=function(x) apply(x,1,eff,c=.84))

plot(1:10,yaplt$icd9[1,],ylim=range(yaplt),type='p',pch=16,axes=FALSE,log='y')
  arrows(1:10,yaplt$icd9[2,],1:10,yaplt$icd9[3,],angle=90,length=.1,code=3)
lines(1:10,yaplt$icd10[1,],ylim=range(yaplt),type='p', pch=1)
  arrows(1:10,yaplt$icd10[2,],1:10,yaplt$icd10[3,],angle=90,length=.1,code=3)

  legend('topright',c('ICD9','ICD10'), pch=c(16,1), bty='n', cex=0.8)
  mtext("All-Cause (RRA)", side=3, outer=F, line=0)
  axis(2)

 #} crosses over in 1999

###YRRDC
  #par(mfrow=c(5,2))
  
#for(year in unique(yrrdcdat$t)){
year = 0  
  
yrrdc.agep= lapply(yrrdc.exp,FUN=function(x)
    aggregate(x[yrrdcdat$t==year,],by=list(yrrdc.agelist[yrrdcdat$t==year]),mean))
  
  ydplt = lapply(yrrdc.agep,FUN=function(x) apply(x,1,eff,c=.84))
  
  plot(1:10,ydplt$icd9[1,],ylim=range(ydplt),type='p', pch=16,log='y',axes=FALSE)
  arrows(1:10,ydplt$icd9[2,],1:10,ydplt$icd9[3,],angle=90,length=.1,code=3)
  lines(1:10,ydplt$icd10[1,],ylim=range(ydplt),type='p', pch=1)
  arrows(1:10,ydplt$icd10[2,],1:10,ydplt$icd10[3,],angle=90,length=.1,code=3)
  
  mtext("Underlying Cause (RRD)", side=3, outer=F, line=0)
  axis(1,at=1:10, labels=c('40-44','45-49','50-54','55-59','60-64','65-69','70-74','75-79','80-84','85+'))
  axis(2)
  #legend('topright',legend=year+1999,bty='n')
#}crosses over at 1997-2000
  
dev.off()

#@@@@@@@@@@@@@@@@@@@@@@@
#unweighted spagheti plot
#@@@@@@@@@@@@@@@@@@@@@@

smp = lapply(yrrdc.m,FUN=function(x) sample(x[,2:ncol(x)],25))

yl=range(c(smp,ob.m.yrrdc$x))
plot(1,type='n',ylim=yl,xlim=c(1,10),xaxt='n')#,log='y')
  apply(smp$icd9,2,FUN=function(x)
      lines(1:5,x[1:5], col=gray(0.75))
      #lines(5:10,x[5:10],col=gray(0.75),lty=2)
    )
  apply(smp$icd10,2,FUN=function(x)
     lines(1:10,x,col=gray(0.75)))
  #lines(1:10,yrrdc.smp$icd9[,1],col=gray(.9))
  lines(1:10,ob.m.yrrdc$x[1:10],type='p',pch=16)

  
  #@@@@@@@@@@@@@@@@@@@
  #Omnibus test of difference [difference in proportions?]
  #@@@@@@@@@@@@@@@@@@@


sink(paste0(outdir,'icd9-diff-pval.txt'))
cat('###File collects info on bayesian p-values\n')  
    
#calculates difference between mean rates under ICD9 and ICD10
tst=function(data){
  mnrt = aggregate(data,by=list(yrracdat$td),mean)
    return(mnrt$x[2]-mnrt$x[1])
  }

cat('\n\n@@@@@@@@@@@\nYRRAC\n@@@@@@@@@@@@@@@@\n\n')  
mn.yrrac=tst(yrracdat$y)
mn.ppdac = apply(model[[2]]$ppd,1,tst)
#bayesian p-values
cat('\nObserved difference in mean rates:\t\t',mn.yrrac)
cat('\nPPD difference (over y):          \t\t')
  printeff(mn.ppdac)
cat('\nBayesian p-val for mean value YRRAC < 0:\t\t',
        sum(mn.ppdac<0)/length(mn.ppdac)) 
cat('\nBayesian p-val for mean value less than observed:\t',
        sum(mn.ppdac<mn.yrrac)/length(mn.ppdac))

#test including undertainty over years
#prepare posterior draws WITHOUT projection and antejection
mn.pptac = ppt.yrrac$icd9
mn.pptac[,yrracdat$td==2] = ppt.yrrac$icd9[,yrracdat$td==2]
mn.pptac=apply(mn.pptac,1,tst)
#bayesian p-values
cat('\n\nPPD difference (over y,tau):\t')
  printeff(mn.pptac)
cat('\nBayesian p-val for mean value YRRAC < 0:\t\t',
    sum(mn.pptac<0)/length(mn.pptac)) 
cat('\nBayesian p-val for mean value less than observed:\t',
    sum(mn.pptac<mn.yrrac)/length(mn.pptac))


cat('\n\n@@@@@@@@@@@\nYRRDC\n@@@@@@@@@@@@@@@@\n\n') 

tst=function(data){
  mnrt = aggregate(data,by=list(yrrdcdat$td),mean)
  return(mnrt$x[2]-mnrt$x[1])
}

mn.yrrdc=tst(yrrdcdat$y)  
mn.ppddc = apply(model[[4]]$ppd,1,tst)
#bayesian p-value
cat('\nObserved difference in mean rates:\t\t',mn.yrrdc)
cat('\nPPD difference (over y):          \t\t')
printeff(mn.ppddc)
cat('\nBayesian p-val for mean value YRRDC < 0:\t\t',
    sum(mn.ppddc<0)/length(mn.ppddc)) 
cat('\nBayesian p-val for mean value less than observed:\t',
    sum(mn.ppddc<mn.yrrac)/length(mn.ppddc))

mn.pptdc = ppt.yrrdc$icd9
mn.pptdc[,yrrdcdat$td==2] = ppt.yrrdc$icd9[,yrrdcdat$td==2]
mn.pptdc=apply(mn.pptdc,1,tst)
#bayesian p-values
cat('\n\nPPD difference (over y,tau):\t')
printeff(mn.pptdc)
cat('\nBayesian p-val for mean value YRRDC < 0:\t\t',
    sum(mn.pptdc<0)/length(mn.pptdc)) 
cat('\nBayesian p-val for mean value less than observed:\t',
    sum(mn.pptdc<mn.yrrac)/length(mn.pptdc))
sink()
  
########
#age standardized plot---
#no good (need individual year changes)
#wrong weights -- would need population
########

#age standardized by weighting according to 1999
wt = dat$tdeaths[dat$Years==0]/sum(dat$tdeaths[dat$Years==0])
wtlim = dat$tdeaths[lim & dat$Years==0]/sum(dat$tdeaths[lim & dat$Years==0])

#calculate weighted mean
yrrac.wtm = lapply(yrrac.exp, FUN= function(x)
  aggregate(x,by=list(dat$Years),FUN=function(y) sum(y*wt)))
yrrac.wtplt = lapply(yrrac.wtm,FUN=function(x)
  apply(x[,2:ncol(x)],1,eff))

yrrdc.wtm = lapply(yrrdc.exp,FUN=function(x)
  aggregate(x,by=list(yrrdcdat$t),FUN=function(y) sum(y*wt)))
yrrdc.wtplt = lapply(yrrdc.wtm,FUN=function(x)
  apply(x[,2:ncol(x)],1,eff))

ob.wtyrrac=aggregate(exp(dat$yrrac),by=list(dat$Years),FUN=
          function(x) sum(x*wt))

ob.wtyrrdc=aggregate(exp(dat$yrrdc),by=list(dat$Years),FUN=
                       function(x) sum(x*wt))

png(paste0(imdir,'wt-series2.png'))
par(mfrow=c(1,1),mar=c(3,3,3,3))


yl=range(c(ob.wtyrrac$x,yrrac.wtplt,ob.wtyrrdc$x,yrrdc.wtplt))

#yl=range(c(ob.yrrac$x,plt.yrrac))
plot(1,type='n',ylim=yl,xlim=c(1,10),xaxt='n')#,log='y')

polygon(c(1:10,rev(1:10)),
        c(yrrac.wtplt$icd9[2,1:10],rev(yrrac.wtplt$icd9[3,1:10])),
        border=NA, col=gray(0.75,alpha=.25)
)

polygon(c(1:10,rev(1:10)),
        c(yrrac.wtplt$icd10[2,1:10],rev(yrrac.wtplt$icd10[3,1:10])),
        border=NA, col=gray(0.75,alpha=.25)
)


polygon(c(1:10,rev(1:10)),
        c(yrrdc.wtplt$icd9[2,1:10],rev(yrrdc.wtplt$icd9[3,1:10])),
        border=NA, col=gray(0.75,alpha=.25)
)

polygon(c(1:10,rev(1:10)),
        c(yrrdc.wtplt$icd10[2,1:10],rev(yrrdc.wtplt$icd10[3,1:10])),
        border=NA, col=gray(0.75,alpha=.25)
)


#lines(1:10,plt.m.yrrac[1,1:10],type="l", lty=1)
lines(1:10,yrrac.wtplt$icd9[1,],type='l',lty=2)
lines(1:10,yrrac.wtplt$icd10[1,],type='l',lty=3)

lines(1:10,yrrdc.wtplt$icd9[1,],type='l',lty=2)
lines(1:10,yrrdc.wtplt$icd10[1,],type='l',lty=3)


lines(1:10,ob.wtyrrac$x[1:10],type="p",pch=15)
lines(1:10,ob.wtyrrdc$x[1:10],type="p",pch=16)

dev.off()


######
#complexity 
#need to fix where it is centered
######

cmplx = matrix(0,4,5)
colnames(cmplx) =  c('complex','icd10','mean','ll','ul')
cmplx[,1] = c(1,1,0,0)
cmplx[,2]=c(0,1,0,1)

cmplx.ppd = list(cmplx,cmplx); names(cmplx.ppd)=c('ac','uc')

#Calculate stuff all-cause
cmplx.ppd$ac[1,3:5] = eff(ppd.yrrac.exp[dat$Complex==1 & dat$Years>=0,])
eff(ppd.yrrac.exp[dat$Complex==1 & dat$Years<0,])

eff(ppd.yrrac.exp[dat$Complex==0 & dat$Years>=0,])
eff(ppd.yrrac.exp[dat$Complex==0 & dat$Years<0,])

#underlying cause
eff(ppd.yrrdc.exp[datuc$Complex==1 & datuc$Years>=0,])
eff(ppd.yrrdc.exp[datuc$Complex==1 & datuc$Years<0,])

eff(ppd.yrrdc.exp[datuc$Complex==0 & datuc$Years>=0,])
eff(ppd.yrrdc.exp[datuc$Complex==0 & datuc$Years<0,])


