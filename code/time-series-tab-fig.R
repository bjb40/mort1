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
  
  d = apply(s,1,FUN=function(x) x[2]-x[1])
  if(missing(p)){p=F}
  if(p==F){return(eff(d,c))}
  if(p==T){return(printeff(d,c))}
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
cat(delta(model[[2]]$beta,p=T),'|\n')

#conditional means of effects - gammas
for(e in 1:16){
  cat('|',gammanames[e],'|')
  cat(printeff(model[[1]]$gamma[,e]), '|')
  #model 2 is three dimensional
  cat(printeff(model[[2]]$gamma[,1,e]), '|')
  cat(printeff(model[[2]]$gamma[,2,e]), '|')
  cat(delta(model[[2]]$gamma[,,e],p=T),'|\n')
  
}

#variances

####need to confirm terminology -- am I getting standard deviations or variances?

cat('| Level 2 Var |')
cat(printeff(model[[1]]$zi),'|')
cat(printeff(model[[2]]$zi[,1]), '|')
cat(printeff(model[[2]]$zi[,2]), '|')
cat(delta(model[[2]]$zi,p=T),'|\n')

cat('| Correlation of Level 2 Var |   -   |    |')
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
cat(delta(model[[4]]$beta,p=T),'|\n')

#conditional means of effects - gammas
for(e in 1:16){
  cat('|',gammanames[e],'|')
  cat(printeff(model[[3]]$gamma[,e]), '|')
  #model 2 is three dimensional
  cat(printeff(model[[4]]$gamma[,1,e]), '|')
  cat(printeff(model[[4]]$gamma[,2,e]), '|')
  cat(delta(model[[4]]$gamma[,,e],p=T),'|\n')
  
}

#variances

####need to confirm terminology -- am I getting standard deviations or variances?

cat('| Level 2 Var |')
cat(printeff(model[[3]]$zi),'|')
cat(printeff(model[[4]]$zi[,1]), '|')
cat(printeff(model[[4]]$zi[,2]), '|')
cat(delta(model[[4]]$zi,p=T),'|\n')

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
    yl=c(0.05,1)
    
    plot(1,type='n',ylim=yl,xlim=c(1,10),log='y')
    o = 0 #offset for visibility
    
    polygon(c((1:10)-o,rev((1:10)-0)),c(yrrac9[2,],rev(yrrac9[3,])),
            border=NA, col=gray(0.9))
    polygon(c((1:10)-o,rev((1:10)-0)),c(yrrac10[2,],rev(yrrac10[3,])),
            border=NA, col=gray(0.9))
    
    lines((1:10)-o,yrrac9['mean',], type='p',pch=10)
    lines(yrrac9[1,])
    #arrows((1:10)-o,yrrac9[2,],(1:10)-o,yrrac9[3,],angle=90,code=3,length=.1)
    
    lines((1:10)+o,yrrac10['mean',], type='p',pch=16)
    lines(yrrac10[1,],lty=2)
    #arrows((1:10)+o,yrrac10[2,],(1:10)+o,yrrac10[3,],angle=90,code=3,length=.1)

}

png(paste0(imdir,'ageplots.png'))

par(mfrow=c(2,1), oma=c(3,3,1,1), mar=c(0.5,1.5,0,0), font.main=1)
ageplot(2)
ageplot(4)
dev.off()

d.all=apply(model[[2]]$gamma[,,1:10],3,delta)
d.under = apply(model[[4]]$gamma[,,1:10],3,delta)

xl = range(c(d.all,d.under))

#need lables for age
png(paste0(imdir,'age_diff.png'))

par(mfrow=c(1,1))
plot(1,type='n',xlim=xl,ylim=c(1,10))

o=0.1 #offset for display

lines(d.all[1,],(1:10)-o, type='p',pch=15)
segments(d.all[2,],(1:10)-o,d.all[3,],(1:10)-o)
segments(d.under[2,],(1:10)+o,d.under[3,],(1:10)+o)
lines(d.under[1,],(1:10)+o, type='p',pch=16)
abline(v=0)
tick=axTicks(1) #get axis ticks from bottom side
#plot logged rate
axis(3,at=tick,labels=rnd(exp(tick),2))


dev.off()

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
        border=NA, col=gray(0.9)
)
lines(6:10,plt.m.yrrac[1,6:10],type="l",lty=2)
lines(6:10,ob.m.yrrac$x[6:10],type="p",pch=15)

abline(v=5.5)

polygon(c(1:5,rev(1:5)),
        c(plt.m.yrrdc[2,1:5],rev(plt.m.yrrdc[3,1:5])),
        border=NA, col=gray(0.9)
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
        c(plt.yrrac[2,1:5],rev(plt.yrrac[3,1:5])),
        border=NA, col=gray(0.9)
)


#lines(1:5,plt.yrrac[1,1:5],type="l")
#lines(1:5,ob.yrrac$x[1:5],type="p",pch=10)
lines(1:5,plt.m.yrrac[1,1:5],type="l", lty=2)
lines(1:5,ob.yrrac$x[1:5],type="p",pch=15)


#icd10
polygon(c(6:10,rev(6:10)),
        c(plt.yrrac[2,6:10],rev(plt.yrrac[3,6:10])),
        border=NA, col=gray(0.9)
)
lines(6:10,plt.yrrac[1,6:10],type="l",lty=2)
lines(6:10,ob.yrrac$x[6:10],type="p",pch=15)

abline(v=5.5)

polygon(c(1:5,rev(1:5)),
        c(plt.yrrdc[2,1:5],rev(plt.yrrdc[3,1:5])),
        border=NA, col=gray(0.9)
)
lines(1:5,plt.yrrdc[1,1:5],type="l",lty=2)
lines(1:5,ob.yrrdc$x[1:5],type="p",pch=16)

#icd10
polygon(c(6:10,rev(6:10)),
        c(plt.yrrdc[2,6:10],rev(plt.yrrdc[3,6:10])),
        border=NA, col=gray(0.9)
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
load(paste0(outdir,'bayesdat.RData'))
bayesdat$t = dat$Years

#copy structure -- need to draw a random year variable,
#and extend through the entire set of periods (based on a1 and a2)

makeppt=function(m){
  #m is model number - 2 or 4
  #returns a list of ppts including implied icd9 and icd10
  
  #make container the size of ppt without periods 
  #(to cover both icd9 and icd10)
  ppt = list(model[[m]]$ppd,model[[m]]$ppd)
  names(ppt) = c('icd9','icd10')
  
  #iterate through posterior
  for(iter in 1:nrow(ppt[[1]])){
    
    if(iter%%500==0){cat('Coding',iter,'of',nrow(ppt[[1]]),'\n')}
    
    #collect previously drawn mu_i's
    mu_i.9=cbind(1:bayesdat$IDS,model[[m]]$mu_i[iter,1,])
    mu_i.10=cbind(1:bayesdat$IDS,model[[m]]$mu_i[iter,2,])
    ieff.9=ieff.10=cbind(bayesdat$id,NA)
    for(i in 1:nrow(mu_i.9)){
      ieff.9[ieff.9[,1]==i,2] = mu_i.9[i,2]
      ieff.10[ieff.10[,1]==i,2] = mu_i.10[i,2]
    }
    
    #calculate E(y) (using only gammas and alpha(drawn from mu))
    yhat.9 = bayesdat$z %*% model[[m]]$gamma[iter,1,] + ieff.9[,2];
    yhat.10 = bayesdat$z %*% model[[m]]$gamma[iter,2,] + ieff.10[,2];
    
    if(iter%%1000==0){
      cat('yhat (mean,min,max)\n')
      print(c(mean(yhat.9),range(yhat.9)))
      print(c(mean(yhat.10),range(yhat.10)))
    }
    
    ##draw tilde{delta}, tilde{y}
    yr=range(bayesdat$t)
    mu_t.9 = cbind(yr[1]:yr[2],rnorm(10,mean=0,sd=model[[m]]$delta[iter]))
    mu_t.10 = cbind(yr[1]:yr[2],rnorm(10,mean=0,sd=model[[m]]$delta[iter]))
    teff.9 = teff.10 = cbind(bayesdat$t,NA)
    for(t in 1:nrow(mu_t.9)){
      teff.9[teff.9[,1]==mu_t.9[t,1],2] = mu_t.9[t,2]
      teff.10[teff.10[,1]==mu_t.10[t,1],2] = mu_t.10[t,2]
    }

    newhat.9 = yhat.9+model[[m]]$beta[iter,1]*bayesdat$t+teff.9[,2]*bayesdat$t
    newhat.10 = yhat.10+model[[m]]$beta[iter,2]*bayesdat$t+teff.10[,2]*bayesdat$t

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

ppt.yrrac = makeppt(2)
#ppt.yrrdc = makeppt(4)

#generate exponentiated ppdt

yrrac.exp = lapply(ppt.yrrac,FUN=function(x) t(exp(x)))
yrrac.log = lapply(ppt.yrrac,FUN=function(x) t(x))
#yrrac.wt = apply(ppd.yrrac.exp,2,FUN=function(x) x*dat$wt)
yrrac.m = lapply(yrrac.log, FUN= function(x)
          aggregate(x,by=list(dat$Years),mean))
yrrac.plt = lapply(yrrac.m,FUN=function(x)
          apply(x[,2:ncol(x)],1,eff))

ob.yrrac.log=aggregate(dat$yrrac,by=list(dat$Years),mean)


#@@@@@@@@@@@@@@@@@@@@@@@
#unweighted ppd plot
#@@@@@@@@@@@@@@@@@@@@@@

#png(paste0(imdir,'logscale-series1.png'))  
par(mfrow=c(1,1),mar=c(3,3,3,3))


yl=range(c(ob.yrrac.log$x,yrrac.plt))

#yl=range(c(ob.yrrac$x,plt.yrrac))
plot(1,type='n',ylim=yl,xlim=c(1,10),xaxt='n')

#polygon(c(1:10,rev(1:10)),
#        c(plt.m.yrrac[2,1:10],rev(plt.m.yrrac[3,1:10])),
#        border=NA, col=gray(0.9)
#)

polygon(c(1:10,rev(1:10)),
        c(yrrac.plt$icd9[2,1:10],rev(yrrac.plt$icd9[3,1:10])),
        border=NA, col=gray(0.9)
)

polygon(c(1:10,rev(1:10)),
        c(yrrac.plt$icd10[2,1:10],rev(yrrac.plt$icd10[3,1:10])),
        border=NA, col=gray(0.9)
)


#lines(1:10,plt.m.yrrac[1,1:10],type="l", lty=1)
lines(1:10,yrrac.plt$icd9[1,],type='l',lty=2)
lines(1:10,yrrac.plt$icd10[1,],type='l',lty=3)

lines(1:10,ob.yrrac.log$x[1:10],type="p",pch=15)

