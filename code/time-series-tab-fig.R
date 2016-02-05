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

plt = apply(model[[4]]$gamma[,1,1:10],2,eff, c=.84)
ext = model[[2]]$gamma[,,12]

plt = apply(model[[2]]$gamma[,,1:10],3,FUN=function(x) x + model[[2]]$gamma[,,12]) 

par(mfrow=c(2,1), oma=c(3,3,1,1), mar=c(0.5,1.5,0,0), font.main=1)

plot(1,type='n',ylim=c(-2,-0.6),xlim=c(1,10))

plot(1,type='n',ylim=c(0,.5),xlim=c(1,10))

#pull effects and 84% intervals
yrrac9 = apply(model[[4]]$gamma[,1,1:10],2,eff, c=.84)
yrrac10 = apply(model[[4]]$gamma[,2,1:10],2,eff, c=.84)

lines(yrrac9['mean',], type='p',pch=10)
lines(yrrac9[1,])
arrows(1:10,yrrac9[2,],1:10,yrrac9[3,],angle=90,code=3,length=.1)
lines(yrrac10['mean',], type='p',pch=16)
lines(yrrac10[1,],lty=2)
arrows(1:10,yrrac10[2,],1:10,yrrac10[3,],angle=90,code=3,length=.1)



ageplot = function(posterior){
  
  rpost = posterior$betas
  rpost[,22:30] = rpost[,22:30] + rpost[,11] + rpost[,2:10] #+ 4*post[,12]
  # adjust for year effect: only for year 1998 (prior to change)
  rpost[,2:10] = rpost[,2:10] + -1*rpost[,12] 
  post = apply(rpost,2,quantile, prob=c(0.08,0.5,0.82))
  post_mean = apply(rpost,2,mean)
  rm(rpost)
  
  plot(1:9+.5,post_mean[2:10], pch=10,
       ylim=c(min(post[,c(2:10,22:30)]), max(post[,c(2:10,22:30)])), 
       xlim=c(1,10), 
       xaxt="n")
  segments(1:9+.5,post[1,2:10],1:9+.5,post[3,2:10])
  segments(1:9+.4,post[1,2:10],2:10-.4,post[1,2:10]) 
  segments(1:9+.4,post[1,2:10],2:10-.4,post[1,2:10])
  segments(1:9+.4,post[3,2:10],2:10-.4,post[3,2:10]) 
  segments(1:9+.4,post[3,2:10],2:10-.4,post[3,2:10])
  
  lines(1:9+.5,post_mean[22:30], type='p', pch=16)
  segments(1:9+.5,post[1,22:30],1:9+.5,post[3,22:30], lty=2)
  segments(1:9+.4,post[1,22:30],2:10-.4,post[1,22:30], lty=2)
  segments(1:9+.4,post[1,22:30],2:10-.4,post[1,22:30], lty=2)
  segments(1:9+.4,post[3,22:30],2:10-.4,post[3,22:30], lty=2)
  segments(1:9+.4,post[3,22:30],2:10-.4,post[3,22:30], lty=2)
}


#figure 3 - saved manually
par(mfrow=c(2,1), oma=c(3,3,1,1), mar=c(0.5,1.5,0,0), font.main=1)
ageplot(yrrac2)
legend('topright',c('ICD9','ICD10'), lty=c(1,2), pch=c(10,16), bty='n', cex=0.8)
mtext("All-Cause (Log RRA)", side=2, outer=F, line=2, cex=0.8)
ageplot(yrrdc2)
mtext("Underlying Cause (Log RRD)", side=2, outer=F, line=2, cex=0.8)
axis(1,at=1:9+.5, labels=colnames(x2[,2:10]))

#PPD year plots for 45, 65, and 85

#construct 84% intervals
ppd = ycrc2$ppd
y = ycrc2$ydata
xdat = ycrc2$xdata
yr = unique(xdat[,"Years"])


plim = list()
pmean = list()
yy = list()

for(i in yr){
  #lim = (xdat[,'65'] == 1 | xdat[,'70'] == 1 | xdat[,'75'] == 1 | xdat[,'80']==1) & xdat[,'Years'] == i & xdat[,'Complex'] ==1
  lim = xdat[,'Years'] == i 
  yy[[as.character(i)]] = mean(y[lim])
  plim[[as.character(i)]] = quantile(ppd[lim,], prob=c(0.08,0.82,0.025,0.975))
  pmean[[as.character(i)]] = mean(ppd[lim,])
}
low=do.call(rbind,plim)[,1]
up=do.call(rbind,plim)[,2]
low2=do.call(rbind,plim)[,3]
up2=do.call(rbind,plim)[,4]

#vertlim = c(min(low2,yy)-.5,max(up2,yy)+.5)
vertlim=c(-3.2,0)

plot(unlist(pmean), ylim=vertlim, pch=15, xaxt="n", ylab="Log All-Cause Ratio", xlab="Year")
lines(unlist(yy), type="p",pch=0)
lines(low,type="l", lty=3)
lines(up,type="l",lty=3)
#lines(up2,lty=2)
#lines(low2,lty=2)

plim = list()
pmean = list()
yy = list()

for(i in yr){
  lim = xdat[,'65'] == 1 & xdat[,'Years'] == i & xdat[,'Complex'] == 0
  yy[[as.character(i)]] = mean(y[lim])
  plim[[as.character(i)]] = quantile(ppd[lim,], prob=c(0.08,0.82,0.025,0.975))
  pmean[[as.character(i)]] = mean(ppd[lim,])
}
low=do.call(rbind,plim)[,1]
up=do.call(rbind,plim)[,2]
low2=do.call(rbind,plim)[,3]
up2=do.call(rbind,plim)[,4]

lines(unlist(pmean), type="p", pch=16)
lines(unlist(yy), type="p",pch=1)
lines(low,type="l", lty=2)
lines(up,type="l",lty=2)
#lines(up2,lty=2)
#lines(low2,lty=2)

abline(v=which(yr==0)-0.5, lty=3)
axis(1,yr+1999,at=1:length(yr))
legend('topleft',c('Complex','Not Complex'), pch=c(15,16), bty='n')
#abline(h=0)

#@@@@
#simple plots
#@@@@

yrrac2$ppd

#@@@@@@@
#PPD ANOVA checks
#@@@@@@@

#are ICD10 ppd's more negative than ICD9?
plot(density(yrrac2$ppd[yrrac2$xdata[,'ICD10']==0,]))
lines(density(yrrac2$ppd[yrrac2$xdata[,'ICD10']==1,]), lty=2)

#omnibus
t.test(yrrac2$ppd[yrrac2$xdata[,'ICD10']==0],yrrac2$ppd[yrrac2$xdata[,'ICD10']==1,])
t.test(yrrdc2$ppd[yrrdc2$xdata[,'ICD10']==0],yrrdc2$ppd[yrrdc2$xdata[,'ICD10']==1,])
t.test(ycrc2$ppd[ycrc2$xdata[,'ICD10']==0],ycrc2$ppd[ycrc2$xdata[,'ICD10']==1,])

age='75'

t.test(yrrac2$ppd[yrrac2$xdata[,'ICD10']==0 & yrrac2$xdata[,age]==1]
       ,yrrac2$ppd[yrrac2$xdata[,'ICD10']==1 & yrrac2$xdata[,age]==1,])

t.test(yrrdc2$ppd[yrrdc2$xdata[,'ICD10']==0 & yrrdc2$xdata[,age]==1],
       yrrdc2$ppd[yrrdc2$xdata[,'ICD10']==1 & yrrdc2$xdata[,age]==1,])

t.test(ycrc2$ppd[ycrc2$xdata[,'ICD10']==0 & ycrc2$xdata[,age]==1],
       ycrc2$ppd[ycrc2$xdata[,'ICD10']==1 & ycrc2$xdata[,age]==1,])




