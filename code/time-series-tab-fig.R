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
    #complex = model[[mod]]$gamma[,,12]
    complex = matrix(0,2,2)
    yrrac9 = apply(model[[mod]]$gamma[,1,1:10],2,FUN=function(x) eff(x+complex[,1], c=.84))
    yrrac10 = apply(model[[mod]]$gamma[,2,1:10],2,FUN=function(x) eff(x+complex[,2], c=.84))
    
    yl=range(c(yrrac9,yrrac10))
    
    plot(1,type='n',ylim=yl,xlim=c(1,10))
    o = .1 #offset for visibility
    
    lines((1:10)-o,yrrac9['mean',], type='p',pch=10)
    lines(yrrac9[1,])
    arrows((1:10)-o,yrrac9[2,],(1:10)-o,yrrac9[3,],angle=90,code=3,length=.1)
    lines((1:10)+o,yrrac10['mean',], type='p',pch=16)
    lines(yrrac10[1,],lty=2)
    arrows((1:10)+o,yrrac10[2,],(1:10)+o,yrrac10[3,],angle=90,code=3,length=.1)

}

par(mfrow=c(2,1), oma=c(3,3,1,1), mar=c(0.5,1.5,0,0), font.main=1)

ageplot(2)
ageplot(4)


d.all=apply(model[[2]]$gamma[,,1:10],3,delta)
d.under = apply(model[[4]]$gamma[,,1:10],3,delta)

yl = range(c(d.all,d.under))

par(mfrow=c(1,1))
plot(1,type='n',ylim=yl,xlim=c(1,10))

polygon(c(1:10, rev(1:10)), c(d.all[2,],rev(d.all[3,])), 
        col="gray90", border=NA)
polygon(c(1:10, rev(1:10)), c(d.under[2,],rev(d.under[3,])), 
        col="gray90", border=NA)

lines(1:10,d.all[1,], type="l",ylim=yl)
lines(1:10,d.under[1,], lty=2)

abline(h=0)



#@@@@@@@@@@@
#time series plot with ppd
#@@@@@@@@@@@

#summarize posterior draws
ppd.yrrac = model[[2]]$ppd






