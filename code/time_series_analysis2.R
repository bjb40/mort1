#dev R 3.1.3
#time series 

rm(list=ls())

#@@@@@@
#Universals
#@@@@@@

#functions
source('funs.R')

rawdir = "H:/projects/mort1/dat~/"
outdir = "H:/projects/mort1/output/"
imdir = "H:/projects/mort1/img~/"

#@@@@
#Load and clean data
#@@@@

#prepared in SAS per notes
raw = read.csv(paste(rawdir,"morttab.csv",sep=''))
#remove 2004-2010 for balanced timeframe before and after transition
raw = raw[raw[,'year']<2004,]

#recodes
a_c = raw$any_c_Sum
a_n = raw$any_n_Sum
d_c = raw$uc_c_Sum
d_n = raw$uc_n_Sum

#@@@
#exploring some dscriptives
#@@@

#ratios sliced by icd only
icd10 = raw$ICD10

#"o" for old = ICD9; "n" for new = ICD10
o = apply(cbind(a_c,a_n,d_c,d_n)[icd10==0,],2,sum)
n = apply(cbind(a_c,a_n,d_c,d_n)[icd10==1,],2,sum)

#desriptives
print(c(o['d_c']/o['d_n'],n['d_c']/n['d_n']))
print(c(o['a_c']/o['a_n'],n['a_c']/n['a_n']))


oCR = (o['d_c']/o['d_n'])/(o['a_c']/o['a_n'])
nCR = (n['d_c']/n['d_n'])/(n['a_c']/n['a_n'])


key = cbind(a_c,a_n,d_c,d_n)
ys = lapply(unique(raw$year), function(y) apply(key[raw$year == y,],2,sum))
CRy = lapply(1:length(ys), function(y) (ys[[y]]['d_c']/ys[[y]]['d_n'])/(ys[[y]]['a_c']/ys[[y]]['a_n']))
RRDy = lapply(1:length(ys), function(y) ys[[y]]['d_c']/ys[[y]]['d_n'])
RRAy = lapply(1:length(ys), function(y) ys[[y]]['a_c']/ys[[y]]['a_n'])

vals = c(unlist(RRDy),unlist(RRAy))
yrs = min(raw$year):max(raw$year)

png(paste(imdir,'rate_plot_new.png',sep=''))
plot(yrs,unlist(RRDy), ylim=c(min(vals),max(vals)+.075), xlab='', ylab='')
  lines(yrs,unlist(RRDy), lty=2)
  lines(yrs,unlist(RRAy), type='p')
  lines(yrs,unlist(RRAy), lty=3)
  abline(v=1998.5, lty=1)
  legend('topright',c('RRD','RRA'), lty=c(2,3), bty='n')
  text(1998.5,max(vals)+.075, labels="ICD 9", pos=2)
  text(1998.5,max(vals)+.075,labels="ICD 10", pos=4)
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@
#Additional data cleaning for model building
#@@@@@@@@@@@@@@@@@@@@@@@@@@


#pool by five year ages
as = seq(40,85,by=5)
raw$agegroup = as.numeric(NA)
for(j in 1:length(as)){
  raw$agegroup[raw$ager >= as[j]] = j 
}

table(raw$ager,raw$agegroup)

ids = c('ICD10','year','female','race','pd','complex','agegroup')
pool=aggregate(raw[,c('uc_c_Sum','uc_n_Sum','any_c_Sum','any_n_Sum')], by=raw[,ids],sum)

#calculate DV concepts and expand race/placdth variables
#can log an calculate RRAc without much trouble; the others, not so much
RRAc = pool[,'any_c_Sum']/pool[,'any_n_Sum']
RRDc = pool[,'uc_c_Sum']/pool[,'uc_n_Sum']
CRc = RRDc/RRAc

#proportion of structural zeros
print(sum(is.finite(CRc)==F)/length(CRc))

#prep analysis dataset

#recodes
oplace = home = ltcare = rep(0,nrow(pool))
ltcare[pool$pd == 1] = 1
home[pool$pd == 2] = 1
oplace[pool$pd == 3] = 1

black = other = rep(0,nrow(pool))
black[pool$race == 1] = 1
other[pool$race == 2] = 1
yrs = pool$year-1999
yrsxicd = yrs*pool$ICD10

agedum = matrix(0,length(pool$agegroup),length(unique(pool$agegroup)))
for(i in 1:length(unique(pool$agegroup))){
  agedum[pool$agegroup==i,i] = 1
}

#@@@@@@@@@@@@@@@@@@@@@
#Prepare Model 1
#@@@@@@@@@@@@@@@@@@@@@

#lots of structural zeros with "other" as the race, so "other race is eliminated"
lim = other != 1 #rep(T,length(CRc))#is.finite(CRc)

ycrc = log(CRc[lim]) 
yrrac = log(RRAc[lim])
yrrdc = log(RRDc[lim])

#sort of normal --- 3 structural zeros leading to infinite /set to maximum negative
for(y in list(yrrac,yrrdc,ycrc)){
  print(paste('Changing',sum(is.finite(y)==F),'structural zeros...'))
  y[y==-Inf] = min(y[y!=-Inf])-.5
}

x1 = cbind(
  #rep(1,length(y)),
  agedum[lim,1:ncol(agedum)],
  yrs[lim],
  pool$female[lim], 
  pool$complex[lim],
  home[lim],
  ltcare[lim],
  oplace[lim],
  black[lim]
)

bnames1=c(#'Intercept',
          #'Ages',
          paste0('a',(as[1:length(as)])),
          'Years',
          'Female',
          'Complex',
          'Home',
          'Ltcare',
          'Oplace',
          'Black'
)


x2 = cbind(
  rep(1,length(y)),
  #pool$agegroup[lim],
  agedum[lim,2:ncol(agedum)],
  pool$ICD10[lim],
  yrs[lim],
  yrsxicd[lim],
  pool$female[lim], 
  pool$complex[lim],
  home[lim],
  ltcare[lim],
  oplace[lim],
  black[lim],
  pool$ICD10[lim]*black[lim],
  pool$ICD10[lim]*pool$female[lim],
  pool$ICD10[lim]*agedum[lim,2:ncol(agedum)],
  pool$ICD10[lim]*pool$complex[lim]
)

colnames(x1) = bnames1
x1=data.frame(x1)


#clean up
rm(key,agedum,pool,raw,black,ltcare,lim,oplace,other,y,yrs,yrsxicd,CRc,RRAc,RRDc,as,home,icd10)


#export to csv for analysis in stata
#prepare 'panel' identifier based on observed characteristics
#dem = numeric(length=length(yrrac))


#output unique demographic slices
cells=unique(x1[,colnames(x1) %in% c(paste0('a',seq(45,85,by=5)),'Female','Black','Complex','Home','Ltcare','Oplace')])
cellsub = x1[,colnames(x1) %in% colnames(cells)]
x1$cell=as.numeric(NA)

for(i in 1:nrow(cells)){
  #test rowise equality
  if(i%%10 == 0){cat('Coding ',i,'of',nrow(cells),'\n')}
  equals = apply(cellsub,1,FUN = function(x) all(x == cells[i,]))
  #print(c(i,sum(equals))) #should be 10 equals each for each of the time periods
  x1[equals,'cell'] = i
}

statadat = cbind(yrrac,yrrdc,x1)

write.csv(statadat,paste0(outdir,'stata-series.csv'))

rm(statadat)

#@@@@@@@@@@@@@@@@@@@
#Run and Collect Models
#@@@@@@@@@@@@@@@@@@@

#load rstan
library('rstan')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(mc.cores = 3) #leave one core free for work

#@@@@@@
#Model 1 yrrac
#@@@@@@

y=yrrac
id=x1[,'cell']
t=x1[,'Years']
z=as.matrix(x1[,!colnames(x1) %in% c('yrrac','yrrdc','Years','cell','Intercept','x')])
#z=z[,1:9]
N=length(y)
IDS=length(unique(id))
P = ncol(z)
td = t #initialize
  td[t<0] = 1 #icd9
  td[t>=0] = 2 #icd10
TDS=length(unique(td))

#iters = 5000
iters = 1000

yrrac1 = stan("bhm.stan", data=c('y','id','t','z','N','IDS','P'),
               #algorithm='HMC',
               chains=3,iter=iters,verbose=T);

sink(paste0(outdir,'stan-m1.txt'))

elapsed = get_elapsed_time(yrrac1)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrac1,pars=c('beta','gamma','zi','sig'))
ieffects=summary(yrrac1,pars='mu_i')
cat('Rhat range:\t\t\t',round(range(summary(yrrac1)$summary[,'Rhat']),3))

cat('\nIterations:\t\t\t',iters)
cat('\nElapsed min:\t\t\t',round(elapsed,3))
cat('\nIters/minute:\t\t\t',round((iters/elapsed),3))
cat('\nn_eff (samp):\t\t\t',round(range(summary(yrrac1)$summary[,'n_eff']),3))
cat('\nn_eff/minutes:\t\t\t',round(range(summary(yrrac1)$summary[,'n_eff'])/elapsed,3))

cat('\n\nParms:\n')
print(round(sum$summary,3))

cat('\n\nFit:\n')
print(dic(yrrac1))
print(waic(yrrac1))


sink()

hist(ieffects$summary[,'mean'])

yrrac2 = stan("bhm-changepoint.stan", data=c('y','id','t','z','N','IDS','P','TDS','td'),
              #algorithm='HMC',
              chains=3,iter=iters,verbose=T);


sink(paste0(outdir,'stan-output-m2.txt'))

elapsed = get_elapsed_time(yrrac2)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrac2,pars=c('beta','gamma','zi1','zi2','sig'))
ieffects=summary(yrrac2,pars='mu_i')
cat('Rhat range:\t\t\t',round(range(summary(yrrac2)$summary[,'Rhat']),3))

cat('\nIterations:\t\t\t',iters)
cat('\nElapsed min:\t\t\t',round(elapsed,3))
cat('\nIters/minute:\t\t\t',round((iters/elapsed),3))
cat('\nn_eff (samp):\t\t\t',round(range(summary(yrrac2)$summary[,'n_eff']),3))
cat('\nn_eff/minutes:\t\t\t',round(range(summary(yrrac2)$summary[,'n_eff'])/elapsed,3))

cat('\n\nParms:\n')
print(round(sum$summary,3))

cat('\n\nFit:\n')
print(dic(yrrac2))
print(waic(yrrac2))

sink()

lim=is.finite(yrrdc)
y=yrrdc[lim]
id=id[lim]
t=t[lim]
z=z[lim,]
td=td[lim]

N=length(y)
P=ncol(z)
IDS=length(unique(id))
TDS=length(unique(td))


yrrdc1 = stan("bhm.stan", data=c('y','id','t','z','N','IDS','P'),
              #algorithm='HMC',
              chains=3,iter=iters,verbose=T);




sink(paste0(outdir,'stan-output-m3.txt'))

elapsed = get_elapsed_time(yrrdc1)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrdc1,pars=c('beta','gamma','zi','sig'))
ieffects=summary(yrrdc1,pars='mu_i')
cat('Rhat range:\t\t\t',round(range(summary(yrrdc1)$summary[,'Rhat']),3))

cat('\nIterations:\t\t\t',iters)
cat('\nElapsed min:\t\t\t',round(elapsed,3))
cat('\nIters/minute:\t\t\t',round((iters/elapsed),3))
cat('\nn_eff (samp):\t\t\t',round(range(summary(yrrdc1)$summary[,'n_eff']),3))
cat('\nn_eff/minutes:\t\t\t',round(range(summary(yrrdc1)$summary[,'n_eff'])/elapsed,3))

cat('\n\nParms:\n')
print(round(sum$summary,3))

cat('\n\nFit:\n')
print(dic(yrrdc1))
print(waic(yrrdc1))

sink()


yrrdc2 = stan("bhm-changepoint.stan", data=c('y','id','t','z','N','IDS','P','TDS','td'),
              #algorithm='HMC',
              chains=3,iter=iters,verbose=T);

sink(paste0(outdir,'stan-output-m4.txt'))

elapsed = get_elapsed_time(yrrdc2)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrdc2,pars=c('beta','gamma','zi1','zi2','sig'))
ieffects=summary(yrrdc2,pars='mu_i')
cat('Rhat range:\t\t\t',round(range(summary(yrrdc2)$summary[,'Rhat']),3))

cat('\nIterations:\t\t\t',iters)
cat('\nElapsed min:\t\t\t',round(elapsed,3))
cat('\nIters/minute:\t\t\t',round((iters/elapsed),3))
cat('\nn_eff (samp):\t\t\t',round(range(summary(yrrdc2)$summary[,'n_eff']),3))
cat('\nn_eff/minutes:\t\t\t',round(range(summary(yrrdc2)$summary[,'n_eff'])/elapsed,3))

cat('\n\nParms:\n')
print(round(sum$summary,3))

cat('\n\nFit:\n')
print(dic(yrrdc2))
print(waic(yrrdc2))

sink()




#delete three structural zeros
yrrdc1 = rungibbs(yrrdc[is.finite(yrrdc)],x1[is.finite(yrrdc),])
yrrdc2 = rungibbs(yrrdc[is.finite(yrrdc)],x2[is.finite(yrrdc),])

#delete three structural zeros
ycrc1 = rungibbs(ycrc[is.finite(ycrc)],x1[is.finite(ycrc),])
ycrc2 = rungibbs(ycrc[is.finite(ycrc)],x2[is.finite(ycrc),])


#saveRDS(yrrac2,file=paste(outdir,'yrrac2.zip',sep=''),compress=T)

#@@@@@@@@@@@@@@@@@@@
#Put together output
#@@@@@@@@@@@@@@@@@@@

summ = function(l,nms){
  colnames(l$betas) = nms
  l$est = apply(l$betas,2,mean)
  l$estsd = apply(l$betas,2,sd)
  l$pval = lapply((l$est/l$estsd),dnorm)
  l$sig = rep('',length(l$pval))
    names(l$sig) = names(l$sig)
    l$sig[l$pval<0.05] = '*'
    l$sig[l$pval<0.01] = '**'
    l$sig[l$pval<0.001] = '***'
 
  return(cbind(round(l$est, digits=3),paste('(',round(l$estsd,digits=3),')',sep=''),l$sig))
}

  m = list(yrrac1=summ(yrrac1,bnames1),
           yrrac2=summ(yrrac2,bnames2),
           yrrdc1=summ(yrrdc1,bnames1),
           yrrdc2=summ(yrrdc2,bnames2),
           ycrc1=summ(ycrc1,bnames1),
           ycrc2=summ(ycrc2,bnames2)
           )

#need to make a two-tail test and make sure  you get the right tails

sink(paste(outdir,'output-',Sys.Date(),'.txt',sep=''))
print(Sys.Date(),quote='F')
for(i in 1:length(m)){
  cat('@@@@@@@@@@@@@@@@@@@\n')
  print(names(m)[i], quote=F)
  cat('\n')
  cat('@@@@@@@@@@@@@@@@@@@\n\n')
  print(cbind(rownames(m[[i]]), apply(m[[i]],2,paste, '|')), quote=F)
  cat('\n\nR-Square\n')
  print(mean(get(names(m)[i])$rsq))
  print(quantile(get(names(m)[i])$rsq,prob=c(0.025,0.08,0.82,0.975)), quote=F)
  cat('\n\nobs:\n')
  print(length(get(names(m)[i])$y))
  cat('\n\n\n\n')
}

sink()

#bayesian pval for fit examples
sum(yrrac2$rsq<mean(yrrac1$rsq))/length(yrrac2$rsq)
sum(yrrdc2$rsq<mean(yrrdc1$rsq))/length(yrrdc2$rsq)
sum(ycrc2$rsq<mean(ycrc1$rsq))/length(ycrc2$rsq)

#@@@@@
#test of whether the combined years + yearsxicd10 = 0
#@@@@@

bpos = which(bnames2 %in% c('Years','YearsxICD10'));

#bayesian p-value
plot(yrrac2$betas[,bpos])

sl10 = apply(yrrac2$betas[,bpos],1,sum)
mean(sl10)
#bayesian pval
sum(sl10>0)/length(sl10)
#interval
quantile(sl10,prob=c(0.025,0.975))

sl10 = apply(yrrdc2$betas[,bpos],1,sum)
mean(sl10)
#bayesian pval
sum(sl10>0)/length(sl10)
#interval
quantile(sl10,prob=c(0.025,0.975))


#plot (SAVED MANUALLY)
y9 = y10 = matrix(NA,nrow(yrrac2$betas),5)

colnames(y9) = -5:-1
colnames(y10) = 0:4

yb = which(bnames2 == 'Years')
ib = which(bnames2 == 'Intercept')
yadd = which(bnames2 == 'YearsxICD10')
iadd = which(bnames2 == 'ICD10')

for (i in 1:5){  
  y9[,as.character(-1*i)] = yrrac2$betas[,yb]*(-1*i) + yrrac2$betas[,ib]
  y10[,as.character(i-1)] = yrrac2$betas[,yb]*(i-1) + yrrac2$betas[,yadd]*(i-1) + yrrac2$betas[,ib] + yrrac2$betas[,iadd]
}

par(mfrow=c(2,1), oma=c(3,3,1,1), mar=c(0.5,1.5,0,0), font.main=1)


plot(-5:-1,apply(y9,2,mean),xlim=c(-5,4), ylim=c(-1,-.4), type="l", lty=1, xaxt='n',ylab='All-Cause (Log RRA)',xlab='' )
    lines(-5:-1,apply(y9,2,quantile, prob=0.08), lty=2)
    lines(-5:-1,apply(y9,2,quantile, prob=0.92), lty=2)
    
lines(0:4,apply(y10,2,mean), lty=1)
    lines(0:4,apply(y10,2,quantile, prob=0.08), lty=2)
    lines(0:4,apply(y10,2,quantile, prob=0.92), lty=2)
    
abline(v=-0.5, lty=3)
mtext("All-Cause (Log RRA)", side=2, outer=F, line=2, cex=0.8)
text(-4.75,-0.45,"ICD9")    
text(0.25,-0.45, 'ICD10')    

    
#underlying cause
    
for (i in 1:5){  
      y9[,as.character(-1*i)] = yrrdc2$betas[,yb]*(-1*i) + yrrdc2$betas[,ib]
      y10[,as.character(i-1)] = yrrdc2$betas[,yb]*(i-1) + yrrdc2$betas[,yadd]*(i-1) + yrrdc2$betas[,ib] + yrrdc2$betas[,iadd]
    }
    
plot(-5:-1,apply(y9,2,mean),xlim=c(-5,4), ylim=c(-.8,.05), type="l", lty=1, xaxt='n',ylab='Underlying-Cause (Log RRD)',xlab='' )
    lines(-5:-1,apply(y9,2,quantile, prob=0.08), lty=2)
    lines(-5:-1,apply(y9,2,quantile, prob=0.92), lty=2)
    
lines(0:4,apply(y10,2,mean), lty=1)
    lines(0:4,apply(y10,2,quantile, prob=0.08), lty=2)
    lines(0:4,apply(y10,2,quantile, prob=0.92), lty=2)
    
abline(v=-0.5, lty=3)
axis(1,-5:4+1999,at=-5:4)
mtext("Underlying Cause (Log RRD)", side=2, outer=F, line=2, cex=0.8)
    
    
#@@@@@@
#(especially for complexity as an explanation)

#cases that are pretty bad are female and black
medest = rowSums(ycrc1$ppd-ycrc1$ydata>0)/ncol(ycrc1$ppd)
hist(medest)
#colnames(ycrc2$xdata) = bnames1
colMeans(x2[medest>.0275,])


#spot checking means looks okay
mean(ycrc1$ydata[ycrc1$xdata[,'Home'] == 1])
quantile(
  colMeans(ycrc1$ppd[ycrc1$xdata[,'Home']==1,]),
  prob=c(.025,.5,.975)
)


hist(ycrc2$ydata[ycrc2$xdata[,'ICD10'] == 1], prob=T, ylim=c(0,.75))
lines(density(ycrc2$ppd[ycrc2$xdata[,'ICD10'] == 1]), type="l", lty=2)

plot(density(yrrdc2$ydata[ycrc2$xdata[,'ICD10'] == 1]), prob=T, ylim=c(0,.75))
lines(density(yrrdc2$ydata[ycrc2$xdata[,'ICD10'] == 0]), prob=T, ylim=c(0,.75), lty=2)

hist(ycrc2$ydata[ycrc2$xdata[,'ICD10'] == 0])

#par(mfrow=1,2)
lim=ycrc2$xdata[,'Years'] < 0 & ycrc2$xdata[,'Black'] == 1 & ycrc2$xdata[,'65'] == 1
hist(ycrc1$ydata[lim],prob=T, ylim=c(0,.75))
lines(density(ycrc1$ppd[lim]), type="l", lty=2)
lines(density(ycrc2$ppd[lim]), type="l", lty=3)

#hist(yrrdc1$ydata,prob=T, ylim=c(0,.75))
#lines(density(yrrdc1$ppd), type="l", lty=2)



#print output---no need for burn-in because it's a gibbs sampler

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




