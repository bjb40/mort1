#dev R 3.2.3; Rstan 2.8.2; RStudio 0.99.0467; Stan 2.8.0
#time series analysis using modified changepoint hierarchical model
#1-2016

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
#a = accidents (external)
#n = noncommunicable (chronic)
#c = communicaable (acute)
#u = unassigned/residual
raw = read.csv(paste(rawdir,"morttab.csv",sep=''))
#remove 2004-2010 for balanced timeframe before and after transition
raw = raw[raw[,'year']<2004,]

#proportions
library(dplyr)
library(reshape2)

##generate descriptive statistics for 40+
raw$agegrp=cut(raw$ager,breaks=c(seq(40,85,by=5),150),right=FALSE)


#include total N, and N "c" and "u"
#uc contains all deaths because each death is assigned one mutually
#exclusive category
describe = raw %>% 
  group_by(ICD10,race,pd,complex,female)  %>%
  mutate(tot = uc_c_Sum + uc_n_Sum + uc_a_Sum + uc_u_Sum,
         uc.acute = uc_c_Sum, uc.chronic = uc_n_Sum,
         uc.external = uc_a_Sum, uc.residual = uc_u_Sum,
         ac.acute = any_c_Sum, ac.chronic = any_n_Sum,
         ac.external = any_a_Sum, ac.residual = any_u_Sum) %>% 
  select(-ager,-year,-matches('Sum'),-agegrp) %>% 
  summarize_each(funs(sum)) %>%
  ungroup %>%
  group_by(ICD10) %>% mutate(all.icddeaths=sum(tot)) %>%
  ungroup

##build descriptives table -- just put in by hand...

props=describe %>% group_by(ICD10) %>% 
  summarize_each(funs(sum),tot,matches('[auc]+\\.')) %>% 
  mutate_each(funs(rnd(./tot)),-ICD10,-tot) %>%
  ungroup

print(t(props))

describe$race = factor(describe$race,labels=c('white','black','other'))
describe$pd= factor(describe$pd,labels=c('Hospital','Long-term Care','Home','Other'))

#helper function building proprortions from the frequencies
describe.props = function(...){
  props = describe %>% group_by_('ICD10',...) %>%
     summarize_each(funs(sum),tot) %>% group_by(ICD10) %>% mutate(tot=rnd(tot/sum(tot)))
  return(props%>%ungroup)
}

for(i in c('race','female','pd','complex')){
  print(describe.props(i))
}

race.describe = describe %>% 
  filter(race %in% c('black','white')) %>%
  group_by(female,race,complex,pd) %>% 
  summarize_each(funs(sum),tot) %>%
  group_by(female,race) %>% mutate(tot=tot/sum(tot))

#confirm proportion sums...
race.describe %>% 
  group_by(female,race) %>% 
  summarize(sum(tot))

race.describe$race = factor(race.describe$race,labels=c('White','Black'))
race.describe$female = factor(race.describe$female,labels=c('Female','Male'))
race.describe$complex = factor(race.describe$complex,labels=c('1 or 2 causes','More than 2'))

library(gridExtra)

p=ggplot(race.describe, aes(y=tot,x=complex,fill=pd)) +
  geom_bar(stat='identity',position='stack') +
  facet_grid(.~race+female) + 
  theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5)) +
  ylab('Proportion') +
  scale_fill_discrete(guide=guide_legend(title='Place of Death'),l=50) +
  ggtitle('Place of Death:\nRace, Gender, and Comorbidities') +
  theme(axis.title.x=element_blank())

jpeg(paste0(imdir,'place_of_death.jpg'),
     height=1554,width=1488,units='px',res=170)
grid.arrange(p,bottom=textGrob('Source: National Center for Health Statistics “Public Use Multiple Cause of Death Mortality File" (1994-2003).
Calculated: Bryce Bartlett (http://www.brycebartlett.com).',
just='left',x=0,y=0.5,gp=gpar(fontsize=8)))
  
dev.off()


#average age:
ages = raw %>% group_by(ager,ICD10) %>%
  mutate(tot = uc_c_Sum + uc_n_Sum + uc_a_Sum + uc_u_Sum) %>% 
  select(tot,ager,ICD10) 

ages %>% group_by(ICD10) %>% 
  summarize(wtmn=funs(weighted.mean(.,tot)),ager)

print(
ages %>% group_by(ICD10) %>% mutate(wt = tot/sum(tot)) %>%
  summarize(wt.mn=weighted.mean(ager,wt),
            wt.sd=sqrt(sum(wt*(ager-wt.mn)^2)),
            mn=mean(ager), sd=sd(ager),
            wts = sum(wt)) #wt should be 1; just double-checks weight
  
)

#proportions by race and gender
cmplx=describe %>% filter(race=='black' | race=='white') %>%
  group_by(female,race,ICD10,pd) %>% 
  summarize_each(funs(sum),tot) %>% 
  group_by(female,race,ICD10) %>% mutate(tot=tot/sum(tot))

#View(cmplx)
cmplx$ICD10=factor(cmplx$ICD10,labels=c('Icd-9','Icd-10')) #RIGHT MAPPING?
cmplx$female=factor(cmplx$female,labels=c('Male','Female')) #RIGHT MAPPING?
#cmplx$complex=factor(cmplx$complex,labels=c('2 or fewer conditions','More than 2')) #RIGHT MAPPING?
cmplx$pd=factor(cmplx$pd,labels=c('Hospital','Home','LT Care','Other'))

ggplot(cmplx, aes(x=pd,y=tot)) + 
  geom_bar(stat='identity') + 
  facet_grid(race~female+ICD10,as.table=TRUE) 


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

#create weight matrix -- WRONG WEIGHTS -- WOULD NEED POPULATION WEIGHTS, NOT DEATH WEIGHTS
ndeaths=cbind(raw[,ids],apply(raw[,c('uc_c_Sum','uc_n_Sum','uc_a_Sum','uc_u_Sum')],1,sum))
ndeaths=aggregate(ndeaths[,length(ndeaths)],by=ndeaths[,ids],sum)
colnames(ndeaths)[length(ndeaths)] = 'tdeaths'
pool$tdeaths = ndeaths$tdeaths

#calculate DV concepts and expand race/placdth variables
#can log an calculate RRAc without much trouble; the others, not so much
RRAc = pool[,'any_c_Sum']/pool[,'any_n_Sum']
RRDc = pool[,'uc_c_Sum']/pool[,'uc_n_Sum']
CRc = RRDc/RRAc
tdeaths=pool$tdeaths


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
tdeaths = pool$tdeaths[lim]

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

colnames(x1) = bnames1
x1=data.frame(x1)

#clean up
rm(ndeaths,key,agedum,pool,raw,black,ltcare,lim,oplace,other,y,yrs,yrsxicd,CRc,RRAc,RRDc,as,home,icd10)

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

statadat = cbind(tdeaths,yrrac,yrrdc,x1)

write.csv(statadat,paste0(outdir,'stata-series.csv'))

rm(statadat,cells,cellsub,equals,i)

#@@@@@@@@@@@@@@@@@@@
#Run and Collect Models
#@@@@@@@@@@@@@@@@@@@
chains=4

#iters = 1500
iters = 2500
#iters = 5000
#iters = 7500
#iters=10000

burn=iters/2
#burn=2500

thin=1
#thin=5


#load rstan
library('rstan')

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(mc.cores = chains) #leave one core free for work

#@@@@@@
#Model 1 yrrac
#@@@@@@

y=yrrac
id=x1[,'cell']
t=x1[,'Years']
z=as.matrix(x1[,!colnames(x1) %in% c('yrrac','yrrdc','Years','cell','Intercept','x','X','tdeaths')])
#z=z[,1:9]
N=length(y)
IDS=length(unique(id))
P = ncol(z)
td = t #initialize
  td[t<0] = 1 #icd9
  td[t>=0] = 2 #icd10
TDS=length(unique(td))
YRS=length(unique(t))
yrctr = 6 #number to recenter to start at value of 1 for indexing

dnames = c('y','id','t','z','N','IDS','P','td','TDS','YRS','yrctr')

#save data used for yrrdc models
yrracdat = lapply(dnames,get,envir=.GlobalEnv) #otherwise defaults to base
names(yrracdat) = dnames
save(yrracdat,file=paste0(outdir,'yrracdat.RData'))
#rm(yrracdat)


yrrac1 = stan("bhm.stan",
              data=yrracdat,
              seed=1404399575,
              warmup=burn,
              chains=chains,
              iter=iters,
              thin=thin,
              verbose=T);

samp = extract(yrrac1,pars=c('beta','gamma','zi','delta','sig','loglik','dev','ppd'))
save(samp,file=paste0(outdir,'m1samp.gz'),compress=T)

sink(paste0(outdir,'stan-m1.txt'))
elapsed = get_elapsed_time(yrrac1)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrac1,pars=c('beta','gamma','zi','delta','sig'))
cat('Seed:      \t\t\t',get_seed(yrrac1))
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

rm(yrrac1,samp)

yrrac2 = stan("bhm-changepoint.stan",
              data=yrracdat,
              seed=1404399575,
              warmup=burn,
              thin=thin,
              chains=chains,
              iter=iters,
              verbose=F);

samp = extract(yrrac2,pars=c('beta','gamma','zi','delta','sig','loglik','dev','ppd','L_Omega','mu_i','mu_t'))
save(samp,file=paste0(outdir,'m2samp.gz'),compress=T)

sink(paste0(outdir,'stan-output-m2.txt'))

elapsed = get_elapsed_time(yrrac2)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrac2,pars=c('beta','gamma','sig','delta','zi','L_Omega'))
cat('Rhat range:\t\t\t',round(range(summary(yrrac2)$summary[,'Rhat']),3))

cat('\nwWarmup:\t\t\t',burn)
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

rm(yrrac2,samp)

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

#save data used for yrrdc models
yrrdcdat = lapply(dnames,get,envir=.GlobalEnv) #otherwise defaults to base
names(yrrdcdat) = dnames
save(yrrdcdat,file=paste0(outdir,'yrrdcdat.RData'))
rm(yrrdcdat)

yrrdc1 = stan("bhm.stan",
              data=yrrdcdat,
              warmup=burn,
              seed=1404399575,
              chains=chains,
              iter=iters,
              verbose=F);

samp = extract(yrrdc1,pars=c('beta','gamma','zi','delta','sig','loglik','dev','ppd'))
save(samp,file=paste0(outdir,'m3samp.gz'),compress=T)

sink(paste0(outdir,'stan-output-m3.txt'))

elapsed = get_elapsed_time(yrrdc1)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrdc1,pars=c('beta','gamma','zi','delta','sig'))
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

rm(yrrdc1,samp)

yrrdc2 = stan("bhm-changepoint.stan",
              data=yrrdcdat,
              seed=1404399575,
              warmup=burn,
              chains=chains,
              iter=iters,
              verbose=F);

samp = extract(yrrdc2,pars=c('beta','gamma','zi','delta','sig','loglik','dev','ppd','L_Omega','mu_i','mu_t'))
save(samp,file=paste0(outdir,'m4samp.gz'),compress=T)

sink(paste0(outdir,'stan-output-m4.txt'))

elapsed = get_elapsed_time(yrrdc2)
elapsed = max(rowSums(elapsed))/60 #minutes elapsed

sum=summary(yrrdc2,pars=c('beta','gamma','zi','delta','sig','L_Omega'))
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
