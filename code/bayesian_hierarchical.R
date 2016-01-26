#exchangeable data across demographic profile and --- can add terms to make more similar similar


#rm(list=ls())

#load data

y=statadat[,'yrrac']
id=statadat[,'dem']
t=statadat[,'Years']
z=statadat[,!colnames(statadat) %in% c('yrrac','yrrdc','Years','dem','Intercept','x')]

N=length(y)
IDS=length(unique(id))
P = ncol(z)
td = t+6
TDS=length(unique(td))

#@@@@@@@@@@@@@@@@@
#call and send to stan 
#@@@@@@@@@@@@@@@@@

library('rstan')

#detect cores to activate parallel procesing
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
options(mc.cores = 3) #leave one core free for work

fit <- stan("bhm-cc.stan", data=c("y","id","td","t","z","P","N","IDS","TDS"),
            #algorithm='HMC',
            chains=3,iter=500,verbose=T);
            #sample_file = paste0(outdir,'diagnostic~/post-samp.txt'),
            #diagnostic_file = paste0(outdir,'diagnostic~/stan-diagnostic.txt'),
            

post_sum = summary(fit)[[1]] #independant of chain

traceplot(fit, pars=c('beta','sig','zi'))
traceplot(fit, pars=c('gamma'))


#sink(paste0(outdir,'bayes.txt'))
print(summary(fit,pars=c('beta','sig','zi','lp__'),digits=3))
print(summary(fit,pars=c('gamma'),digits=3))

#sink()
