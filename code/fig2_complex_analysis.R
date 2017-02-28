####
#adding complex
###

ageplot = function(mod){
  #input approppiate model number; must be 2 or 4
  #pull effects and 84% intervals; include predicted effects for complex = 1
  complex = model[[mod]]$gamma[,,12]

  yrrac9 = apply(model[[mod]]$gamma[,1,1:10],2,FUN=function(x) eff(exp(x+complex[,1]), c=.84))
  yrrac10 = apply(model[[mod]]$gamma[,2,1:10],2,FUN=function(x) eff(exp(x+complex[,2]), c=.84))
  deltas=apply(model[[mod]]$gamma[,,1:10],3,FUN=function(x) 
                 delta(exp(x+complex),pval=TRUE))
  
  sigs=unlist(lapply(deltas['pval',],sig))
  sigs[sigs=='+  '] = '   ' #get rid of 0.1 level
  sigs=trimws(sigs) #trim whitespace for centering
  

  yl=c(range(c(yrrac9,yrrac10)))
  
  plot(1,type='n',ylim=yl,xlim=c(0.5,10.5),log='y',axes=FALSE)

  
  segments(0.5:9.5,yrrac9['mean',],1.5:10.5,yrrac9['mean',])
  text(1:10,
       (yrrac10['mean',])-(deltas['mean',]/2),
       labels=sigs,cex=.6)
  
 segments(0.5:9.5,yrrac10['mean',],1.5:10.5,yrrac10['mean',],lty=2)

}

print('Complex comparison 1999')
par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(2,2,2,0), font.main=1)
ageplot(2)
legend('topright',bty='n',lty=c(1,2),legend=c('ICD-9','ICD-10'))
title(main='All Cause',outer=FALSE)
axis(2)

ageplot(4)
title(main='Underlying Cause')
axis(2)
axis(1,labels=paste(seq(40,85,by=5),seq(44,89,by=5),sep='-'), 
     at=1:10)  

#Sys.sleep(10)
####
#comparing years of change (not complex)
###

ageplot = function(mod,icd9year,icd10year){
  #input approppiate model number; must be 2 or 4
  #pull effects and 84% intervals; include predicted effects for complex = 1
  complex = model[[mod]]$gamma[,,12]
  #because comparing 1998 yrrac (add year*-1) 
  #and 1999 yrac10 (year*0) -- only need to adjust year for icd9
  year = model[[mod]]$beta[,1]*icd9year
  year = cbind(year,model[[mod]]$beta[,2]*icd10year) + complex
  #complex = matrix(0,2,2)
  yrrac9 = apply(model[[mod]]$gamma[,1,1:10],2,FUN=function(x) eff(exp(x+year[,1]), c=.84))
  yrrac10 = apply(model[[mod]]$gamma[,2,1:10],2,FUN=function(x) eff(exp(x+year[,2]), c=.84))
  deltas=apply(model[[mod]]$gamma[,,1:10],3,FUN=function(x) 
    delta(exp(x+year),pval=TRUE))
  sigs=unlist(lapply(deltas['pval',],sig))
  sigs[sigs=='+  '] = '   ' #get rid of 0.1 level
  sigs=trimws(sigs) #trim whitespace for centering
 
  yl=c(range(c(yrrac9,yrrac10)))
  
  plot(1,type='n',ylim=yl,xlim=c(0.5,10.5),log='y',axes=FALSE)

  
  segments(0.5:9.5,yrrac9['mean',],1.5:10.5,yrrac9['mean',])
  text(1:10,
       (yrrac10['mean',])-(deltas['mean',]/2),
       labels=sigs,cex=.6)
  

  segments(0.5:9.5,yrrac10['mean',],1.5:10.5,yrrac10['mean',],lty=2)

  
}

#print('Complex comparison 1998 and 1999 (not complex)')

pdf(paste0(imdir,'full-ageplots-complex.pdf'))

par(mfrow=c(5,2), oma=c(0,0,0,2), mar=c(2,2,2,2), font.main=1)

for(y in -5:-1){
ageplot(2,y,0)
  y.text=1999+y
if(y==-5){  

title(main=paste0('All Cause'),outer=FALSE)
}
axis(2)
mtext(paste0(y.text,' (ICD-9)'),
      side=4,
      las=2,
      cex=.6,
      adj=-0.2)

if(y==-1){

  axis(1,labels=paste(seq(40,85,by=5),seq(44,89,by=5),sep='-'), 
       at=1:10)  
  
}
  
  
ageplot(4,y,0)

if(y==-5){
  title(main='Underlying Cause')
  legend('topright',bty='n',lty=c(1,2),
         legend=c('ICD-9','ICD-10'))
  }
axis(4)

if(y==-1){
axis(1,labels=paste(seq(40,85,by=5),seq(44,89,by=5),sep='-'), 
     at=1:10)  
}
}
dev.off()