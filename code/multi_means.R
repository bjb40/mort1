
rm(list=ls())

source('funs.R')

#number of groups
j=5

#balanced design -- equal numbers in group
#should loosen this assumption
n=50

#draw means relative to group 1
grps=matrix(0,j,3); colnames(grps)=c('mean','sd','se')
grps[,1] = runif(j,-1,1)
grps[,2] = runif(j,.5,1.5)
grps[,3] = grps[,'sd']/sqrt(n)


#Draw samples of equal n from groups
#Place in dataframe 
grpsim=data.frame(x=numeric(),group_no=integer())

for(g in 1:j){
  grpsim=rbind(grpsim,
    data.frame(
      x=rnorm(n,mean=grps[g,'mean'],sd=grps[g,'sd']),
      group_no=rep(g,n))
  )
}


#TEST: for checking to make sure draw is okay
grpsim.sum = aggregate(grpsim$x,
                  by=list(grpsim$group_no),
                  function(g)
                    c(mean=mean(g),
                      sd=sd(g),
                      se=se(g))
                  )
"
print('Should be close to equal.')
print(grpsim.sum)
print(grps)
"

"
#TEST calculation of gamma
#real, i.e. gamma constructed from sample parameters
r.gamma=matrix(0,j,j)
#should set out to go by ratio of s.e. -- Figures 1 and 2 range from 0 to 100
#these are nonlinear but mirrored in scaled effect (Fig 2)
for(r in 1:j){r.gamma[r,]=grps[,3]/grps[r,3]}
r.gamma.bar = mean(c(r.gamma[lower.tri(r.gamma)],r.gamma[upper.tri(r.gamma)]))

print('real gamma')
print(r.gamma)
print(r.gamma.bar)
print('test function')
print(calc_gamma(grps[,3]))
#"


#revise thise lines into into adj_ci function


adj_ci = naive_ci=as.data.frame(aggregate(grpsim$x,
                   by=list(grpsim$group_no),
                   ci)[[2]])

gamma = calc_gamma(naive_ci$se)
print('Ratios of SE (row-wise comparisons.')
print(gamma)
adj_ci$gamma=calc_gamma(adj_ci$se)$ref

adj=function(crit,gamma){
  
  return(
    crit*((sqrt(1+gamma^2)-1)/gamma)
    )
  
}

adj_ci$adj_crit[adj_ci$gamma != 1] = 
  adj(crit=adj_ci$crit[adj_ci$gamma != 1]
      ,gamma=adj_ci$gamma[adj_ci$gamma != 1])

adj_ci$adj_crit[adj_ci$gamma==1] = adj_ci$crit[adj_ci$gamma==1]

#print(adj_ci)

#ref -- lynch & bartlett
ref_ci=
    lapply(unique(grpsim$group_no),
       FUN=function(g)
         ci(
           var=grpsim$x[grpsim$group_no==g],
           crit=adj_ci$adj_crit[g]
         )
         )#end lapply

#mean -- alt recommendation
adj_ci$adj_mean = adj(adj_ci$crit,gamma$mean)

mean_ci=
  lapply(unique(grpsim$group_no),
         FUN=function(g)
           ci(
             var=grpsim$x[grpsim$group_no==g],
             crit=adj_ci$adj_mean[g]
           )
  )#end lapply

#simplify cis
naive_ci=simplify2array(naive_ci)
ref_ci=t(simplify2array(ref_ci))
mean_ci=t(simplify2array(mean_ci))

#print(naive_ci)
#print(ref_ci)
#print(mean_ci)

##end lines to place in function--with ref return


###
#make some plots
###
cinms=c('lower','upper')
lims = range(c(
  naive_ci[,cinms],
  mean_ci[,cinms],
  ref_ci[,cinms]
))

alphas = sprintf('%.2f',round(ref_ci[,'alpha'],2))


plot(naive_ci[,'mean'],ylab='Mean',xlab='Group',
     ylim=lims,xlim=c(0.5,j+0.5),xaxt='n')

  #95% c.i.
arrows(1:j,naive_ci[,'lower'],
       1:j,naive_ci[,'upper'],
       angle=90,code=3,length=.05)

segments((1:j)+.0625,ref_ci[,'lower'],
         (1:j)+.0625,ref_ci[,'upper'],lty=2)

segments((1:j)-.0625,mean_ci[,'lower'],
         (1:j)-.0625,mean_ci[,'upper'],lty=3)

for(a in (1:j)){
text(a,
     naive_ci[a,'mean'],
     labels=bquote(paste(alpha,'=',.(alphas[a]) )),
     pos=4,
     cex=.6)
}

legend('topleft',
       legend=c(
         expression(paste('Naive (',alpha,'=0.95)')),
         expression(paste('Ref (','variable ',alpha,')')),
         bquote(paste('Mean (', alpha, '=', .(round(mean_ci[1,'alpha'],2)),')'))
         ),
       bty='n',
       lty=c(1,2,3),
       cex=.7
       )

axis(1,at=1:j)

#plot sim to fig 1
##another plot
ratios=c(seq(.01,100,by=.01))
#select 1 becuase then this erases the mult of crit val to recover error
plot(ratios,adj(1,ratios),type='l',xlab=expression(gamma),ylab='Z Adjustment Ratio')
