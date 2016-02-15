

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

#snd = density(rnorm(1000,mean=0,sd=1))
#std = density(rt(1000,df=30))  
#plot(snd); lines(std,lty=2)

#snd = density(rnorm(1000,mean=5,sd=1))
#std = density(rt(1000,df=30,ncp=5))  
#plot(snd); lines(std,lty=2)

#snd = density(rnorm(1000,mean=0,sd=.5))

#scaled student_t -- see bda3, pp. 294,576
#draw as a mixture of normals
#alpha=30/2
#beta=alpha*(0.5)
#v_i = 1/rgamma(1000,shape=alpha,rate=beta)
#std = density(rnorm(1000,mean=0,sd=v_i))
#plot(std,lty=2);lines(snd) 


nd = density(rnorm(1000,mean=0,sd=delta))
alpha=9/2
beta=alpha*delta
v_i = 1/rgamma(1000,shape=alpha,rate=beta)
td = density(rnorm(1000,mean=0,sd=v_i))  
funt=density(rscaled_t(n=1000,df=9,s=delta))

plot(nd); lines(td,lty=2); lines(funt,lty=3)