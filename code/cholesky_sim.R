#test cholesky

require(MASS)

#variance of x and y
sd.x = .230
sd.y = .239
#correlaton between -1 and 1
cor.xy =.5

sg = matrix(0,2,2)
sg[1,1] = sd.x^2
sg[1,2] = sg[2,1]= cor.xy*(sd.x*sd.y)
sg[2,2] = sd.y^2

#note that R does upper triangular decomp; whereas Stan does lower triangular
#lower triangular in R is simply the transpose of the upper triangular, i.e. t(R_chol)

R_chol = chol(sg); L_chol = t(R_chol)
mu=matrix(c(3.2,-.4),1,2)

mvdraw=mvrnorm(n=1000,mu=mu,Sigma=sg)

#Right cholesky
chdraw=matrix(rnorm(2000,mean=0,sd=1),1000,2)%*% R_chol + t(matrix(mu,2,1000))

#left cholesky

#correlation cholesky
cormat = matrix(c(1,cor.xy,cor.xy,1),2,2)
R_Omega=chol(cormat); L_Omega = t(R_Omega)
s = c(sd.x,sd.y)

eta=matrix(rnorm(2000,mean=0,sd=1),2,1000)
beta = t(matrix(mu,2,1000)) + t(diag(s) %*% L_Omega %*% eta)

#NOTE THE EQUALITY FOR WHY THE BETA WORKS L_Omega is chol of correlation; L_chol is chol of cov
print(diag(s) %*% L_Omega)
print(L_chol)

plot(mvdraw)
lines(chdraw,type='p',pch=3)
lines(beta,type='p',pch=4)
abline(h=mean(mvdraw[,2]),v=mean(mvdraw[,1]))
abline(h=mean(chdraw[,2]),v=mean(chdraw[,1]), col='gray')
abline(h=mean(beta[,2]),v=mean(beta[,1]), col='gray')

print(apply(mvdraw,2,FUN= function(x) c(mean=mean(x), sd=sd(x))))
print(apply(chdraw,2,FUN= function(x) c(mean=mean(x), sd=sd(x))))
print(apply(beta,2,FUN= function(x) c(mean=mean(x), sd=sd(x))))

#print("Chol' * Chol")
#print(t(cholsg) %*% cholsg)
#print("Sigma")
#print(sg)
#print("diag(s) * cholcor' * cholcor * diag(s)")
#diag(s) %*% t(cholcor) %*% cholcor %*% diag(s)