
#bryce bartlett
#converting correlation and covariances

sd.x = 7
sd.y = .239
sd.z = .5
#correlaton between -1 and 1
cor.xy =.7
cor.xz = -.3
cor.yz = 0

cor=matrix(c(1,cor.xy,cor.xz,
             cor.xy,1,cor.yz,
             cor.xz,cor.yz,1)
           ,3,3)

cov = matrix(0,3,3)
cov[1,1] = sd.x^2
cov[1,2] = cov[2,1]= cor.xy*(sd.x*sd.y)
cov[1,3] = cov[3,1] = cor.xz*(sd.x*sd.z)
cov[2,3] = cov[3,2] = cor.yz*(sd.y*sd.z)
cov[2,2] = sd.y^2
cov[3,3] = sd.z^2

#variance test
print('variance test')
var=sqrt(diag(cov))
print(var[1]==sd.x)
print(var[2]==sd.y)
print(var[3]==sd.z)

#corellation correlation with variance
print('making covariance with correlation and variance')
print(diag(var)%*%cor%*%diag(var) == cov)

#covariance to correlation with variance
print('making correlation out of covariance')
dc=solve(diag(var)) %*% cov %*% solve(diag(var))
print(dc==cor)

#working with the cholesky is computationally
#more efficient, because it is triangular -- same derivations
cov_R=chol(cov)
cov == t(cov_R) %*% cov_R
cor_R=chol(cor)
cor == t(cor_R) %*% cor_R


cov_R %*% solve(diag(var)) == cor_R; 

cor_R %*% diag(var) == cov_R