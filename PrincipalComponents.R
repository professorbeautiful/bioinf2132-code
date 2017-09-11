v1 = 1
v2 = 2
ro = .2
v12 <- ro*sqrt(v1*v2)
v = matrix(c(v1, v12, v12, v2), nrow=2)

pc = princomp(covmat=v)
pc[1:2]


svd.v = svd(v)
sqrt.v = svd.v$u %*% diag(sqrt(svd.v$d)) %*% svd.v$v
theta = seq(0, 2*pi, length=10000)
x = cbind(sin(theta), cos(theta)) %*% sqrt.v
plot(x[,1], x[,2], type="l")
require("mvtnorm")
summary(dmvnorm(x, c(0,0), v))

#abline(a=0, b=pc[[1]][1]/pc[[1]][2])
#plot(x[,1], x[,2], type="l", xlim=c(.8,1.05), 
#	ylim=c(1,1.4))
sqrtinv.v = svd.v$u %*% diag(1/sqrt(svd.v$d)) %*% svd.v$v
abline(a=0, b =  pc[[1]] %*% solve(sqrtinv.v), col="blue")
all.equal(sqrt(svd.v$d), pc[[1]])

##### pc$loadings de-correlates the data
 plot(x%*%pc$loadings)
 abline(v=pc$sdev[1]*c(-1,1))
 abline(h=pc$sdev[2]*c(-1,1))
##### pc$loadings are related to svd matrices
 print(svd.v$u)
 print(pc$loadings)

##### plot an equi-scalar plot!!
 myrange = c(-1,1)*max(x)
 plot(x, xlim=myrange, ylim=myrange, type="l", col="red")
 l2norm = function(x, dim=1) apply(x^2, dim, sum)
 thetamax = theta[l2norm(x)==max(l2norm(x))]
 thetamin = theta[l2norm(x)==min(l2norm(x))]
# plot(theta, l2norm(x))
# abline(v=thetamax)
# abline(v=thetamin)
 maxpoint = x[l2norm(x)==max(l2norm(x)), ]
 abline(a=0, b=maxpoint[2]/maxpoint[1], col="blue", lwd=4)					### major axis
abline(a=0, b=pc$loadings[2,1]/pc$loadings[1,1],col="green")  ### major axis
# lines(rbind(c(0,0), pc$sdev[2:1]), lwd=6)
 		### major axis--- NOT!
 minpoint = x[l2norm(x)==min(l2norm(x)), ]
 lines(rbind(c(0,0), minpoint))
 lines(rbind(maxpoint,minpoint))
 # abline(a=0, b=sqrt(v2/v1), col="green")  ## NOT the major axis
 
 ### draw a circle of radius max(l2norm)
 #symbols(0,0,circles=max(l2norm(x)), add=T)
 lines(sqrt(max(l2norm(x))) * cbind(cos(theta), sin(theta)))
 abline(h=0,v=0,lty=2)
 