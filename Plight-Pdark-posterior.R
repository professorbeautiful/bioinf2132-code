logit = function(P)  log(P/(1-P))
antilogit = function(Z) 1 - 1/(1+exp(Z))
plotPlightPdarkPosterior = function() {
  sig11 = sig12 = sig21 = matrix(c(tau+phi,tau,tau,tau+phi),nrow=2)  #### <-- error was here, off-diagonals.
  varhat = apply(DLdata, 1, function(r)sum(1/(r+1/2)))
  sig22 = diag(phi+varhat) + tau   ## marginal variance of the data?
  logit.prior.mean = c(mu0,mu0)
  postmean.logit = logit.prior.mean + sig12%*%solve(sig22) %*% (logit.hat-logit.prior.mean) 
  postmean.p = antilogit(postmean.logit)
  postvar.logit = sig11 - sig12%*%solve(sig22)%*%sig21
  
  ####################
  
  svdSig11 = svd(sig11)
  sqrtSig11 = svdSig11$v %*% diag(sqrt(svdSig11$d)) %*% svdSig11$u
  svdPostvar.logit = svd(postvar.logit)
  sqrtPostvar.logit = svdPostvar.logit$v %*% diag(sqrt(svdPostvar.logit$d)) %*% svdPostvar.logit$u
  
  angles = seq(0,2*pi,length=100)
  qlevel = qchisq(p=0.99, df=2)
  circlepoints = cbind(cos(angles),sin(angles))
  plot(0:1, 0:1, xlab = "P dark", ylab = "P light", pch=" ")
  for (plevel in seq(.1,.9,.1)) {
  	qlevel = qchisq(p=plevel, df=2)
  	contourlevel1 = antilogit(logit.prior.mean 
  			+ (circlepoints*qlevel) %*% sqrtSig11)
  	lines(contourlevel1, col=ColorForPrior)
  	contourlevel2 = antilogit(matrix(c(postmean.logit),nr=100,nc=2, byrow=T)
  		+ (circlepoints * qlevel ) %*% sqrtPostvar.logit)
  	lines(contourlevel2, col=ColorForPosterior, lty=1, lwd=1)
  	if(plevel==.5){
      contourlevel1 = antilogit( logit.prior.mean 
      	+ (cbind(cos(angles),sin(angles)) * qlevel ) %*% sqrtSig11)
      xorder = order(contourlevel1[ ,1])
      lines(contourlevel1[xorder,1], contourlevel1[xorder,2],col = ColorForPrior, lty=2)
      contourlevel2 = antilogit(matrix(c(postmean.logit),nr=100, nc=2, byrow=T)
      	+ (cbind(cos(angles),sin(angles)) * qlevel ) %*% sqrtPostvar.logit)
      xorder = order(contourlevel2[ ,1])
      lines(contourlevel2[xorder,1], contourlevel2[xorder,2], col = ColorForPosterior,lty=2)
    }
  }
  points(antilogit(mu0), antilogit(mu0),pch = "M", cex=3, col=ColorForPrior)
  points(antilogit(postmean.logit[1,1]), antilogit(postmean.logit[2,1]), pch = "M", cex=3, col = ColorForPosterior)
  points(antilogit(logit.hat)[1], antilogit(logit.hat)[2],	pch = "X", cex=3)
}


DLdata = matrix(c(3,5,2,90),nrow=2)
dimnames(DLdata) = list(c("dark","light"),c("R","N"))
logit.hat = logit(apply(DLdata, 1, function(r)r[1]/sum(r)))
mu0 = logit(0.50)  ### prior mean.
#
#  Now, some spacing to make "Compile notebook.." come out nice.
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

ColorForPrior="green";     ColorForPosterior="red"

  ### tau = var of deviation shared by D and L.
  ### phi = var of individual deviations for D and L.




### Split:
tau = 0; phi = 1   ### Split:  D unconnected to L.
plotPlightPdarkPosterior()
title("Split: pL and pD are independent\n tau = 0; phi = 1")
#
#
#
#
#
#
### Lump:
tau = 1; phi = 0.002   ### Lump:  no individual variation:   D is same as L.
plotPlightPdarkPosterior()
title("Lump: pL = pD\n  tau = 1; phi = 0.002 ")
#
#
#
#
#
#
### Maybe:  a DIFFERENT compromise between Lump and Split
tau = 1/2;   phi = 1/2   ### D and L are somewhat connected.
plotPlightPdarkPosterior()
title("A compromise between Lump and Split:\n  tau = 1/2;   phi = 1/2")

