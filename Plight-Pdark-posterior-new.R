plotPlightPdarkPosterior = function(
  DLdata = matrix(c(3,5,2,90),nrow=2),
  tau, phi, mu0,
  showPrior = TRUE, showPosterior=TRUE, showLikelihood=TRUE,
  fudgeFactor = 0.001,
  addFudge = TRUE
) {

  logit.hat = logit(apply(DLdata, 1, function(r)r[1]/sum(r)))
  varhat = apply(DLdata, 1, function(r)sum(1/r))
  if(any(abs(logit.hat) == Inf) | addFudge){
    DLdataFudged = DLdata + fudgeFactor  # if needed, if any zero cells.
    logit.hat.fudged = logit(apply(DLdataFudged, 1, function(r)r[1]/sum(r)))
    varhat.fudged = apply(DLdataFudged, 1, function(r)sum(1/r))
    logit.hat = logit.hat.fudged
    varhat = varhat.fudged
  }
  ### deltat method with protection from zero's.
  sig11 = sig12 = sig21 = matrix(c(tau+phi,tau,tau,tau+phi),nrow=2)  
  sig22 = sig11 + diag(varhat)   ## marginal variance of the data?
  logit.prior.mean = c(mu0, mu0)
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
  plot(0:1, 0:1, xlab = "Pr(R | D)", ylab = "Pr(R | L)", pch=" ",
       cex=2)
  if(showPrior) {
    for (plevel in seq(.1,.9,.1)) {
      qlevel = qchisq(p=plevel, df=2)
      contourlevel1 = antilogit(logit.prior.mean 
                                + (circlepoints*qlevel) %*% sqrtSig11)
      lines(contourlevel1, col=ColorForPrior)
      if(plevel==.5){
        contourlevel1 = antilogit( logit.prior.mean 
                                   + (cbind(cos(angles),sin(angles)) * qlevel ) %*% sqrtSig11)
        xorder = order(contourlevel1[ ,1])
        lines(contourlevel1[xorder,1], contourlevel1[xorder,2],col = ColorForPrior, lty=2)
      }
    }
    points(antilogit(mu0), antilogit(mu0),pch = "M", cex=3, col=ColorForPrior)
  }
  if(showPosterior) {
    for (plevel in seq(.1,.9,.1)) {
      qlevel = qchisq(p=plevel, df=2)
      contourlevel2 = antilogit(matrix(c(postmean.logit),nrow=100,ncol=2, byrow=T)
                                + (circlepoints * qlevel ) %*% sqrtPostvar.logit)
      lines(contourlevel2, col=ColorForPosterior, lty=1, lwd=1)
      if(plevel==.5){
        contourlevel2 = antilogit(matrix(c(postmean.logit),nrow=100, ncol=2, byrow=T)
                                  + (cbind(cos(angles),sin(angles)) * qlevel ) %*% sqrtPostvar.logit)
        xorder = order(contourlevel2[ ,1])
        lines(contourlevel2[xorder,1], contourlevel2[xorder,2], col = ColorForPosterior,lty=2)
      }
    }
    points(antilogit(postmean.logit[1,1]), antilogit(postmean.logit[2,1]), pch = "M", cex=3, col = ColorForPosterior)
  }
#  if(showLikelihood) {
#    argmin = function(v, target=0) which(abs(v-target) == min(abs(v-target))[1])
#     probs = seq(0.01, 0.99, by = 0.01)
#     normalizedLikelihoodDark = dbinom(x = DLdata["dark", "R"], size=sum(DLdata["dark", ]), 
#                                       prob=probs)
#     normalizedLikelihoodDark = normalizedLikelihoodDark/sum(normalizedLikelihoodDark)/0.01
#     cumLikelihoodDark = cumsum(normalizedLikelihoodDark)
#     normalizedLikelihoodLight = dbinom(x = DLdata["light", "R"], size=sum(DLdata["light", ]), 
#                                       prob=probs)
#     normalizedLikelihoodLight = normalizedLikelihoodLight/sum(normalizedLikelihoodLight)/0.01
#     cumLikelihoodLight = cumsum(normalizedLikelihoodLight)
#     for (plevel in seq(.1,.9,.1)) {
#       D1 = cumLikelihoodDark[argmin(cumLikelihoodDark, plevel)]
#       D2 = cumLikelihoodDark[argmin(cumLikelihoodDark, 1 - plevel)]
#       L1 = cumLikelihoodLight[argmin(cumLikelihoodLight, plevel)]
#       L2 = cumLikelihoodLight[argmin(cumLikelihoodLight, 1 - plevel)]
#       cat(paste("Likelihood ", plevel, D1, D2, L1, L2, "\n"))
#       lines(c(D1, D1, D2, D2), c(L1, L2, L2, L1), col= ColorForLikelihood)
#      }
#  }
  points(antilogit(logit.hat)[1], antilogit(logit.hat)[2],	pch = "X", cex=3)
  mtext(text = paste("posterior mean for Pr(R | D) = ",
                     round(digits=4, postmean.p[1])),
        side=3, cex=2,
     col = ColorForPosterior)
  abline(v=postmean.p[1],col=ColorForPosterior, lty=2, lwd=2)
}


