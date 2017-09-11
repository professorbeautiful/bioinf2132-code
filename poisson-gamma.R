plot(x<-(221-30):(221+30), dgamma(x,221,1), type="l", ylim=c(0,0.028))
abline(v=221, lty=2)
toleranceInterval = qgamma(c(0.05,0.95), 221, 1)
lines(rep(toleranceInterval[1], each=2),
	c(0, dgamma(toleranceInterval[1],221,1)), lty=2)
lines(rep(toleranceInterval[2], each=2),
	c(0, dgamma(toleranceInterval[2],221,1)), lty=2)
lines(x, dnorm(x, 221, sqrt(221)), lty=3)
