####  2011-02-01
###  The binomial CDF and quantile functions
###  Harry and Tom play tennis ---- or else flipping a coin.

###  the binomial CDF
pbinom(0:12, 12, 0.5)
plot(stepfun(0:12, c(0, .Last.value)),	main="CDF = F(x), pbinom", ylab="F(x)", col="darkgreen")
pnbinom(0:12, 3, 0.5)
lines(stepfun(0:12, c(0, .Last.value)), col="red")
legend(locator(1), legend=c("binomial", "negbin"), col=c("darkgreen","red"), lty=c(1,1),
	text.col=c("darkgreen","red"))

###  a close-up of the binomial CDF.  Note how the CDFs cross.  The red line never reaches 1.
dev.new()
pbinom(0:12, 12, 0.5)
plot(stepfun(0:12, c(0, .Last.value)),	main="CDF = F(x), pbinom", ylab="F(x)", col="darkgreen",
	xlim=c(7,12), ylim=c(0.9,1))
pnbinom(0:12, 3, 0.5)
lines(stepfun(0:12, c(0, .Last.value)), col="red")
legend(locator(1), legend=c("binomial", "negbin"), col=c("darkgreen","red"), lty=c(1,1),
	text.col=c("darkgreen","red"))

###  1 minus the CDF, to get the P-values.
1 - pbinom(0:12, 12, 0.5)
plot(stepfun(0:12, c(0, .Last.value)),	main="CDF = F(x), pbinom", ylab="F(x)",
	xlim=c(7,12), ylim=c(0,0.1), col="darkgreen")
1 - pnbinom(0:12, 3, 0.5)
lines(stepfun(0:12, c(0, .Last.value)), col="red")
text.default(c(8,8), c(0.075, 0.035), c(0.075, 0.035), col=c("darkgreen","red"), pos=4)
legend(locator(1), legend=c("binomial", "negbin"), col=c("darkgreen","red"), lty=c(1,1),
	text.col=c("darkgreen","red"))
	
####  the binomial quantile function
dev.new(width=10, height=5)    ### new graphics device 
par(mfrow=c(1,2))

plot(stepfun(0:12, c(0,pbinom(0:12, 12, 0.4))), main="CDF function = F(q) = pbinom")
lines(0:12, pbinom(0:12, 12, 0.4), col="red")
ptemp<-seq(0,1,length=1000)
plot(ptemp, qbinom(ptemp, 12, 0.4), pch=".", main="Quantile function = F^(-1)(p) = qbinom")
notDup = ! duplicated(qbinom(ptemp, 12, 0.4))
points(pbinom(0:12, 12, 0.4), 0:12, col="red", cex=1.1, type="b")
### Notice that the two graphics are essentially mirror images.

####  quantile function
qbinom(
	pbinom(0:12, 12, 0.4),
	12, 0.4)
### This shows that pbinom=F(X) and qbinom=F^(-1)(X) are inverse functions.
qnbinom(
	pnbinom(0:12, 12, 0.4),
	12, 0.4)
########################



##  Bayesian analysis ######

dev.new(width=10, height=5)
par(mfrow=c(1,2))
theta.vector  =  seq(0,1,length=500)

###  (A)  flipping coins:   prior is heavily near 1/2.

## likelihood  (normalized)
plot(theta.vector, normalize(dbinom(9, 12, theta.vector)) * 500, col="darkgreen", type="l",
	ylim=c(0,10), lwd=3)
## prior
a0 = 50; b0 = 50
lines(theta.vector, dbeta(theta.vector, a0, b0), lwd=3, lty=2)
## posterior
lines(theta.vector, dbeta(theta.vector, a0+9, b0+3), lwd=3, lty=3, col="darkgreen")
title("highly informative prior a0=50 b0=50\n(coin flipping)")
legend(0, 8, legend=c("likelihood", "prior", "posterior"),
	text.col=c("darkgreen", "black", "darkgreen"),
	col=c("darkgreen", "black", "darkgreen"),
	lty=c(1,2,3),
	lwd=c(3,3,3))
	
	
###  (B)  tennis matches:   prior is nearly uninformative.

## likelihood  (normalized)
plot(theta.vector, normalize(dbinom(9, 12, theta.vector)) * 500, col="darkgreen", type="l",
	ylim=c(0,10), lwd=3)
## prior
a0 = 2; b0 = 2
lines(theta.vector, dbeta(theta.vector, a0, b0), lwd=3, lty=2)
## posterior
lines(theta.vector, dbeta(theta.vector, a0+9, b0+3) 
, lwd=3, lty=3, col="darkgreen")
title("almost uninformative prior a0=2 b0=2\n(tennis matches)")
legend(0, 8, legend=c("likelihood", "prior", "posterior"),
	text.col=c("darkgreen", "black", "darkgreen"),
	col=c("darkgreen", "black", "darkgreen"),
	lty=c(1,2,3),
	lwd=c(3,3,3))

######   EXERCISE:  re-run these plots (A) and (B), but using the negative binomial instead.
######  Compare results to the binomial.   Discuss.

