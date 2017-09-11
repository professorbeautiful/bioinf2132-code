###### Simulate the t distribution by monte carlo "method of composition"

###### R features & methods used:
# Creating empty placeholders for loops
# apply
# dt, dnorm, rgamma
# Graphics:  rug, lines.density, legend
# cq()  from mvbutils.

df <- 4
xRange <- seq(-3, 3, by=.1)
nsigma <- 500
plot(xRange, dt(xRange,df), type="l")
lines(xRange, dnorm(xRange), col=4, lty=2)
###Placeholders:
densities <- matrix(NA, nrow=nsigma, ncol=length(xRange))
sigsq <- vector ("numeric", length=nsigma)
xValues <- vector ("numeric", length=nsigma)
for ( isigma in 1:nsigma) {
	sigsq[isigma] <- 2/rgamma(1, df/2) 
		 ####  Inverse chisq random number (marginal).
	densities[isigma,] <- dnorm(xRange, sd=sqrt(sigsq[isigma])) 		 ### Conditional density [Yma|sig].
	xValues[isigma] <- rnorm(1,mean=0, sd=sqrt(sigsq[isigma]))
	    ### Conditional random number
}
marginalDensity <- apply(densities,2,mean)
lines(xRange, marginalDensity, col=2, lwd=2, lty=2)
rug(xValues, col="green")
lines(density(xValues), col="green")
legend(locator(1), legend=c("dt", "dnorm", "mean of dnorms", "nonparametric\nfrom sample"),
	col=cq(black,blue,red,green), lty=c(1,2,2,1))

#######  Demonstrating the use of quantile-quantile plots:
probs = seq(0, 1, length=length(xValues)+2)[-length(xValues)][-1]
qqplot(xValues, qt(probs, df=df))
abline(a=0, b=1, col="red")