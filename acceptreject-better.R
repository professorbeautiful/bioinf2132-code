### acceptance-rejection sampling example
### Modified 2014-03-19

a <- 2
b <- 4
par(mfrow=c(1,1))
p.seq <- seq (0.01, 0.99, 0.01)
plot(p.seq, dbeta(p.seq, a, b), type="l", lwd=2, ylim=c(0,5))
#M <- max(dbeta(p.seq, a, b))

n <- 5000

theMean = a/(a+b)
theVar = a*b/(a+b)/(a+b+1)
samplingDensity = dnorm(p.seq, theMean, sqrt(theVar))
lines(p.seq, samplingDensity, col="orange", lty=2, lwd=3)

M = max(dbeta(p.seq, a, b)/samplingDensity)  * 1.1
#M=6
lines(p.seq, samplingDensity * M, col="orange", lwd=2)

theta <- rnorm(n,theMean, sqrt(theVar))
u <- runif(n)

heights = u*M*dnorm(theta, theMean, sqrt(theVar))
points (theta, heights, pch=".", cex=4, col="red")


accept <- (heights < dbeta(theta, a, b))

n.accept <- sum(accept)
points (theta[!accept], heights[!accept], col="black", cex=4, pch=".")   
points (theta[ accept], heights[ accept], col="green", cex=4, pch=".")   
cat("Proportion accepted is ", n.accept / n, "\n")

qqplot( theta[accept],
  qbeta(seq(0,1,length=n.accept), a, b))
abline(a=0,b=1, col="green", lwd=3)
