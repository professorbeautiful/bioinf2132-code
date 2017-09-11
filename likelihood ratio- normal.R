plot(theta<-seq(-4, 4,length=500), dnorm(theta, mean=-1, sd=1), type="l", col="red")
lines(theta, dnorm(theta, mean=+1, sd=1), type="l", col="blue")

likelihoodRatio = dnorm(theta, mean=-1, sd=1) / dnorm(theta, mean=+1, sd=1)

lines(theta, likelihoodRatio / max(likelihoodRatio) * par("usr")[4], col="green")

axis(4, labels = T)  ####  Error message: cannot coerce type 'closure' to vector of type 'double'

