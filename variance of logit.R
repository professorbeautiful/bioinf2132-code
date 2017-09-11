#quartz(width=8, height=4)
par(mfrow=c(1,2))
#### Sampling model:  binomial
p <- 1/4
n <- 1000
size <- 100
k <- rbinom(n, size, p)
logit = function(P) log(P/(1-P))
table(theLogits <- logit(k/size))
notExtreme <- k!=size & k!=0
meanHat <- logit(mean(k)/size)
varHat <- 1/mean(k) + 1/(size-mean(k))
qqnorm(
	(theLogits[notExtreme]  -  meanHat) 
		/ sqrt(varHat),
	main=paste("Binomial\n",
			"size=", size, "  p=", p))
abline(a=0, b=1)

#### Sampling model:  independent Poissons
k1 <- rpois(n, size * p)
k2 <- rpois(n, size * (1-p))
table(theLogits <- logit(k1/(k1+k2)))
notExtreme <- k1!=0 & k2!=0
meanHat <- logit( mean(k1)/mean(k1+k2) )
varHat <- 1/mean(k1) + 1/mean(k2)
qqnorm(
	(theLogits[notExtreme]  -  meanHat) 
		/ sqrt(varHat),
	main=paste("Independent Poisson\n",
			"size=", size, "  p=", p))
abline(a=0,b=1)
