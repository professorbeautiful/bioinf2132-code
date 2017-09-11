# Data
a=9;b=3
m = a/(a+b)

# future sample size
n=10;

## binomial at the mean a/(a+b)
## "estimate and plug in"
plot(0:n, dbinom(0:n, n, m), pch="F", type="b",
	ylab="probability", xlab="#heads", col="red")


### Compare to the beta-binomial
Beta = function(a, b) 
	exp(lgamma(a) + lgamma(b) - lgamma(a+b))
dbetabinomial=function(x,a,b,n){
	choose(n,x) * Beta(a+x,b+n-x) / Beta(a,b)
}

n=10;a=9+1;b=3+1
lines(0:n, dbetabinomial(0:n,a,b,n), type="b",
	col="green", lty=1, pch="B")

legend(0.1, .25, legend=
		c("\"frequentist\": binomial",
		 "Bayesian: beta-binomial"), 
		pch=c("F","B"), 
		col=c("red","green"))
title("Predictive distribution for #heads out of 10")
mtext("prior data:  9 heads, 3 tails",cex=1.5)
pngSave("betabinomial.png")