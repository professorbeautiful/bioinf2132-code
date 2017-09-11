#### Question 3.3
genlink.data.big   <- c(125,18,20,34)
genlink.data.small <- c(14,0,1,5)

# use importance sampling from normal distribution to estimate posterior mean, sd
imp.samp <- function(x,iter=10000,I.random, I.dens, ...){
# 	print(sys.call())
# 	print(length(sys.call()))
# 	print(names(sys.call()))
# 	for(i in 1:length(sys.call()))
# 		print(sys.call()[[i]])
	theta.star <- I.random(iter,...) #sample from normal dist
	I.theta <- I.dens(theta.star,...)  #density from normal dist
	g.theta <- ((2+theta.star)/4)^x[1]*((1-theta.star)/4)^(x[2]+x[3])*(theta.star/4)^x[4] #likelihood
	g.theta <- g.theta*(theta.star>=0 & theta.star<=1) #weight of 0 if theta.star <0 or >1
	weights <- g.theta/I.theta
	imp.samp.mean <- sum(theta.star*weights)/sum(weights) # Posterior E[theta]
	imp.samp.x2 <- sum(theta.star^2*weights)/sum(weights) # Posterior E[theta^2]
	imp.samp.sd <- sqrt(imp.samp.x2-imp.samp.mean^2) # Posterior sd[theta]
	plot(sort(weights))   ##### Quantiles of the weights.
	return(list(post.mean=imp.samp.mean,post.sd=imp.samp.sd,sd.weights=sd(weights)))
}

###################################################################################################
#3a) Posterior mean+sd using importance sampling + normal importance function (data set 1)
	imp.samp(c(125,18,20,34), iter=10000, rnorm, dnorm, mean=0.5,sd=0.25)
#3b) Posterior mean+sd using importance sampling + normal importance function (data set 2)
	imp.samp(c(14,0,1,5), iter=10000, rnorm, dnorm, mean=0.5,sd=0.25)

#3c) Posterior mean+sd using importance sampling + beta importance function
	imp.samp(c(125,18,20,34), iter=10000, rbeta, dbeta, shape1=0.5, shape2=0.5)
	imp.samp(c(14,0,1,5), iter=10000, rbeta, dbeta, shape1=0.5, shape2=0.5)
###################################################################################################


# likelihood function
	lik.fn <- function(theta,x) ((2+theta)/4)^x[1]*((1-theta)/4)^(x[2]+x[3])*(theta/4)^x[4]

# importance sampling estimate (from normal) for denominator of normalized likelihood 
imp.samp.sum.norm <- function(x,iter=10000,mean=0.5,sd=0.25){
	theta.star <- rnorm(iter,mean,sd)
	I.theta <- dnorm(theta.star,mean,sd)
	g.theta <-  ((2+theta.star)/4)^x[1]*((1-theta.star)/4)^(x[2]+x[3])*(theta.star/4)^x[4]
	weights <- 1/I.theta*(theta.star>=0 & theta.star<=1) #weight of 0 if theta.star <0 or >1
	imp.samp.sum <- sum(g.theta*weights)/sum(weights)
	imp.samp.sum
}

#################################################################################################
#3d) Normalized likelihood for first data set using importance sampling
data.1 <- (c(125,18,20,34))
norm.lik.1 <- lik.fn(theta.plot,data.1)/imp.samp.sum.norm(data.1)
plot(theta.plot,norm.lik.1,type="l")

# "true" normalized likelihood; denominator uses integral
real.norm.lik.1 <- lik.fn(theta.plot,data.1)/integrate(lik.fn,0,1,x=data.1)$value 
lines(theta.plot,real.norm.lik.1,lty=2)
legend(0,8,c("importance sampling","true"),lty=c(1,2))

#3e) Normalized likelihood for second data set using importance sampling
data.2 <- c(14,0,1,5)
norm.lik.2 <- lik.fn(theta.plot,data.2)/imp.samp.sum.norm(data.2)
plot(theta.plot,norm.lik.2,type="l")

real.norm.lik.2 <- lik.fn(theta.plot,data.2)/integrate(lik.fn,0,1,x=data.2)$value
lines(theta.plot,real.norm.lik.2,lty=2)
legend(0,4,c("importance sampling","true"),lty=c(1,2))
