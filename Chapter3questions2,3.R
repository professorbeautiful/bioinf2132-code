genlink.data.big <- c(125,18,20,34)
genlink.data.small <- c(14,0,1,5)

# uses Laplace approximation to estimate posterior mean+variance
genlink.lap <- function(the.data=genlink.data.big){
	# assess log-likelihood
	log.lik.neg <- function(theta,the.data)
		-(the.data[1]*log((2+theta)/4)+(the.data[2]+the.data[3])*log((1-theta)/4)+the.data[4]*log(theta/4))
	# assess log-likelihood with g(theta)=theta
	log.lik.star.neg <- function(theta,the.data)
		-(the.data[1]*log((2+theta)/4)+(the.data[2]+the.data[3])*log((1-theta)/4)+the.data[4]*log(theta/4)+log(theta))
	# assess log-likelihood with g(theta)=theta^2
	log.lik.star2.neg <- function(theta,the.data)
		-(the.data[1]*log((2+theta)/4)+(the.data[2]+the.data[3])*log((1-theta)/4)+the.data[4]*log(theta/4)+log(theta^2))

	# obtain MLEs for above log-likelihoods
	theta.hat <- optimize(log.lik.neg, 0:1, the.data=the.data)$minimum
	theta.star <- optimize(log.lik.star.neg, 0:1, the.data=the.data)$minimum
	theta2.star <- optimize(log.lik.star2.neg, 0:1, the.data=the.data)$minimum

	# obtain second derivatives for above log-likelihoods at MLEs
	sigma.fn <- deriv(~-(x1*log((2+theta)/4)+(x2+x3)*log((1-theta)/4)+x4*log(theta/4))/4,
						"theta",
						hessian=TRUE,
						function(theta,x1,x2,x3,x4) NULL)
	sigma.hessian <- function(theta,x1,x2,x3,x4) attr(sigma.fn(theta,theta,x1,x2,x3,x4),									"hessian")
	sigma <- sqrt(1/sigma.hessian(theta.hat,the.data[1],the.data[2],the.data[3],the.data[4]))

	sigma.star.fn <- deriv(~-(x1*log((2+theta)/4)+(x2+x3)*log((1-theta)/4)+x4*log(theta/4)+log(theta))/4,"theta",hessian=TRUE,function(theta,x1,x2,x3,x4) NULL)
	sigma.star <- sqrt(1/attr(sigma.star.fn(theta.star,the.data[1],the.data[2],the.data[3],the.data[4]),"hessian")[1])

	sigma2.star.fn <- deriv(~-(x1*log((2+theta)/4)+(x2+x3)*log((1-theta)/4)+x4*log(theta/4)+log(theta^2))/4,"theta",hessian=TRUE,function(theta,x1,x2,x3,x4) NULL)
	sigma2.star <- sqrt(1/attr(sigma2.star.fn(theta2.star,the.data[1],the.data[2],the.data[3],the.data[4]),"hessian")[1])

	#Expected value of g(theta)=theta
	post.mean.lap <- (sigma.star/sigma)*(exp(-log.lik.star.neg(theta.star,the.data))/exp(-log.lik.neg(theta.hat,the.data)))
	#Expected value of g(theta)=theta^2
	post.x2.lap <- (sigma2.star/sigma)*(exp(-log.lik.star2.neg(theta2.star,the.data))/exp(-log.lik.neg(theta.hat,the.data)))
	post.sd.lap <- sqrt(post.x2.lap-post.mean.lap^2)

	return(list(post.mean=post.mean.lap,post.sd=post.sd.lap))
}

###################################################################################################
#2a) Posterior mean+sd for first data set
	genlink.lap(genlink.data.big)
#2b) Posterior mean+sd for second data set
	genlink.lap(genlink.data.small)
###################################################################################################

# use importance sampling from normal distribution to estimate posterior mean, sd
imp.samp <- function(x,iter=10000,I.random, I.dens, ...){
	sum.x <- sum(x)
	theta.star <- I.random(iter,...) #sample from normal dist
	I.theta <- I.dens(theta.star,...)  #density from normal dist
	g.theta <- ((2+theta.star)/4)^x[1]*((1-theta.star)/4)^(x[2]+x[3])*(theta.star/4)^x[4] #likelihood
	g.theta <- g.theta*(theta.star>=0 & theta.star<=1) #weight of 0 if theta.star <0 or >1
	weights <- g.theta/I.theta
	imp.samp.mean <- sum(theta.star*weights)/sum(weights) # Posterior E[theta]
	imp.samp.x2 <- sum(theta.star^2*weights)/sum(weights) # Posterior E[theta^2]
	imp.samp.sd <- sqrt(imp.samp.x2-imp.samp.mean^2) # Posterior sd[theta]
	options("device")[[1]]()  ## creates a new graph window, regardless of platform.
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
sum.x <- sum(x)
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
