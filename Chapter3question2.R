genlink.data.big <- c(125,18,20,34)
genlink.data.small <- c(14,0,1,5)

# uses Laplace approximation to estimate posterior mean+variance
genlink.lap <- function(the.data=genlink.data.big, plot.do=TRUE){
	# assess log-likelihood
	log.lik.neg <- function(theta,the.data)
		-(the.data[1]*log((2+theta)/4)+(the.data[2]+the.data[3])*log((1-theta)/4)+the.data[4]*log(theta/4))
	# assess log-likelihood with g(theta)=theta
	log.lik.star.neg <- function(theta,the.data)  log.lik.neg(theta, the.data) - log(theta)
	# assess log-likelihood with g(theta)=theta^2
	log.lik.star2.neg <- function(theta,the.data)   log.lik.neg(theta, the.data) - log(theta^2)

	# obtain MLEs for above log-likelihoods
	theta.hat <- optimize(log.lik.neg, 0:1, the.data=the.data)$minimum
	theta.star <- optimize(log.lik.star.neg, 0:1, the.data=the.data)$minimum
	theta2.star <- optimize(log.lik.star2.neg, 0:1, the.data=the.data)$minimum

	# obtain second derivatives for above log-likelihoods at MLEs
	sigma.fn.formula = ~ (-1)*(y1*log((2+theta)/4)+(y2+y3)*log((1-theta)/4)+y4*log(theta/4))
	"+.formula" = function(f1, f2) {
			if(class(f1)=="formula") f1=as.character(f1)[-1]
			if(class(f2)=="formula") f2=as.character(f2)[-1]
			formula(paste("~", f1, "+", f2, collapse=""))
	}
	sigma.fn <- deriv(sigma.fn.formula,
						"theta",
						hessian=TRUE,
						function.arg=function(theta,y1,y2,y3,y4) NULL)
	sigma.hessian <- function(theta,y1,y2,y3,y4) 
								attr(sigma.fn(theta,y1,y2,y3,y4), "hessian")
	sigma <- sqrt(1/sigma.hessian(theta.hat,the.data[1],the.data[2],the.data[3],the.data[4]))
	sigma <- c(sigma)  ### convert to a simple scalar

	sigma.star.fn.formula = "+"(sigma.fn.formula, "-log(theta)")
	sigma.star.fn <- deriv(sigma.star.fn.formula,
						"theta",
						hessian=TRUE,
						function.arg=function(theta,y1,y2,y3,y4) NULL)
	sigma.star.hessian <- function(theta,y1,y2,y3,y4) 
								attr(sigma.star.fn(theta,y1,y2,y3,y4), "hessian")
	sigma.star <- sqrt(1/sigma.star.hessian(theta.hat,the.data[1],the.data[2],the.data[3],the.data[4]))
	sigma.star <- c(sigma.star)  ### convert to a simple scalar

	sigma2.star.fn.formula = "+"(sigma.fn.formula, "-log(theta^2)")
	sigma2.star.fn <- deriv(sigma2.star.fn.formula,
						"theta",
						hessian=TRUE,
						function.arg=function(theta,y1,y2,y3,y4) NULL)
	sigma2.star.hessian <- function(theta,y1,y2,y3,y4) 
								attr(sigma2.star.fn(theta,y1,y2,y3,y4), "hessian")
	sigma2.star <- sqrt(1/sigma2.star.hessian(theta.hat,the.data[1],the.data[2],the.data[3],the.data[4]))
	sigma2.star <- c(sigma2.star)  ### convert to a simple scalar

	#Expected value of g(theta)=theta
	post.mean.lap <- (sigma.star/sigma)*
		(exp(-log.lik.star.neg(theta.star,the.data))/
		 exp(-log.lik.neg(theta.hat,the.data)))
	#Expected value of g(theta)=theta^2
	post.x2.lap <- (sigma2.star/sigma)*
		(exp(-log.lik.star2.neg(theta2.star,the.data))/
		 exp(-log.lik.neg(theta.hat,the.data)))
	post.sd.lap <- sqrt(post.x2.lap-post.mean.lap^2)
	if(plot.do) {
		npoints = 100
		theta.seq <- seq(0,1,length=npoints)
		logliks  <- -log.lik.neg(theta.seq, the.data)
		normalizer <- sum(exp(logliks)) / npoints
		options("device")[[1]]()  ## creates a new graph on any platform
		plot(theta.seq, (exp(logliks)/normalizer), type="l")
		lines(theta.seq, dnorm(theta.seq, mean=post.mean.lap, sd=post.sd.lap), col="red")
		options("device")[[1]]()  ## creates a new graph on any platform 
		plot(theta.seq, log(exp(logliks)/normalizer), type="l")
		lines(theta.seq, log(dnorm(theta.seq, mean=post.mean.lap, sd=post.sd.lap)), col="red")
	}
	return(c(post.mean=post.mean.lap,post.sd=post.sd.lap))
}

###################################################################################################
#2a) Posterior mean+sd for first data set
	genlink.lap(genlink.data.big)
#2b) Posterior mean+sd for second data set
	genlink.lap(genlink.data.small)
###################################################################################################
	⁃	
