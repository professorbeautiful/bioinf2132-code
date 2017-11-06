"make.dataset.contam" <- 
function(phi, n, savedseed, print=TRUE, plot=TRUE, xvalues=seq(-5,5,0.1))
{
	###   phi is a vector with elements
	#	p1, sig1, sig2, mu1, mu2
	# or
	#	p1, sig1, sig2, mu
	# or
	#	p1, sig, mu1, mu2
	"%includes%" = function(a,b) b %in% a
	p1 = phi["p1"]
	if(names(phi) %includes% "sig1"
	 & names(phi) %includes% "sig2") {
		sig1 = phi["sig1"]
		sig2 = phi["sig2"]
	}
	else {
		sig1 = sig2 = phi["sig"]
	}
	sig1 = phi["sig1"]
	sig2 = phi["sig2"]
	if(names(phi) %includes% "mu1"
	 & names(phi) %includes% "mu2") {
		mu1 = phi["mu1"]
		mu2 = phi["mu2"]
	}
	else {
		mu1 = mu2 = phi["mu"]
	}
	if(! missing(savedseed))
		.Random.seed <- savedseed
	else
		savedseed <- .Random.seed
	sig1sq <- sig1^2
	sig2sq <- sig2^2
	Uind <- runif(n)
	Jind <- rep(1, n)
	Jind[Uind >= p1] <- 2
	n1 <- sum(Jind == 1)
	n2 <- n - n1
	y <- rep(NA, n)
	y[Jind == 1] <- rnorm(n1, mu1, sig1)
	y[Jind == 2] <- rnorm(n2, mu2, sig2)

	if(print == TRUE) {
		print(sapply(split(y, Jind),summary))
		print(tapply(y, Jind, summary))
		cat (" Std Devs\n")
		print(tapply(y, Jind, sd))
		cat (" Proportion of NONoutliers = ", 
			mean(Jind==1),
			" compare with p=", phi["p1"], "\n")
	}
	if(plot == TRUE) {
		plot(xvalues, 
			dnorm(xvalues,0,phi["sig1"]), col="green", type="l")
		lines(xvalues, 	
			dnorm(xvalues,0,phi["sig2"]), col="red", type="l")
		lines(xvalues, 
			phi["p1"] * dnorm(xvalues,0,phi["sig1"])
			 + (1-phi["p1"]) * dnorm(xvalues,0,phi["sig2"]), type="l",
			lwd=2)
		points(y[Jind==1], rep(0,sum(Jind==1)), pch="1", col="green")
		points(y[Jind==2], rep(0.01,sum(Jind==2)), pch="2", col="red")
		lines(density(y), new=FALSE, col="blue")
	####   END OF PLOTTING  ############################
	}

	dtemp <- data.frame(y, Uind, Jind)
	attr(dtemp, "call") <- sys.call()
	attr(dtemp, "defaults") <- args(make.dataset.contam)
	attr(dtemp, "initial.seed") <- savedseed
	return(dtemp)
}

