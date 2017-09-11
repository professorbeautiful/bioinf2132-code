# 
# odl.rao <- function(y,p)
# 	{ (1/2+p/4)^y[1] * ((1-p)/4)^(y[2]+y[3]) * (p/4)^y[4] }
# y.rao <- c(125,18,20,34)  #Original data.
# 
# ####  Compare with true posterior  (flat prior)
# Npoints <- 600
# odl.rao.args <- seq(0.01,0.99,length=Npoints)
# odl.rao.values <- odl.rao(y.rao, p=odl.rao.args)
# odl.rao.pdf <- odl.rao.values/sum(odl.rao.values) * Npoints
# odl.rao.cdf <- cumsum(odl.rao.pdf) 
# 
# logit <- function(P) log(P) - log(1-P)
# antilogit <- function(Z) 1 - 1/(1+exp(Z))
# 
#   ### Metropolis-Hastings algorithm
# #do.call.with.matching.dots = function(f, args, ...) do.call(f, args, ...)
# metrop <- function (start, jumpsig,
# 		dens=function(p)odl.rao(y.rao,p), 
# 		rtrans=function(p,jumpsig)antilogit(logit(p)+rnorm(1,0,jumpsig)),   ### randomly chosen proposal
# 		Jtrans=function(p1,p2,jumpsig)  ### prob density of proposal
# 			{
#         jacobean <- p2*(1-p2); #1/(dlogit(p2)/dp2);
#         return(dnorm(logit(p2)-logit(p1),sd=jumpsig) * jacobean/jumpsig)
# 		  },
#     chainlength=100, 
# 		nchains=ifelse(is.list(start), length(start), 1), 
# 		do.print=FALSE, do.plot=FALSE, ...) 
# {
# 	print(is.list(start))
# 	if(is.list(start) | (nchains > 1)) {
# 		nchains <- ifelse(is.list(start), length(start), nchains)
# 		cat("length(start) = ", length(start), "\n")
# 		cat("nchains = ", nchains, "\n")
# 		# We assume each entry in start is a different starting value.
# 		## TODO: extract the dots args, check if they should be distributed or repeated to the 
# 		#  calls to metrop, then send the appropriate args with do.call().
# 		chains <- list()
# 		for (j in 1:length(start)) {
# 			chains[[j]] <- metrop(start[[j]], jumpsig=jumpsig, dens=dens,rtrans=rtrans,
# 				Jtrans=Jtrans, chainlength=chainlength, 
# 				do.print=do.print, do.plot=do.plot, 
# 				...)
# 		}
# 		oldClass(chains) <- "metrop.output"
# 		return(chains)
# 	}
# 	jumpcount <- 0
# 	values <- rep(NA,chainlength)
# 	proposals <- rep(NA,chainlength)
# 	densities <- rep(NA,chainlength)
# 	jumped <- rep(0,chainlength)
# 	values[1] <- current.value <- start
# 	densities[1] <- current.dens <- dens(current.value)
# 	for (i in 2:chainlength) {
# 		proposal.value <- proposals[i] <- 
# 			rtrans(p=current.value, jumpsig=jumpsig)
# 		proposal.dens <- dens(proposal.value)
# 		Jtrans.jump <- Jtrans(p1=current.value,p2=proposal.value,jumpsig=jumpsig)
# 		Jtrans.back <- Jtrans(p1=proposal.value,p2=current.value,jumpsig=jumpsig)
# 		w.p1p2 <- (proposal.dens * Jtrans.back) /
# 					(current.dens * Jtrans.jump)
# 		pr.jump <- min (1, w.p1p2)
# 		if (runif(1)<pr.jump) {
# 			current.value <- proposal.value
# 			current.dens <- proposal.dens
# 			jumpcount <- jumpcount+1
# 			jumped[i] <- 1
# 		}
# 		values[i] <- current.value
# 		densities[i] <- current.dens
# 		if(do.print) 
# 			print(c(proposal.value,proposal.dens,
# 				w.p1p2,pr.jump,current.value))
# 	}
# 	cat("Jump acceptance proportion is ", jumpcount/chainlength, "\n")
# 	#cat("Lag 1 correlation is ", acf(ts(1:5),1), "\n")
# 	chain <- list(values=values,densities=densities,proposals=proposals,
# 		jumped=jumped,start=start,jumpsig=jumpsig,chainlength=chainlength,call=sys.call())
# 	oldClass(chain) <- "metrop.output"
# 	attr(chain, "start") = start
# 	attr(chain, "jumpsig") = jumpsig
# 	if(do.plot)  plot(chain)
# 	return(invisible(chain))
# }
# 
# plot.metrop.output <- function(chain, burnstop=1, newplot=TRUE, ...) {
#   if(inherits(chain[[1]],"metrop.output")) {
#   	## arg is a list of chains
# 		for (i in 1:length(chain)) {
# 			plot(chain[[i]], burnstop=burnstop, newplot=TRUE, ...)  ##not great-- temporary measure.
# 		}
# 		return(invisible(NULL))
# 	}
# 	#if(newplot) quartz(width=12, height=6)
#   par(mfrow=c(1,2))
#   plot.metrop.output.1(chain, burnstop=1, newplot=TRUE, ...) 
#   plot.metrop.output.2(chain, burnstop=1, newplot=TRUE, ...) 
#   par(mfrow=c(1,1))
# }
# 
# plot.metrop.output.1 <- function(chain, burnstop=1, newplot=TRUE, ...) {
# 	plot.default(burnstop:chain$chainlength, chain$values[burnstop:chain$chainlength], ...)
# 	title("chains")
# }
# plot.metrop.output.2 <- function(chain, burnstop=1, newplot=TRUE, ...) {
#   plot.default(odl.rao.args,odl.rao.pdf,col=3, type="l", ...,
#                main="true and estimated \nposterior density")
# 	lines(density(chain$values[burnstop:chain$chainlength]), col="red")
# 	rug(chain$values[burnstop:chain$chainlength], col="red")
# } 
# 
# metrop.out.1 = metrop(start = .1,chainlength=200,jumpsig=1, do.print=FALSE, do.plot=TRUE); 
# lag.plot(metrop.out.1$values,
#          main="start=0.1, jumpsig=1", ask = F)
# 
# metrop.out.2 = metrop(start = .1, chainlength=50000,
#                       jumpsig=1, do.print=FALSE, do.plot=TRUE); 
# lag.plot(metrop.out.2$values, main="start=0.1, jumpsig=1"
#          , ask = F)
# 
# metrop.out.3 = metrop(start = .95, chainlength=50000,
#                       jumpsig=1, do.print=FALSE, do.plot=TRUE); 
# lag.plot(metrop.out.3$values, main="start=0.95, jumpsig=1")
# 
# metrop.out.4 = metrop(start = .1, chainlength=50000,
#                       jumpsig=5, do.print=FALSE, do.plot=TRUE); 
# ###  jumpsig is too big
# ### zoom in to see what's wrong.  too sticky.
# plot.metrop.output.1(metrop.out.4, xlim=c(5000,5500)) 
# lag.plot(metrop.out.4$values, main="start=0.1, jumpsig=5", ask = F)
# 
# metrop.out.5 = metrop(start = .1, chainlength=10000,
#                       jumpsig=0.01, do.print=FALSE, do.plot=TRUE); 
# ###  Too accepting!!!   Jumps too small.   Almost a random walk.
# lag.plot(metrop.out.5$values,   ### very high autocorrelation
#          main="start=0.1, jumpsig=0.01"
#          , ask = F)
# 
# #######  DENSITY ESTIMATION
# plotrange = 500:2000
# plotrange = 8000:10000
# 
# density.estimate = density(metrop.out.2$values[plotrange])
# plot(density.estimate, ylim=c(0,9),
# 		main=paste("density estimate and truth", 
# 			plotrange[1], rev(plotrange)[1]))
# lines(odl.rao.args, 
# 		odl.rao.values
# 		/sum(odl.rao.values[-1]*diff(odl.rao.args)), 
# 		col="red", lwd=2, lty=2)
# 
# ###Now, multiple chains at once:
# ### four chains with different initial points.
# ### This enables a diagnostic indicator: Convergence monitoring:

R.hat <- function(metrop.out, howfar) {
  chains = metrop.out$chains
  n <- chains[[1]]$chainlength  
  ## Assumes all chains have same length.
  if(missing(howfar))
    howfar <- seq(100, n, by=100)
  matrix.of.values <- t(plyr::laply(metrop.out, '[[', "values"))
  R.hat.out <- rep(NA,length(howfar))
  for (i in 1:length(howfar)) {
    howfarthistime <- howfar[i]
    Between <- var(apply(matrix.of.values[1:howfarthistime,], 2, mean)) * howfarthistime 
    Within <- mean(apply(matrix.of.values[1:howfarthistime,], 2, var))
    shrinker <- 1/howfarthistime 
    var.hat <- Within * (1-shrinker) + Between * shrinker
    R.hat.out[i] <- R.hat <- sqrt(var.hat/Within)
    #cat(" For chains of length ", howfarthistime, "  R.hat is ", R.hat, "\n")
  }
  names(R.hat.out) = as.character(howfar)
  unclass(R.hat.out)
}
testConvergence = function(start, jumpsig, 
                           chainlength=2000, nchains=4, 
                           howfar) {
  if(missing(start))
    start <- as.list(c(.1,.2,.9,.99))
  if(missing(howfar))
    howfar <- seq(100, chainlength, by=100)
  metrop.out <- metrop(start, jumpsig=jumpsig,
                       chainlength=chainlength,
                       do.plot=F,do.print=F, nchains=nchains)
  plot(metrop.out) 
  #require(plyr)
  values <- plyr::laply(metrop.out, '[[', "values")
  plot(rep(1:chainlength, nchains), c(values), pch="")
  for(chain in 1:nchains) points(1:chainlength, 
                           values[ chain, ], pch=".", cex=2, col=chain+1)
  rowMeans(values)
  plot(as.numeric(names(R.hat(metrop.out, howfar=howfar))), 
       R.hat(metrop.out, howfar=howfar),
       xlab = "chainlength", ylab="R.hat",
       main = paste0("jumpsig = ", jumpsig) )
}

testConvergence(jumpsig=5)
testConvergence(jumpsig=0.01)
