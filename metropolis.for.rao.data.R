source("metropolis-functions.R")

odl.rao <- function(y,p)
	{ (1/2+p/4)^y[1] * ((1-p)/4)^(y[2]+y[3]) * (p/4)^y[4] }
y.rao <- c(125,18,20,34)  #Original data.

####  Compare with true posterior  (flat prior)
Npoints <- 600
odl.rao.args <- seq(0.01,0.99,length=Npoints)
odl.rao.values <- odl.rao(y.rao, p=odl.rao.args)
odl.rao.pdf <- odl.rao.values/sum(odl.rao.values) * Npoints
odl.rao.cdf <- cumsum(odl.rao.pdf)

metrop.out.1 = metrop(start = .1,chainlength=200,jumpsig=1, 
                      do.print=FALSE, do.plot=T);
lag.plot(metrop.out.1$values,
         main="start=0.1, jumpsig=1", ask = F)
str(metrop.out.1)

metrop.out.2 = metrop(start = .1, chainlength=50000,
                      jumpsig=1, do.print=FALSE, do.plot=TRUE);
lag.plot(metrop.out.2$values, main="start=0.1, jumpsig=1"
         , ask = F)

metrop.out.3 = metrop(start = .95, chainlength=50000,
                      jumpsig=1, do.print=FALSE, do.plot=TRUE);
lag.plot(metrop.out.3$values, main="start=0.95, jumpsig=1")

metrop.out.4 = metrop(start = .1, chainlength=50000,
                      jumpsig=5, do.print=FALSE, do.plot=TRUE);
###  jumpsig is too big
### zoom in to see what's wrong.  too sticky.
plot.metrop.output.1(metrop.out.4, xlim=c(5000,5500))
lag.plot(metrop.out.4$values, main="start=0.1, jumpsig=5", ask = F)

metrop.out.5 = metrop(start = .1, chainlength=10000,
                      jumpsig=0.01, do.print=FALSE, do.plot=TRUE);
###  Too accepting!!!   Jumps too small.   Almost a random walk.
lag.plot(metrop.out.5$values,   ### very high autocorrelation
         main="start=0.1, jumpsig=0.01"
         , ask = F)

#######  DENSITY ESTIMATION
plotrange = 500:2000
plotrange = 8000:10000

density.estimate = density(metrop.out.2$values[plotrange])
plot(density.estimate, ylim=c(0,9),
		main=paste("density estimate and truth",
			plotrange[1], rev(plotrange)[1]))
lines(odl.rao.args,
		odl.rao.values
		/sum(odl.rao.values[-1]*diff(odl.rao.args)),
		col="red", lwd=2, lty=2)

###Now, multiple chains at once:
### four chains with different initial points.
### This enables a diagnostic indicator: Convergence monitoring:

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
  R.hat.output = R.hat(metrop.out, howfar=howfar)
  plot(as.numeric(names(R.hat.output)), 
       R.hat.output,
       xlab = "chainlength", ylab="R.hat",
       main = paste0("jumpsig = ", jumpsig) )
  return(R.hat.output)
}

testConvergence(jumpsig=5)  ## convergence looks good
testConvergence(jumpsig=0.01)  ## convergence looks bad
