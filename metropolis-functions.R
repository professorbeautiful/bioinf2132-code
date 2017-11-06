#### MCMC metropolis function

logit <- function(P) log(P) - log(1-P)
antilogit <- function(Z) 1 - 1/(1+exp(Z))

### Metropolis-Hastings algorithm
#do.call.with.matching.dots = function(f, args, ...) do.call(f, args, ...)
metrop <- function (start, jumpsig,
                    dens=function(p)odl.rao(y.rao,p),
                    rtrans=function(p,jumpsig)antilogit(logit(p)+rnorm(1,0,jumpsig)),   ### randomly chosen proposal
                    Jtrans=function(p1,p2,jumpsig)  ### prob density of proposal
                    {
                      jacobean <- p2*(1-p2); #1/(dlogit(p2)/dp2);
                      return(dnorm(logit(p2)-logit(p1),sd=jumpsig) * jacobean/jumpsig)
                    },
                    chainlength=100,
                    nchains=ifelse(is.list(start), length(start), 1),
                    do.print=FALSE, do.plot=FALSE, ...)
{
  print(is.list(start))
  if(is.list(start) | (nchains > 1)) {
    nchains <- ifelse(is.list(start), length(start), nchains)
    cat("length(start) = ", length(start), "\n")
    cat("nchains = ", nchains, "\n")
    # We assume each entry in start is a different starting value.
    ## TODO: extract the dots args, check if they should be distributed or repeated to the
    #  calls to metrop, then send the appropriate args with do.call().
    chains <- list()
    for (j in 1:length(start)) {
      chains[[j]] <- metrop(start[[j]], jumpsig=jumpsig, dens=dens,rtrans=rtrans,
                            Jtrans=Jtrans, chainlength=chainlength,
                            do.print=do.print, do.plot=do.plot,
                            ...)
    }
    oldClass(chains) <- "metrop.output"
    return(chains)
  }
  jumpcount <- 0
  values <- proposals <- matrix(NA, nrow=chainlength, ncol=length(start))
  densities <- rep(NA,chainlength)
  jumped <- rep(0,chainlength)
  values[1, ] <- current.value <- start
  densities[1] <- current.dens <- dens(current.value)
  for (i in 2:chainlength) {
    proposal.value <- proposals[i, ] <-
      rtrans(p=current.value, jumpsig=jumpsig)
    proposal.dens <- dens(proposal.value)
    Jtrans.jump <- Jtrans(p1=current.value,p2=proposal.value,jumpsig=jumpsig)
    Jtrans.back <- Jtrans(p1=proposal.value,p2=current.value,jumpsig=jumpsig)
    w.p1p2 <- (proposal.dens * Jtrans.back) /
      (current.dens * Jtrans.jump)
    pr.jump <- min (1, w.p1p2)
    tryResult = try({
      if (runif(1)<pr.jump) {
        current.value <- proposal.value
        current.dens <- proposal.dens
        jumpcount <- jumpcount+1
        jumped[i] <- 1
      }
    })
    if('try-error' == class(tryResult)) browser()
    values[i, ] <- current.value
    densities[i] <- current.dens
    if(do.print)
      print(c(proposal.value,proposal.dens,
              w.p1p2,pr.jump,current.value))
  }
  cat("Jump acceptance proportion is ", jumpcount/chainlength, "\n")
  #cat("Lag 1 correlation is ", acf(ts(1:5),1), "\n")
  chain <- list(values=values,densities=densities,proposals=proposals,
                jumped=jumped,start=start,jumpsig=jumpsig,chainlength=chainlength,call=sys.call())
  oldClass(chain) <- "metrop.output"
  attr(chain, "start") = start
  attr(chain, "jumpsig") = jumpsig
  if(do.plot)  plot(chain)
  return(invisible(chain))
}

plot.metrop.output <- function(chain, whichdim=1, burnstop=1, newplot=TRUE, ...) {
  if(inherits(chain[[1]],"metrop.output")) {
    ## arg is a list of chains
    for (i in 1:length(chain)) {
      plot(chain[[i]] [, whichdim], burnstop=burnstop, newplot=TRUE, ...)  ##not great-- temporary measure.
    }
    return(invisible(NULL))
  }
  #if(newplot) quartz(width=12, height=6)
  par(mfrow=c(1,2))
  plot.metrop.output.1(chain, burnstop=1, newplot=TRUE, whichdim=whichdim, ...)
  plot.metrop.output.2(chain, burnstop=1, newplot=TRUE, whichdim=whichdim, ...)
  par(mfrow=c(1,1))
}

plot.metrop.output.1 <- function(chain, burnstop=1, newplot=TRUE, whichdim=1, ...) {
  plot.default(burnstop:chain$chainlength, 
               chain$values[burnstop:chain$chainlength, whichdim], ...)
  title("chains")
}

plot.metrop.output.2 <-
plot.metrop.output.compare.to.true.posterior <- 
  function(chain, burnstop=1, newplot=TRUE, whichdim=1, ...) {
   # browser()
    argrange=odl.rao.args; posterior=odl.rao.pdf
    plot.default(argrange, posterior, col=3, type="l", ...,
                 xlab='parameter', ylab='posterior density',
                 main="true and estimated \nposterior density")
    lines(density(chain$values[burnstop:chain$chainlength, whichdim]), col="red")
    rug(chain$values[burnstop:chain$chainlength, whichdim], col="red")
}
