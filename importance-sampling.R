### Simple example of importance sampling.
###  g = triangular distribution from -1 to 1.
###  h  = x^2   
###  Our goal is E(h(x)) - INT (g(x)h(x))
###  I  = Unif on (-1,1)
par(mfrow=c(1,1))
g = dtriang = function(x) pmax(0, 1-abs(x))
plot(x<-seq(-1,1,.1), g(x), pch="g", type="b", col="red")
h = square = function(x) x^2
lines(x, g(x)*h(x))
text(x = c(-0.5, 0.5), c(0.06, 0.06), "g*h")

####  Unweighted sampling
samplesize = 1000
rtriang = function(samplesize)
  (1-sqrt(runif(samplesize))) * (2*(rbinom(samplesize,1,1/2))-1)
lines(density(rtriang(samplesize)),  
      col="red", lty=2, lwd=2)

samplesize = 10;    nreps = 1000
answers.unweighted = sapply(1:nreps,function(i,samplesize)
  mean(h(rtriang(samplesize))), samplesize=samplesize)
est.unweighted = mean(answers.unweighted)
sd.unweighted = sd(answers.unweighted)
###  sd, not stdev, in R

# Now we introduce a density that mimics g*h fairly well.
mean.I = 0.6
sd.I = 0.3
dI = function(x) dnorm(x, mean.I, sd.I)/2 + 
                dnorm(x, -mean.I, sd.I)/2
lines(x<-seq(-1,1,.01), dI(x), pch="I", col="blue")
rI = function(n, m=mean.I, sd=sd.I) {
  rnorm(n, m*(2*rbinom(n,1,1/2)-1), sd)
}
lines(density(rI(samplesize)), new=F, col="blue", lty=2)

### Weighted average this time.
### We normalize by the sample size.

answers.weighted = sapply(1:nreps,
                          function(i,samplesize) {
                            rpoints = rI(samplesize) ### Sampling from I
                            weights = g(rpoints)/dI(rpoints)
                            return(sum(weights*h(rpoints))/samplesize)
                          }, samplesize=samplesize)
est.weighted = mean(answers.weighted)
sd.weighted = sd(answers.weighted)
data.frame(row.names = c("unweighted", "weighted"),
           est=c(est.unweighted, est.weighted),
           sd=c(sd.unweighted, sd.weighted) )

plot(density(answers.unweighted), col="red", ylim=c(0,200))
rug(answers.unweighted, col='red')
lines(density(answers.weighted), col="blue")
rug(answers.weighted, col="#0000ff22", ticksize = 0.05)
### 22 is opacity.

cat("Relative Efficiency = ", (sd.unweighted/sd.weighted)^2, "\n")

### Note that E(sum(weights) | I) = 1,
###  but in general, sum(weights) does not = 1.
### A true weighted average would normalize by sum(weights).

answers.weighted2 = sapply(1:nreps,
                           function(i,samplesize) {
                             rpoints = rI(samplesize) ### Sampling from I
                             weights = g(rpoints)/dI(rpoints)
                             return(sum(weights*h(rpoints))/sum(weights))
                           }, samplesize=samplesize)
est.weighted2 = mean(answers.weighted2)
sd.weighted2 = sd(answers.weighted2)
data.frame(row.names = c("unweighted", "weighted; norm by samplesize", 'weighted; norm by sum(w)'),
           est=c(est.unweighted, est.weighted, est.weighted2),
           sd=c(sd.unweighted, sd.weighted, sd.weighted2) )
hist(answers.weighted2, add=T, xlim=c(0.1, 0.2), breaks = 15, col="green")

#' imp.samp
#' 
#'   This approach shows how to write a general function that takes functions as args.
#' @param  h.fun is the function of X to take the expected value of.
#' @param  g.fun is proportional to the density of X. Need not be normalized.
#' @param  dI.fun is the importance sampling density to use instead.
#' @param  rI.fun produces the random values of X according to dI.fun
#' @param  ... is extra args to g.fun.
imp.samp <- function(h.fun, g.fun, dI.fun, rI.fun, n, divide.by.n=T, ...) {
  
  x <- rI.fun(n)
  if (is.character(h.fun))
    h.values <- do.call(h.fun, list(x))
  else
    h.values <- do.call("h.fun", list(x))
  g.values <- do.call("g.fun", list(x, ...))
  dI.values <- do.call("dI.fun", list(x))
  weights <- g.values/dI.values
  if(divide.by.n==F)
    result <- sum(weights*h.values)/sum(weights)
  else
    result <- sum(weights*h.values)/n
  if(divide.by.n==F)
    stderr <- sqrt(sum( ((h.values-result)*weights)^2 ))/sum(weights)
  else
    stderr <- sd(weights*h.values)/sqrt(n)
  #cat(" Using call ", sys.call())
  cat("  \nResult is ", format(result, digits=5), "  (", format(stderr,digits=5), ")\n")
  return(invisible(list(result=result, stderr=stderr, x=x, h.values=h.values, weights=weights)))
}

IDENTITY <- function(x)x
ONE <- function(x)rep(1,length(x))
CONSTANT = function(x, Constant=1)rep(Constant,length(x))
## In this example, we want the area under the triangle:  
imp.samp.1 <- imp.samp(h.fun = "IDENTITY" , g.fun = ONE , 
                       dI.fun = function(x)dbeta(x,1,1), rI.fun = function(n)rbeta(n,1,1), 
                       n = 100)
imp.samp.2 <- imp.samp(h.fun = IDENTITY , g.fun = ONE , 
                       dI.fun = function(x)dbeta(x,1,2), rI.fun = function(n)rbeta(n,1,2), 
                       n = 100)
### Perfect match:  I = g * h.
imp.samp.3 <- imp.samp(h.fun = "IDENTITY" , g.fun = ONE , 
                       dI.fun = function(x)dbeta(x,2,1), rI.fun = function(n)rbeta(n,2,1), 
                       n = 100)
## A true weighted sum is not as good.
imp.samp.3.w <- imp.samp(h.fun = "IDENTITY" , g.fun = ONE , 
                       dI.fun = function(x)dbeta(x,2,1), rI.fun = function(n)rbeta(n,2,1), 
                       n = 100, divide.by.n=F)
sum(imp.samp.3$weights)
sum(imp.samp.3.w$weights)
imp.samp.4 <- imp.samp(h.fun = "IDENTITY" , g.fun = CONSTANT , 
                       dI.fun = function(x)dbeta(x,2,1), rI.fun = function(n)rbeta(n,2,1), 
                       n = 100, Constant=2)
imp.samp.5 <- imp.samp(h.fun = "IDENTITY" , g.fun = CONSTANT , 
                       dI.fun = function(x)dbeta(x,2,1), rI.fun = function(n)rbeta(n,2,1), 
                       n = 100, Constant=2,  divide.by.n=F)
###  So, if g is not normalized, we'd better divide by the sum of the weights, not by n.