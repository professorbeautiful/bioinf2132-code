# uses Laplace approximation to estimate posterior mean, variance, skewness, and Pr(theta < P0)

### CHOICES HERE

genlink.data <- c(125,18,20,34);    label="original"  ### original
genlink.data <- c(14,0,1,5);        label="small"   ### smaller
genlink.data <- c(125,18,20,34) *3; label="big"  ### bigger

y1 = genlink.data[1]; y2 = genlink.data[2]; y3 = genlink.data[3]; y4 = genlink.data[4]

log.post = function(theta)
  y1*log((2+theta)/4)+(y2+y3)*log((1-theta)/4)+y4*log(theta/4)
#### Inefficient, but cute.
deriv.log.post = deriv(body(log.post), "theta", log.post, hessian=TRUE)
deriv.log.post
D.log.post = function(theta) c(attr(deriv.log.post(theta), "gradient"))
DD.log.post = function(theta) c(attr(deriv.log.post(theta), "hessian"))

moments = sapply(1:3,	function(Power) {
  log.g = function(theta) Power*log(theta)
  log.post.g <- function(theta)  log.post(theta) + log.g(theta)
  deriv.log.g = deriv(body(log.g), "theta", log.g, hessian=TRUE)
  D.log.g = function(theta) c(attr(deriv.log.g(theta), "gradient")) 
  DD.log.g = function(theta) c(attr(deriv.log.g(theta), "hessian")) 
  
  D.log.post.g = function(theta) D.log.post(theta) +  D.log.g(theta) 
  DD.log.post.g = function(theta) DD.log.post(theta) + DD.log.g(theta)
  
  theta.hat = print(optimize(log.post, 0:1, maximum=TRUE)$maximum)
  theta.hat.g = print(optimize(log.post.g, 0:1, maximum=TRUE)$maximum)
  
  sig.hat = print(DD.log.post(theta.hat))
  sig.hat.g = print(DD.log.post.g(theta.hat))
  sig.hat.g/sig.hat * exp(log.post.g(theta.hat.g) - log.post(theta.hat))
}
)
names(moments) = c("first", "second", "third")
moments
post.mean.lap = moments[1]
post.variance.lap = moments[2] - moments[1]^2
post.sd.lap = sqrt(post.variance.lap)
skewness = print((moments[3] + 2*moments[1]^3 - 3*moments[2]*moments[1])
                 /post.variance.lap^(3/2))

#### Plot results

#quartz(width=10, height=6)  ## creates a new graph on any platform 
par(mfrow=c(1,2))
npoints = 5000
theta.seq <- seq(0,1,length=npoints)
logliks  <- log.post(theta.seq)
normalizer <- sum(exp(logliks)) / npoints

plot(theta.seq, (exp(logliks)/normalizer), 
     type="l", col="darkgreen", ylab="likelihood")
lines(theta.seq, dnorm(theta.seq, mean=post.mean.lap, sd=post.sd.lap), col="red")
mtext("red: Laplace", col="red", side = 3)
mtext("green: exact", col="darkgreen", side = 3, line = 1)

plot(theta.seq, log(exp(logliks)/normalizer), 
     type="l", col="darkgreen", ylab="log likelihood")
lines(theta.seq, log(dnorm(theta.seq, mean=post.mean.lap, sd=post.sd.lap)), col="red")
mtext("red: Laplace", col="red", side = 3)
mtext("green: exact", col="darkgreen", side = 3, line = 1)


#  Not necessary here: log.post.vectorized = Vectorize(log.post)
rbind(laplace=moments,
      exact=sapply(1:3, function(Power)
  integrate(function(theta) 
    exp(log.post((theta))) * theta^Power
    , lower=0, upper=1)$value/
    integrate(function(theta) 
      exp(log.post(theta))
      , lower=0, upper=1)$value 
))

