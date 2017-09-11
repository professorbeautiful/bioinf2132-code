## confirm that the jacobean for lognormal does make the density integrate to 1.

integrate(dnorm, lower=-Inf,upper=Inf)

integrate(function(u) dnorm(log(u)), lower=0,upper=Inf)
#Bad!

integrate(function(u) dnorm(log(u))  / u, lower=0,upper=Inf)
#  Good!