### guilty or innocent?  bayes factor calculation

BFnormal = function(n, tau, sigma, z)
	sqrt(n*tau^2 + sigma^2) * 
		exp(-1/2 * z^2/(n*tau^2 + sigma^2) / n*tau^2)
###  confusiong.  tau is the prior variance in the equation, but the error variance in the example.
### the example only gives zG=3 and zI = 4.


BFnormal = function(zG, zI, V, tau) {
		
		
zG = 3
zI = 4
dnorm(zG)/dnorm(zI)  ###  BF = 33.11   
