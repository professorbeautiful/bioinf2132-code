### DATA
Xpois = 221
# Xrange = ( Xpois -30):( Xpois +30)
XrangeGap = 1
XrangeUpper = 2.0
Xrange = seq( 0, Xpois +  sqrt(Xpois) * XrangeUpper, XrangeGap)

normalize = function(x) x/sum(x)

plot(x<-Xrange, normalize(dpois(x, Xpois)) , 
     type="l", ylim=c(0,0.028), lwd=2, ylab="Density")
abline(v= Xpois, lty=2)

### PRIOR
delta = 2
phi = 1/50
lines(Xrange, dgamma(Xrange, delta, rate=phi), col="blue")

### POSTERIOR
deltaStar  = delta + Xpois
phiStar = phi + 1
toleranceInterval = qgamma(c(0.05,0.95), deltaStar, phiStar)
lines(rep(toleranceInterval[1], each=2),
      c(0, dgamma(toleranceInterval[1], deltaStar, phiStar)), lty=3)
lines(rep(toleranceInterval[2], each=2),
      c(0, dgamma(toleranceInterval[2], deltaStar, phiStar)), lty=3)
lines(Xrange, dgamma(Xrange, deltaStar, rate=phiStar), col="red")

## normal approx
lines(x, dnorm(x, deltaStar/phiStar, sqrt(deltaStar)/phiStar), lty=2)
legend("topleft", legend = c("likelihood (normalized)",
                              "prior", "posterior",
                              "normal approximation"), col = c("black", "blue", "red", "black"),
       lty=c(1,1,2,1), lwd=c(6,3,3,3))