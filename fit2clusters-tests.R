require(IdMappingAnalysis)
args(fit2clusters)
?fit2clusters
### First, with essentially no measurement noise.
fit = fit2clusters(testMe=TRUE, seed=NA,  
                   estimatesOnly=F,
                   simAlpha=1e3, simBeta=1e10)
showFit2 = function(fit) {
  head(fit)
  names(attributes(fit))
  attr(fit, 'estimates')
  plot(fit$Y, fit$posteriorOdds, 
       log="y")
  abline(h=c(1/100, 1, 100), col='red')
  text(-0.5, 100, "odds=100/1", adj=c(0,0), col='red')
  text(0.5, 1/100, "odds=1/100" , adj=c(0,1), col='red')
  pi1 = attr(fit, "estimates")["estimated", "pi1"]
  piVector = c(0, pi1)
  posteriorMeans = rowSums(outer(fit$postProb, piVector ))
  plot(fit$Y, posteriorMeans)
}
showFit2(fit)

###Now, add measurement noise.
fit = fit2clusters(testMe=TRUE, seed=NA, simV=c(0.01,0.01), 
                   estimatesOnly=F)
showFit2(fit)

#### Add constraint
fit = fit2clusters(testMe=TRUE, seed=NA, simV=c(0.01,0.01), psi0Constraint = 0,
                   estimatesOnly=F)
showFit2(fit)



fit = fit2clusters(testMe=TRUE, seed=NA, simV=c(0.1,0.1), 
                   estimatesOnly=F)
showFit2(fit)

