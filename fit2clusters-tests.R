#require(IdMappingAnalysis)
#args(fit2clusters)
# fit2clusters = IdMappingAnalysis::fit2clusters

### First, with essentially no measurement noise.
fit = fit2clusters(testMe=TRUE,   
                   estimatesOnly=F,
                   simAlpha=1e3, simMeanSD=1e-10)
### Add measurement noise.
fit = fit2clusters(testMe=TRUE, Ntest = 1000,
                   simAlpha=5, simMeanSD=0.1, #simV=c(1,1), 
                  estimatesOnly=F)
showFit2(fit)

#### Add constraint, when it is wrong
fit = fit2clusters(testMe=TRUE,  simV=c(0.01,0.01), 
                   psi0Constraint = 0, simPsi=c(-0.5, 0.5),
                   estimatesOnly=F)
showFit2(fit)

fit = fit2clusters(testMe=TRUE,  simV=c(0.1,0.1), 
                   estimatesOnly=F)
showFit2(fit)

