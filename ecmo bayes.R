ecmo2Data = matrix(c(9, 0, 6, 4), nrow=2)
ab = c(1,2)

lik2 = function(ab) {
	a = ab[1]
	b = ab[2]
	p = antilogit(a + b*(1:0))
	prod(
	  p^ecmo2Data[1,] * (1-p)^ecmo2Data[2,]
	  )
	#  or , dbinom(ecmo2Data[1,], ecmo2Data[1,] + ecmo2Data[2,], p)
	# (proportional),   choose(10,6)*0.0007986267 = 0.1677116
} 
# comparing 9/9 to 6/10.
fisher.test(ecmo2Data)

require("mvtnorm")

prior2 = function(ab, a0=0, b0=0, va=1, vab=0, vb=1)
		dmvnorm(ab, c(a0, b0), 
			matrix(c(va, vab, vab, vb), nrow=2))

avec = seq(-1, 3, length.out=20)
bvec = seq(-1, 10, length.out=20)
abgrid = expand.grid(avec, bvec)

par(mfrow=c(1,2))

prior2Matrix = matrix(prior2(abgrid), nrow=length(avec))
contour(avec, bvec, prior2Matrix)

lik2Matrix = matrix(apply(abgrid, 1, lik2), 
	nrow=length(avec))
contour(avec, bvec, lik2Matrix,
	add=T)

posterior2Matrix = prior2Matrix * lik2Matrix
contour(avec, bvec, posterior2Matrix)

#########   MCMC metropolis

source("/Users/Roger/Box Sync/teaching/bioinf2132,bios2063/bioinf2132-code/Plight-Pdark-posterior-new.R", local=TRUE)
twobytwodata = matrix(c(90, 100, 60, 40), nrow=2)
twobytwodata = ecmo2Data
jumpsigSD = .10
jumpsigMatrix = diag(jumpsigSD, 4)
# prior = function(p) {
#   prior2()
# }
metrop.out = metrop(chainlength = 2000,
                    start = rep(1/4, 4), 
                   jumpsig = jumpsigMatrix,
                   dens=function(p) 
                     return(
                       #prior(p) * 
                         dmultinom(x=c(twobytwodata), prob = p)),
                   rtrans=function(p,jumpsig) { ### randomly chosen proposal
                     proposal = antilogit(rmvnorm(n = 1,mean = logit(p), sigma = jumpsigMatrix))
                     proposal = proposal/sum(proposal)
                     TINY = 1e-15
                     proposal = pmax(proposal, TINY)
                     proposal = pmin(proposal, 1-TINY)
                     # cat('.')
                     return(proposal)
                   },
                   Jtrans=function(p1,p2,jumpsig)  ### prob density of proposal
                   {
                     jacobean <- diag(c(p2*(1-p2))); #1/(dlogit(p2)/dp2);
                     return(dmvnorm(x = logit(p2), 
                                    mean=logit(p1), sigma = jumpsigMatrix) 
                            * det(jacobean %*% solve(jumpsigMatrix)) )
                   }
)

pEvalues = metrop.out$values[ , 1] / (metrop.out$values[ , 1]+ metrop.out$values[ , 2])
summary(pEvalues)
pCvalues = metrop.out$values[ , 3] / (metrop.out$values[ , 3]+ metrop.out$values[ , 4])
summary(pCvalues)
### Convert to a and b
# Recall p = antilogit(a + b*(1:0))
aValues = logit(pCvalues)
bValues = logit(pEvalues) - aValues
summary(aValues)
summary(bValues)
points()

par(mfrow=c(2,2))
for(whichdim in 1:4)
  plot.metrop.output.1(metrop.out, whichdim = whichdim, ylab='value')
summary(metrop.out$values[ , ])
apply(metrop.out$values, 2, mean)
twobytwodata/sum(twobytwodata)
