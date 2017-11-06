ecmo2Data = matrix(c(9, 0, 6, 4), nrow=2)
antilogit = function(z) exp(z)/(1+exp(z))
lik2 = function(ab) {
  a = ab[1]; b = ab[2]
  p = antilogit(a + b*(1:0))
  prod(p^ecmo2Data[1,] * (1-p)^ecmo2Data[2,])
} 
require("mvtnorm")
prior2 = function(ab, a0=0, b0=0, va=5, vab=0, vb=5)
  dmvnorm(ab, c(a0, b0), 
          matrix(c(va, vab, vab, vb), nrow=2))
par(mfrow=c(1,2))
avec = seq(-1, 3, length.out=20)
bvec = seq(-1, 10, length.out=20)
abgrid = expand.grid(avec, bvec)
prior2Matrix = matrix(prior2(abgrid), nrow=length(avec))
contour(avec, bvec, prior2Matrix, col="red")
lik2Matrix = matrix(apply(abgrid, 1, lik2), nrow=length(avec))
contour(avec, bvec, lik2Matrix, add=T, col="blue")

posterior2Matrix = prior2Matrix * lik2Matrix
contour(avec, bvec, posterior2Matrix, col="green")

