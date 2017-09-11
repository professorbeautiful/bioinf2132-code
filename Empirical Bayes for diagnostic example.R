# 

N = 100
TP = 30
eta0 = 0
eta1 = 0.2
eta = c(eta0, eta1)
p0 = 1/2
## Empirical Bayes for diagnostic example.
## 

a = 0.9
b = 0.8

prT = function() {
  Tvec = 0:N
  prOneEta = function(et) {
    dbinom(0:N, size = N, prob = 
    a*et + (1-b)*(1-et))
  }
  sapply(eta, prOneEta)
}

prTvalues = prT()
colSums(prTvalues)  ### Check.

plot(0:N, prTvalues[,1], col=NULL, yaxs='i')
lines(0:N, prTvalues[,1],col='green')
lines(0:N, prTvalues[,2],col='red')
abline(v=30)
# Posterior odds that eta = eta0

( (1-b)^TP * b^(N-TP) ) / 
  ((a*eta1+(1-b)*(1-eta1))^TP * ( (1-a)*eta1 + b*(1-eta1))^(N-TP))
dbinom(TP, N, 1-b) / dbinom(TP, N, a*eta1+(1-b)*(1-eta1))
### odds of eta0 = 0.086, Pr(eta=eta0 | T) =  0.079.
0.921*0.18/0.34
# 0.488