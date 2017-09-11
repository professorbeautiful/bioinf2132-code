###  t distribution:   Checking the lemma

# generate data via PART A

nu1 = 3
S1 = 30
mu1 = 20
n1 = 10
Nrep = 1000000

phiA = (rchisq(Nrep, nu1)  / S1)   ^(-1)
muA = rnorm(Nrep, mu1, sqrt(phiA/n1))

muB = rt(Nrep, nu1) * sqrt(S1/nu1/n1) + mu1
phiB = (rchisq(Nrep, nu1+1) / (S1+n1*(muB-mu1)^2))   ^(-1)

	quartz(width=10, height=5)
	par(mfrow=c(1,2))
	qqplot(phiA, phiB, log="xy"); diagonal()
	qqplot(muA, muB); diagonal()

var(muB)
S1/nu1/n1 * nu1/(nu1-2)  #### Perfect.

### Compare with E(var) + var(E)
mean(phiA/n1) + 0