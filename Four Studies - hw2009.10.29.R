#Let  = the probability of responding.  
#Assume a flat prior for the log odds, that is a Beta(0,0) prior for .
#Also repeat c,d,e below with Beta(1,1) and Beta(2,2).
# 
#For each study, calculate:
#
#a)	The observed probability and odds.
#b)	The P-value for testing H0: =1/2 versus HA: != 1/2.
#c)	The posterior distribution of the probability of responding.
#d)	The posterior distribution of the odds of responding.
#e)	Bayes factor, considering the prior Beta(0,0) for the alternative hypothesis.

scan(sep=",", what=list(character(0), numeric(0), numeric(0)))Study_1,15,5Study_2,115,86Study_3,1046,954Study_4,1001445,998555

###  A blank line ends the input.
hw2009.10.29 = as.data.frame(.())
names(hw2009.10.29) = c("study", "resp", "nonresp")

hw2009.10.29 = within(hw2009.10.29, {
		N <- resp+nonresp;
		proportion <- resp/(resp+nonresp)  ;
		odds <- resp/nonresp; 
		Pvalue <- 2*(1-pbinom(resp-1, resp+nonresp, 0.5))
		post.a.00 <- resp
		post.b.00 <- nonresp
		post.a.11 <- resp + 1
		post.b.11 <- nonresp + 1
		post.a.22 <- resp + 2
		post.b.22 <- nonresp + 2
	}
)  
with(hw2009.10.29, {
	postprobrespond.00 =  post.a.00/(post.a.00+post.b.00) # expectation of beta(a, b)
	postprobrespond.11 =  post.a.11/(post.a.11+post.b.11) # expectation of beta(a, b)
	postprobrespond.22 =  post.a.22/(post.a.22+post.b.22) # expectation of beta(a, b)
	print(cbind(postprobrespond.00, postprobrespond.11, postprobrespond.22))
})
##

## BAYES FACTOR, relative to prior 1,1
# #e)	Bayes factor, considering the prior Beta(0,0) for the alternative hypothesis.

# pr(data | p=1/2) 
##  "marginal"  probability under H0.
hw2009.10.29$m0 = with(hw2009.10.29, dbinom(resp, resp+nonresp, 1/2))

###  We will need this to get the marginal under HA.
# pr(data | uniform) = beta-binomial
Beta = function(a,b) exp(lgamma(a)+lgamma(b)-lgamma(a+b))
dbetabinomial = function(a,b,x,n)
	Beta(a+x,b+n-x)/Beta(a,b)*choose(n,x)
##Claim:  Beta(a,b) = 1/choose(a+b-2,b-1) / (a+b-1).   Validated.
###  a=15;b=17;Beta(a,b) * choose(a+b-2,b-1) * (a+b-1)

hw2009.10.29$mA00 = with(hw2009.10.29, dbetabinomial(0, 0, resp, resp+nonresp))
hw2009.10.29$mA11 = with(hw2009.10.29, dbetabinomial(1, 1, resp, resp+nonresp))
hw2009.10.29$mA22 = with(hw2009.10.29, dbetabinomial(2, 2, resp, resp+nonresp))

print(cbind(m0, mA))
m0/mA

# where's the problem?
with(hw2009.10.29, Beta(resp+nonresp, resp))
with(hw2009.10.29, choose(resp+nonresp, resp))

### How about normal approximation?
hw2009.10.29$m0approx1 = with(hw2009.10.29, dnorm(resp, N/2, sqrt(N/4)))  
### A good approximation for m0.

#### Marginal variance of X.
varAlist = with(hw2009.10.29, 
	sapply(list(c(1e-12,1e-12),c(1,1),c(2,2), c(1e10,1e10)),
		function(ab) {
				a=ab[1]; b=ab[2]; 
				Mean=a/(a+b);
				N*Mean + N*(N-1)*Mean*(a+1)/(a+1+b) - N^2*Mean^2
		}
))
hw2009.10.29$varX_A_00 = varAlist[ , 1]
hw2009.10.29$varX_A_11 = varAlist[ , 2]
hw2009.10.29$varX_A_22 = varAlist[ , 3]
hw2009.10.29$varX_A_BigBig = varAlist[ , 4]

varAlist  ### Confirmed via      var(rbinom(1e6, hw2009.10.29$N[3], rbeta(1e6, 2,2)))
m0approx1 = with(hw2009.10.29, dnorm(resp, N/2, sqrt(N/4)))  ### A good approximation for m0.
hw2009.10.29$m0approx2 = with(hw2009.10.29, dnorm(resp, N/2, sqrt(varX_A_BigBig)))  
### Another good approximation for m0.  Validated.
hw2009.10.29[cq(m0,m0approx1,m0approx2)]



hw2009.10.29$mAapprox00 = with(hw2009.10.29, dnorm(resp, N/2, sqrt(varX_A_00)))  ### A good approximation for mA?
hw2009.10.29$mAapprox11 = with(hw2009.10.29, dnorm(resp, N/2, sqrt(varX_A_11)))  ### A good approximation for mA.
hw2009.10.29$mAapprox22 = with(hw2009.10.29, dnorm(resp, N/2, sqrt(varX_A_22)))  ### A good approximation for mA.
### Validation: with(hw2009.10.29, qqplot(rbinom(1e6, N[3], rbeta(1e6, 2,2)), rnorm(1e6,N[3]/2, sqrt(varX_A_22[3]))));    diagonal()

hw2009.10.29[cq(mA00,mAapprox00,mA11,mAapprox11,mA22,mAapprox22)]

####  Bayes factors
hw2009.10.29$BF00 = hw2009.10.29$m0 / hw2009.10.29$mAapprox00
hw2009.10.29$BF11 = hw2009.10.29$m0 / hw2009.10.29$mAapprox11
hw2009.10.29$BF22 = hw2009.10.29$m0 / hw2009.10.29$mAapprox22

