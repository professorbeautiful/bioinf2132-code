library(rstan)
# "This model is not directly supported by Stan because it involves discrete
# parameters"
###  ERROR.  Check the manual how to do mixture models.

stan_fit2_modelcode <- "
data {
  int<lower=1> Ntest;
  real Y[Ntest];
  real<lower=0> Ysigsq[Ntest];
} 

parameters {
 real simPi1;     // Pr(G = 1)
// real simAlpha; //known
// real simBeta; //known
 vector[2] simV;  // dispersion from the group mean.
 vector[2] simPsi; // the means for the 2 groups
 vector[Ntest] G; // group, 0 or 1. 
// stan cannot handle discrete parameters!
} 

model {
 for(j in 1:2) {
   simPsi[j] ~ normal(0, 10);
 }
 for(i in 1:Ntest) {
   G[i] ~ bernoulli(simPi1);
   mu[i] = simPsi[G[i]+1];
   Vtotal[i] = (simV[G[i]+1] + Ysigsq[i]);
   Y[i] ~ normal(mu[i], Vtotal[i]);
 }
} 
"
# RStan accepts data as a list, an environment, or a vector of characters
# containing object names in the working space.

Ntest = 50
simPi = c(2/3, 1/3)
simV = c(0,0)
simPsi = c(0, 0.4)
simV = c(0.05^2, 0.05^2)
simAlpha = 5
simBeta = 400

G = rbinom(Ntest, 1, simPi[2])
Ysigsq = rgamma(Ntest, simAlpha, simBeta)
Vtotal = simV[G+1] + Ysigsq
Y = simPsi[G+1] + rnorm(Ntest)*sqrt(Vtotal)

fit <- stan(data = c("Y","Ntest","Ysigsq"),
            model_code = stan_fit2_modelcode, model_name = "example-fit2c", 
             iter = 2012, chains = 3, sample_file = 'norm-fit2c.csv',
            verbose = F) 
###  ERROR.  Check the manual how to do mixture models.