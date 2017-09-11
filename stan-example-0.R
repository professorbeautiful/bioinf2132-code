stanmodelcode <- "
data {
  int<lower=0> N;
real y[N];
} 

parameters {
real mu;
} 

model {
mu ~ normal(0, 10);
y ~ normal(mu, 1); 
} 
"

y <- rnorm(20) 
dat <- list(N = 20, y = y); 
require(rstan)
fit <- stan(model_code = stanmodelcode, model_name = "example", 
            data = dat, iter = 2012, chains = 3, sample_file = 'norm.csv',
            verbose = TRUE) 
print(fit)

# extract samples 
find("extract")  ### Careful.  
e <- rstan::extract(fit, permuted = TRUE) # return a list of arrays 
mu <- e$mu 

m <- rstan::extract(fit, permuted = FALSE, inc_warmup = FALSE) # return an array 
print(dimnames(m))

# using as.array directly on stanfit objects 
m2 <- as.array(fit)
