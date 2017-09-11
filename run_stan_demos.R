library("rstan")

stan_demo(0)

grep(value=T, "cluster", stan_demo(0))

rm(list=ls())

model <- stan_demo(
  grep("naive-bayes.stan", stan_demo(0))
  )
