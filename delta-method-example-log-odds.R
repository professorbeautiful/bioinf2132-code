logOR = log((10/14)/(12/64))

varLogOR = sum(1/c(10, 14, 12, 64))
alpha = 0.05
exp(
  logOR + c(-1,0,1) * qnorm(1 - alpha/2) * sqrt(varLogOR)
)
