alpha = 4
beta = 1

theThetas = rgamma(1000, shape = alpha, scale = beta)
plot(density(theThetas))
rug(theThetas)

theCounts = rpois(10000, theThetas)
plot(theThetas, theCounts)

hist(theCounts)
temp = table(theCounts)
names(temp)
temp[as.character(setdiff(0:max(theCounts), 
             as.numeric(names(temp))))] = 0
temp = temp[order(as.numeric(names(temp)))]
rbind(temp,
  round(dnbinom(x = 0:17, size = alpha, prob = 1/2) * 1000)
)
