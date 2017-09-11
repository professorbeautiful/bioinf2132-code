ht = c(h=9, t=3)
ab = c(a=1/2, b=1/2)
yMax = 4
thGap = 0.01
th <- seq(0,1,by=thGap)
plot(th,
     xlab="theta", ylab="dbeta",
     dbeta(th, shape1=ab[1], shape2=ab[2]),
     ylim=c(0, yMax),
     type="l",
     lwd=3, col="blue"
)

likelihood = dbinom(
  x = ht["h"], size = sum(ht), p=th)
normalizedLikelihood = likelihood/sum(likelihood)
lines(th, normalizedLikelihood/thGap
      , col="green")

abstar = c(a=ab["a"] + ht["h"],
           b=ab["b"] + ht["t"])
lines(th, 
      dbeta(th, 
            shape1=abstar[1], shape2=abstar[2]),
      col="red"
)

legend(x = 0.2, y = 3, text.col = c("blue","green", "red"),
       legend = c("prior", "lik (normed)", "posterior"))
