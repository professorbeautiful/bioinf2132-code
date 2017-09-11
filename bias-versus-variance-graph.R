plot(P<-seq(-1,1,length=1000), (P+1)^2 / 2 + (-P+1)^2 / 2, col="black", type="l", 
     xlab="", ylab="", axes=F,
     ylim=c(0,2), cex=1.7, lwd=4)
lines(P, (P+1)^2 / 2 , col="blue", lwd=2)
lines(P, (-P+1)^2 / 2 , col="darkorange", lwd=2)


axis(1, at=c(-1,1), col.ticks="black", labels=c("",""))
mtext(side=1, expression(complexity), cex=1.7, line=1)

axis(2, at=c(0,2), col.ticks="white", labels=c("",""))
mtext(side=2, expression(bias^2), cex=1.7, line=1, col="darkorange")

axis(4, at=c(0,2), col.ticks="white", labels=c("",""))
mtext(side=4, "variance", cex=1.7, line=1, col="blue")

text(0, 1.6, "mean squared error\n=", cex=1.3)
text(-0.2, 1.34, expression(bias^2), cex=1.3, col="orange" )
text(0, 1.34, " + ", cex=1.3 )
text(0.3, 1.34, expression(variance), cex=1.3, col="blue" )


  # "simple\nhighly smooth\nlumped", 
  # "complex\nhighly rough\nsplit"
  # ))
