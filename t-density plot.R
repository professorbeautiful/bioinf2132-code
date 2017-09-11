###   t-density plot

plot(x<-seq(-5,5,length=200), dt(x, 1), type="l", col="red", ylim=c(0,.4), ylab="density")
title("Density of the t distribution")
for(df in c(2:10))
	lines(x, dt(x,df))
lines(x, dt(x,30), col="blue")
lines(x, dnorm(x), col="green")
dflist<-c(1:10,30,Inf)
legend(3, 0.3, 
	legend = paste("df =", dflist),
	lty=rep(1, length(dflist)),
	text.col=c("red", rep("black", length(dflist)-3), "blue", "green")
	)
pngSave("t-density")

plot(x<-seq(1.9, 3 ,length=200), pt(x, 1), type="l", col="red", ylab="cum prob", ylim=c(.7,1))
title("Cum prob of the t distribution")
for(df in c(2:10))
	lines(x, pt(x,df))
lines(x, pt(x,30), col="blue")
lines(x, pnorm(x), col="green")
dflist<-c(1:10,30,Inf)
legend(2.8, 0.85, 
	legend = paste("df =", dflist),
	lty=rep(1, length(dflist)),
	text.col=c("red", rep("black", length(dflist)-3), "blue", "green")
	)
pngSave("t-cum")

