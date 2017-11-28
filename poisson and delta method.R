plot(lambda.vec<-exp(seq(log(0.01), 
	log(1000), length=50)), 
		sapply(lambda.vec, 
		function(lambda)
			var(sqrt(rpois(1e5, lambda)))),
  log="x")
abline(h=1/4)


###  Freeman-Tukey

plot(lambda.vec<-exp(seq(log(0.01), 
	log(1000), length=50)), 
		sapply(lambda.vec, 
		function(lambda) {
			Xvec <- rpois(1e5, lambda);
			var(sqrt(Xvec) + sqrt(Xvec+1))
		}),
  log="x",
	xlab=expression(lambda ~ "(Poisson mean)"),
	ylab=expression(var(sqrt(lambda)+sqrt(lambda+1)))
)
title("Freeman-Tukey stabilization")
abline(h=1)

plot(c(-3,3), c(-3,3), col="white")
lambda.vec <- 10^(seq(log10(0.01), log10(1000), length=6)) 
colorOffset = length(lambda.vec) - min(1, log10(lambda.vec[1]))
colors
invisible(sapply(lambda.vec, 
	function(lambda) {
			Xvec <- rpois(1e4, lambda);
			qq <- qqnorm(plot.it=FALSE,
				(sqrt(Xvec) + sqrt(Xvec+1) - sqrt(4*lambda + 1)))
			#print(str(qq))
			lines(qq$x[order(qq$x)], qq$y[order(qq$x)], 
			      col= which(lambda.vec==lambda)+colorOffset)#  (round(log10(lambda) + colorOffset)))
		}
		))
legend(-3, 3, legend=as.character(round(lambda.vec,3)), 
	lty=rep(1,length(lambda.vec)), text.col=colorOffset+1:length(lambda.vec) ,
	col=colorOffset+1:length(lambda.vec) ) #round(log10(lambda.vec) + colorOffset))
title("comparing Freeman-Tukey \nnormal approx \nto random Poisson")
