#####    Simple decision theory problem

simple.bayes.vectorized <- function (
	prev =  0.001, 
	sens = 0.99,	spec = 0.99,
	L.SW = 1,    #### inaction  (wait)
	L.HT = 1,    #### treatment
	data = "P"  #### or "N"
) {
	prior = c(prev, 1-prev)
	prior.odds = prev/(1-prev) 
	if (data=="P") {
		model.S = sens
		model.H = 1-spec
	}
	else if (data=="N"){
		model.S = 1-sens
		model.H = spec
	}
	else stop("data should be either P or N")   #### Use "stop" when a "fatal" error is detected.
	lik.ratio = model.S/model.H
	posterior.odds = 	prior.odds * lik.ratio
	posterior.S = posterior.odds/(1+posterior.odds)
	posterior.H = 1/(1+posterior.odds)
	E.loss.W = posterior.S * L.SW
	E.loss.T = posterior.H * L.HT
	return (data.frame(	E.loss.W, E.loss.T))
}

plot.simple.bayes= function (
	prev =  0.001, 
	sens = 0.99,
	spec = 0.99,
	L.SW = 1,    #### inaction  
	L.HT = 1,    #### treatment
	data = "P",  #### or "N"
	...
) {
	params = list(prev, sens, spec, L.SW, L.HT)
	paramnames = c("prevalence","sensitivity","specificity","Loss_given_S_and_W","Loss_given_H_and_T")
	which.varies = which(unlist(lapply(params, length)) > 1)
	print(paste(which.varies, paramnames[which.varies]) )
  
	if (length(which.varies) != 1)
		stop("Exactly one parameter must vary")
	default.value = formals(plot.simple.bayes)[[which.varies]]
	horiz.values = params[[which.varies]]
	horiz.label = paramnames[which.varies]
	cat(horiz.label, default.value, "\n")
	E.loss = simple.bayes.vectorized(prev, sens, spec, L.SW, L.HT, data)
	E.loss.W = E.loss[ , 1]
	E.loss.T = E.loss[ , 2]
	plot (c(horiz.values, horiz.values), c(E.loss.W, E.loss.T), pch=" ", xlab=horiz.label, 
		ylab="Bayes expected loss", ...)
	lines (horiz.values, E.loss.W, pch= "W", col="red", lwd=3)
	lines (horiz.values, E.loss.T, pch= "T", col="green", lwd=3)
	points(horiz.values, E.loss.W, pch= "W", col="red", lwd=3)
	points(horiz.values, E.loss.T, pch= "T", col="green", lwd=3)
	abline(v=default.value, lty=3, lwd=5, col="blue")
	title (paste("Varying", horiz.label))
	mtext(paste("data = ", data))
	cat("Click on the graph where you want the legend's upper left point to go.\n")
	legend(locator(1),       #### NOTE:  "locator()" is interactive.
		c("wait",  "treat", "default value"), col=c("red","green","blue"),
		pch="WT.",
		lty=c(1,1,3), lwd=c(3,3,3) 
		)
	dev.copy(png, file=paste("simplebayes", horiz.label, "png", sep="."))  ### png = Portable Network Graphics
	dev.off()     ####  Closes the output graphic file.
	return(E.loss)
}

plot.simple.bayes(L.SW=seq(0.5,20,1))
plot.simple.bayes(L.HT=seq(0,1,0.05))
plot.simple.bayes(prev=1/10^seq(0.5,4,0.1), log="x")  ##### NOTE: use of log and "..."
plot.simple.bayes(sens=seq(0.1,1,.01))
plot.simple.bayes(spec=seq(0.99,1,.0002))
invisible( plot.simple.bayes(L.HT=seq(0,1,0.05))  )  #### NOTE: Suppresses automatic printing the return value.

plot.simple.bayes(data="N",    #### NOTE:  overriding the default arg value.
	prev=1-1/10^seq(0.5,4,0.1))

