#####    A simple decision theory problem
simple.bayes <- function (
   ### Any argument except data can be a vector (one at a time!).  Input error checking is needed.
	prev =  0.001, 	### prevalence
	sens = 0.99,		### sensitivity
	spec = 0.99, 	### specificity
	L.SW = 1,    	### action = Wait  
	L.HT = 1,    	### action = Treat
	data = "P"  		### data="positive", or "N" for "negative"
) {
	prior = c(prev, 1-prev)
	prior.odds = prev/(1-prev) 
	if (data=="P") model.S = sens
	if (data=="N") model.S = 1-sens
	if (data=="P") model.H = 1-spec
	if (data=="N") model.H = spec
	lik.ratio = model.S/model.H
	posterior.odds = 	prior.odds * lik.ratio
	posterior.S = posterior.odds/(1+posterior.odds)
	posterior.H = 1/(1+posterior.odds)
	E.loss.W = posterior.S * L.SW
	E.loss.T = posterior.H * L.HT
	return (data.frame(	E.loss.W, E.loss.T))
	### Return value = Bayes expected losses for the two actions W and T.
}

plot.simple.bayes= function (
	prev =  0.001, 
	sens = 0.99,
	spec = 0.99,
	L.SW = 1,      
	L.HT = 1,    
	data = "P",  
	 ...				#### extra arguments for the graphics
) {
	params = list(prev, sens, spec, L.SW, L.HT)
	paramnames = c("prevalence", "sensitivity", "specificity", "Loss_given_S_and_W", "Loss_given_H_and_T")
	which.varies = which(unlist(lapply(params, length)) > 1)
	print(which.varies)
	if (length(which.varies) != 1)
		stop("Exactly one parameter must vary")
	default.value = formals(plot.simple.bayes)[[which.varies]]
	horiz.values = params[[which.varies]]
	horiz.label = paramnames[which.varies]
	cat(horiz.label, default.value, "\n")
	E.loss = simple.bayes(prev, sens, spec, L.SW, L.HT, data)
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
	legend(locator(1), 
		c("wait",  "treat", "default value"), col=c("red","green","blue"),
		pch="WT.",
		lty=c(1,1,3), lwd=c(3,3,3) 
		)
	dev.copy(png, file=paste("simplebayes", horiz.label, "png", sep="."))
	dev.off()
	return(E.loss)
}

plot.simple.bayes(L.SW=seq(0.5,20,1))
plot.simple.bayes(L.HT=seq(0,1,0.05))
plot.simple.bayes(prev=1/10^seq(0.5,4,0.1), log="x")
plot.simple.bayes(sens=seq(0.1,0.99,.01))
plot.simple.bayes(spec=seq(0.99,1,.0002))
plot.simple.bayes(L.HT=seq(0,1,0.05))



#####		Solving for the break-even value.  
#####   EXERCISE:   Install R. Run the following code.  
#####		Write a one-page commentary on the results---- especially the errors!

cat("If we let L.SW vary, its break-even value where rho(W)=rho(T) is L.SW = ",
	uniroot(function(L.SW) diff(unlist(simple.bayes(L.SW=L.SW))), interval=c(0, 100))$root, "\n")
cat("If we let L.HT vary, its break-even value where rho(W)=rho(T) is L.HT = ",
	uniroot(function(L.HT) diff(unlist(simple.bayes(L.HT=L.HT))), interval=c(0, 100))$root, "\n")
cat("If we let spec vary, its break-even value where rho(W)=rho(T) is spec = ",
	uniroot(function(spec) diff(unlist(simple.bayes(spec=spec))), interval=c(0.50, 1.0))$root, "\n")
cat("If we let spec vary, its break-even value where rho(W)=rho(T) is spec = ",
	uniroot(function(spec) diff(unlist(simple.bayes(spec=spec))), interval=c(0.50, 0.999999999999))$root, "\n")
cat("If we let sens vary, its break-even value where rho(W)=rho(T) is sens = ",
	uniroot(function(sens) diff(unlist(simple.bayes(sens=sens))), interval=c(0.50, 1))$root, "\n")
cat("If we let prev vary, its break-even value where rho(W)=rho(T) is prev = ",
	uniroot(function(prev) diff(unlist(simple.bayes(prev=prev))), interval=c(1e-6, 0.99))$root, "\n")

intervals = list(
	L.SW=c(0,100),
	L.HT=c(0,100),
	spec=c(0.50, 1-1e-7),
	sens=c(0.50, 1-1e-7),
	prev=c(1e-6,0.99))
sapply(cq(L.SW,L.HT,spec,sens,prev), function(par) {
	temp=try(uniroot(function(X)diff(unlist(
		eval(parse(text="simple.bayes(" %&% par %&% "=" %&% X %&% ")")))), 
			interval=intervals[[par]], tol=.Machine$double.eps^0.75)$root); 
		ifelse(class(temp)=="try-error", NA, temp)
		}
)
