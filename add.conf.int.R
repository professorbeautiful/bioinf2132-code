#CONFIDENCE INTERVALS

add.conf.int <- function(f, p, mle, plot=TRUE, ...) {
	maxloglik <- f(mle, ...)
	is.rejected <-  ! (f(p, ...) 
				> maxloglik - qchisq(0.95, df=1)/2)
	confint.lower <- min(p[!is.rejected])
	confint.upper <- max(p[!is.rejected])
	if(plot) {
		abline(h=max(f(p, ...)) - qchisq(0.95, 1), col="red")
		ylim<-par()$usr[3:4]
		lines(rep(confint.lower,2), ylim, col="blue")
		lines(rep(confint.upper,2), ylim, col="blue")
	}
	return(c(confint.lower, confint.upper))
}