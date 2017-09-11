e.m.algorithm <- function(data, e.suff.stat, cdl.maximizer, phi.start, ...) {
	#assign("cdl.maximizer.temp", cdl.maximizer, frame=0, immediate=T)
	#assign("e.suff.stat.temp", e.suff.stat, frame=0, immediate=T)
	#on.exit(try(remove(c("cdl.maximizer.temp","e.suff.stat.temp"),frame=0)))
	### The previous lines are necessary in splus, 
	###  with use of e.suff.stat.temp and cdl.maximizer.temp below. 
	### The reason: Scoping issue in Splus. 
	Qmaximizer <- function(data, phi) {
		e.suff.stat.result <- (e.suff.stat(data,phi))
		cdl.maximizer(e.suff.stat.result)
	}
	em.algorithm (data, Qmaximizer, phi.start, ...)
}

### To test it: 

### Rao genetics data
e.suff.stat.rao <- function(data, phistar) {
	y<-data
	fracstar <- phistar/(2+phistar)
	x <- c(c(1-fracstar, fracstar)*y[1], y[-1])
	return (x) 
}

cdl.maximizer.rao <- function(complete.data) {
	x <- complete.data
	phi.opt <- (x[2]+x[5])/(x[2]+x[3]+x[4]+x[5])
	return (phi.opt) 
}

rao.data = c(125,18,20,34)
cdl.maximizer.rao(e.suff.stat.rao(rao.data , .62682))

e.m.result.rao <- e.m.algorithm (rao.data, e.suff.stat.rao, cdl.maximizer.rao, phi.start = 0.1, logodl.rao,
		Q=Q.rao, H=H.rao
)

options(digits=15)
e.m.result.rao 
options(digits=8)
e.m.result.rao 

e.m.algorithm(NULL, 
		e.suff.stat=function(data,phi){cat("I am e.suff.stat\n"); return(NULL)}, 
		cdl.maximizer=function(suffstat){cat("I am cdl.maximizer.\n"); return(suffstat)}, 
		NULL)

