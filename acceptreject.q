# Goal: simulate binomial variates,
# using an accept-reject method.

# (1) find an "envelope" that can be sampled from.

p <- 0.20
n <- 25

m <- 5

plot(0:25,dbinom(0:25,25,p),type="b",cex=2)
lines(0:25,dpois(0:25,m),col=2)
f <- max( dbinom(0:25,25,p)/dpois(0:25,m))
envelop <- function(j) {dpois(j,m)*f}
lines(0:25,envelop(0:25), lty=2, lwd=2, col=2)


acceptprob <- dbinom(0:25,25,p)/envelop(0:25)
howmany <- rep(0,25+1)
howmany.ok <- rep(0,25+1)
for ( i in 1:1000){
	x <- rpois(1,m)
	if ( x <= 25) {
		u <- runif(1)
		wheretoput <- u*envelop(x)
		ok <- (u < acceptprob[x+1])
		howmany[x+1]<-howmany[x+1]+1
		if (ok)  howmany.ok[x+1]<-howmany.ok[x+1]+1
		if (ok)  color<-3
		else color <- 4
		points( x, wheretoput,col=color)
	} else
		ok <- F
}

