
MarkovProcess <- function(P, start=1, steps=10, reps=1, print.P=F) {
### "start" is an integer between 1 and dim(P), which is square.
### Each row of P should sum to 1, all element non-negative.
###  To do: check inputs, implement reps, implement for a start vector, 
###         make sure "rmultinomial" works if by chance runif returns 0 or 1,
###         generalize "rmultinomial".
	if(print.P) print(P)
	if(reps==1) {
		state <- start
		states <- rep(NA, steps+1)
		states[1] <- state
		for (step in 1:steps) {
			state <- states[step+1] <- rmultinom(1, 1, P[state, ])
		}
		return(states)
	}
	states <- matrix(NA, nrow=reps, ncol=steps+1)
	for (irep in 1:reps)
		states [irep, ] <- MarkovProcess(P, start, steps, 1)
	return(states)
}

MarkovProcess(P=matrix(c(1-(ptrans<-.999),ptrans,0,0,1-ptrans,ptrans,ptrans,0,1-ptrans),byrow=T, nrow=3), 1, 10, print.P=T)
plot(0:10, .Last.value)

MarkovProcess(P=random.walk.matrix(), start=5, steps=10, reps=5)


many.chains <- MarkovProcess(P=random.walk.matrix(), 
	start=5, steps=100, reps=100)
for (j in c(1:5,10,15,20,30)) 
	hist(many.chains[,j], breaks=0:10, main=paste(j))

many.chains <- MarkovProcess(P=random.walk.matrix(), 
	start=8, steps=100, reps=100)
for (j in c(1:5,10,15,20,30)) 
	hist(many.chains[,j], breaks=0:10, main=paste(j))

matpow <- function(P, n, expressionMethod=TRUE) {
	if(expressionMethod == FALSE) {
		temp <- P; for (i in 1:(n-1)) temp <- temp%*%P  #inefficient!!!
	}
	else {
		PpowerString = print(paste(rep("P", pow), collapse=" %*% "))
		PpowerExpression = print(parse(text=PpowerString))
		temp = eval(PowerpowerExpression)
	}
	return (temp)
} 