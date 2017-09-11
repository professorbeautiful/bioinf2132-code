###########################################################################################################
#####   	EXERCISE:   STUDY THIS CODE.  ANYTHING YOU DON'T UNDERSTAND, BE PREPARED TO ASK ABOUT.
#####		EXERCISE:	CALCULATE THE POSTERIOR DISTRIBUTION ALGEBRAICALLY FOR THE DISCRETE PRIOR WE USED HERE.  CONFIRM AGREEMENT WITH THE LIMITING DISTRIBUTION DERIVED HERE.
#####		EXERCISE:	CALCULATE THE POSTERIOR DISTRIBUTION ALGEBRAICALLY FOR THE FLAT CONTINUOUS PRIOR WE USED HERE.  CONFIRM AGREEMENT WITH THE LIMITING DISTRIBUTION DERIVED HERE.  HOW FINE DOES THE GRID HAVE TO BE?
###########################################################################################################


Y = y.original   ### Rao's original data.

####  The "Mapping" loop--  doing the "integrals" in the case of a discrete prior.

MappingLoopRao = function(
	valuesForTheta=seq(0, 1, by=0.2),   ### where we put prior point probabilities.
	priorForTheta=rep(1, length(valuesForTheta))/length(valuesForTheta),  # the prior point probabilities.
	startingThetas=NULL,  #### If actual chains are desired, set this equal to a vector of starting thetas.
							####  or set   startingThetas="sample prior"
	nIterations=10,
	verbosity=1)
{
	Jtheta = length(valuesForTheta)
	distrXgivenThetaAndY = sapply(0:Y[1], 
		function(X1) dbinom(X1, Y[1], valuesForTheta))
	if(verbosity>1) catn("dim(distrXgivenThetaAndY) = ", dim(distrXgivenThetaAndY))
	###  This is a matrix. Each row is [X1 | theta_j, Y1]
	apply(distrXgivenThetaAndY, 1, sum)	### Check that the probabilities add to 1.
	
	### recall that
	# cdl.rao =
	# function(p, x)  dmultinom(x, prob=c(1/2, p/4, (1-p)/4, (1-p)/4, p/4))
	distrThetaGivenX1AndY = sapply(0:Y[1],
		function(X1) normalize(
			priorForTheta * sapply(valuesForTheta, cdl.rao,
				x=c(X1, Y[1]-X1, Y[2], Y[3], Y[4]))
			)
	)
	if(verbosity>1) catn("dim(distrThetaGivenX1AndY) = ", dim(distrThetaGivenX1AndY))
	apply(distrThetaGivenX1AndY, 2, sum)  ### Check that the probabilities add to 1.
	#### Sum over X1 to obtain K
	Karray = lapply(0:Y[1],
		function(X1) outer(distrXgivenThetaAndY[ , X1+1], distrThetaGivenX1AndY[ , X1+1]))
	if(verbosity>1) catn("length(Karray) = ", length(Karray))
	if(verbosity>1) print(str(Karray[1:3]))
	
	### Now we have to add up the component matrices in this list.
	Karray = pancake(Karray)  ### My function to convert a list into an array.
	if(verbosity>1) catn("dim(Karray) = ", dim(Karray))
	K = apply(Karray, 1:2, sum)  ## add up the component slices (matrices)
	if(verbosity>1) catn("dim(K) = ", dim(K))

	####  Beginning of iteration loop.
	
	currentGuessDistribution = priorForTheta
	allDistributions = matrix(NA, nrow=nIterations, ncol=Jtheta)   ### Initialize
	if(!is.null(startingThetas)) {
		if(length(startingThetas) == 1 & startingThetas == round(startingThetas)) 
			###   startingThetas is a single integer:  the NUMBER of starting thetas, sampled from the prior.
			startingThetaIndices = sample(1:Jtheta, startingThetas, replace=TRUE, prob=priorForTheta)
		else
			###	   startingThetas is a vector of actual Theta values to start at.
			startingThetaIndices = which(startingThetas == thetaValues)
		currentThetaIndices = startingThetaIndices   ### Initial values for each chain.
		####  Allocate space for all values of all chains.
		allThetaIndices = matrix(NA, nrow=nIterations, ncol=length(startingThetaIndices))
	}
	for(iterationNumber in 1:nIterations) {
		### Integrate the kernel K against the currentGuessDistribution
		currentGuessDistribution = currentGuessDistribution %*% K
		if(verbosity>1) catn(currentGuessDistributioncurrentGuessDistribution)
		allDistributions[iterationNumber, ] = currentGuessDistribution
		### And do the corresponding sampling if desired.
		if(!is.null(startingThetas)) {
			currentThetaIndices = sapply(currentThetaIndices, 
					function(thetaIndex)
						sample(1:Jtheta, 1, prob=K[thetaIndex, ])) 
			allThetaIndices[iterationNumber, ] = currentThetaIndices
		}
	}
	allDistributions = rbind(priorForTheta, allDistributions)
	result = allDistributions
	if(!is.null(startingThetas)) {
		allThetaIndices = rbind(startingThetaIndices, allThetaIndices)
		result = list(allDistributions=allDistributions, allThetaIndices=allThetaIndices)
	}
	attr(result, "call") = match.call()
	return(result)
}

# MappingLoopRao()
# loopResult = MappingLoopRao(startingThetas=rep(1:6, 1000))
loopResult = MappingLoopRao(startingThetas=10)
# MappingLoopRao(priorForTheta=c(0.95,0.01,0.01,0.01,0.01,0.01), nIterations=100)

loopResult = MappingLoopRao(priorForTheta=c(0.95,0.01,0.01,0.01,0.01,0.01), nIterations=4,
	startingThetas="sample prior"
)
#### This shows that the distributions are correct.
for(iterNum in 1:4) 
	print(rbind(
		loopResult[[1]][iterNum, ], 
		sapply(1:ncol(loopResult[[1]]), 
			function(i) mean(loopResult[[2]][iterNum,] == i))
	))

#####  Here's a cute way of multiplying many copies of K without a loop:
pow=2
KpowerString = print(paste(rep("K", pow), collapse=" %*% "))
KpowerExpression = print(parse(text=KpowerString))
eval(KpowerExpression)



 