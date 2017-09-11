
##The model should be in terms of Pr(X1,X2|Y).

twoTestData.original = array(
	c(	   
	700,221,139,467,		 
	238,33,141,106, 
	938,254,280,573		 
	)
	,
	dim=c(2,2,2),
	dimnames=list(
	c("EST+", "EST–" ),	
	c("CAD","noCAD"),
	c("EKG-","EKG+" )
	)	
)

twoTestData = aperm(twoTestData.original, c(1,3,2))
twoTestData = twoTestData[ 2:1, , ]

apply(twoTestData, 3, sum)
CADmodel = 
  sweep(twoTestData, 3, apply(twoTestData, 3, sum), "/")

priorArray = CADmodel
priorArray[ , , ] = rep(prior, each=4)

priorArray * CADmodel
rep(prior, each=4) * CADmodel

posterior.from.2.tests = function(
	prevalence=0.58,		# this the the prevalence, friend
	model=CADmodel,		# this is condit'l distr of alltest results, as an array.
	marginalizeOver,		#
	verbose=TRUE			#
) {
	prior = c(prevalence, 1-prevalence)
	priorArray = model
	#priorArray[ , , ] = 1:8  ## to see where entries go.
	priorArray[ , , ] = rep(prior, each=4)
	jointAll = priorArray * model
	if(verbose){
		cat("\n=======jointAll\n")
		print(jointAll)
	}
	dataDims = 1:2
	if(!missing(marginalizeOver)) {
		conditionOn = 1   ### whether is was dim 1 or 2, it's 1 now!
		joint = apply(jointAll, (1:3)[-marginalizeOver], sum)
	} else {
		conditionOn = dataDims
		joint = jointAll
	}
	if(verbose){
		cat("\n=======joint\n")
		print(joint)
	}
	thetaDim = length(dim(joint))  ##  = 3 if no marginalizing
	marginalOfData = apply(joint, conditionOn, sum)
	if(verbose){
		cat("\n=======marginalOfData\n")
		print(marginalOfData)
	}
	posterior = sweep(joint, conditionOn, 
		marginalOfData, "/")
	if(verbose)  cat("\n\n=======posterior\n")
	return(posterior)
}


posterior.from.2.tests()
posterior.from.2.tests(marginalizeOver=1)
posterior.from.2.tests(marginalizeOver=2, verbose=FALSE)
posterior.from.2.tests(prevalence=0.05, marginalizeOver=2, verbose=FALSE)

BayesRule.from.2.tests = function(
	posterior,
	lossRatio=1
) {
	posteriorDF = as.data.frame(aperm(posterior))
	posteriorOdds = posteriorDF[1,] / posteriorDF[2,]
	rule = ifelse(posteriorOdds * lossRatio > 1, "do", "don't")
	dimnames(rule)[[1]] = "action"
	return(rule)
}

BayesRule.from.2.tests(posterior.from.2.tests(verbose=F))
BayesRule.from.2.tests(posterior.from.2.tests(marginalizeOver=2, verbose=F))
BayesRule.from.2.tests(posterior.from.2.tests(marginalizeOver=2, verbose=F, prevalence=0.05))

loglin(twoTestData, list(1:3), fit=TRUE)
loglin(twoTestData, 1:3, fit=TRUE)
loglin(twoTestData, list(1:2, 2:3, c(1,3)), fit=TRUE)

twoTestData.IndependentFit =
	loglin(twoTestData, list(1:2, 2:3, c(1,3)), fit=TRUE)$fit
CADmodelIndependent = 
  sweep(twoTestData.IndependentFit, 3, apply(twoTestData.IndependentFit, 3, sum), "/")

posterior.from.2.tests(verbose=FALSE, model=CADmodelIndependent)
posterior.from.2.tests(verbose=FALSE, model=CADmodel)


expectedlossRatio = function(
	posterior,
	lossRatio=1
) {
	posteriorDF = as.data.frame(aperm(posterior))
	posteriorOdds = posteriorDF[1,] / posteriorDF[2,]
	posteriorOdds * lossRatio 
}

expectedlossRatio(posterior.from.2.tests(verbose=FALSE, 	marginalizeOver=2, model=CADmodel))

cutoff = uniroot(
	f=function(prev)
		1 - expectedlossRatio(
			posterior.from.2.tests(
				prevalence=prev, verbose=FALSE,
			 	marginalizeOver=2, model=CADmodel))[1,2],
	interval=c(0.05, 0.58)
	)$root
	
	