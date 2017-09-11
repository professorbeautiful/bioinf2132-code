
### DOCUMENTATION:   in the style of Roxygen
#'	@param prevalence  The prevalence of disease. Default=0.58.
#' @param model This is condit'l distr of all test results, as an array. The first two dimensions are tests with binary outcomes.  The third dimension is the disease state dimension; the first disease state corresponds to "disease present", the second to "disease absent". Default=CADmodel.
#' @return The posterior probability of disease given the data.

####### ANSWER ########

posterior.from.2.tests = function(
	X1="EKG+",
	X2="EST+",
	prevalence=prior[1],		
	model=CADmodel,
	verbose=TRUE
) {
	if(length(prevalence) > 1)
		return(sapply(prevalence, posterior.from.2.tests, X1=X1, X2=X2, model=model, verbose=verbose))
	if(length(X1) != 1) stop("length of X1 must equal 1")
	if(length(X2) != 1) stop("length of X2 must equal 1")
	if(is.na(X1) & is.na(X2))	return(prevalence)   #### no data!
	priorArray = model  ### Just to get the right dimensions and dimnames.
	priorArray[ , , ] = rep(c(prevalence, 1-prevalence), each=4)  ### to multiply prior times models.
	joint = priorArray * model  #### prior times models.  
	if(verbose) {cat("Joint Distribution:\n"); print(joint)}
	if( ! abs(1-sum(joint))<1e-8) stop("probabilities should add to 1")  ### = 1.
	if(verbose)  cat("data: X1,X2 = ", X1, ",", X2, "\nprevalence: ", prevalence, "\n")
	DiseaseIsPresent = 1  ### Required by documentation.
	if(is.na(X1))  	 		{		## X1 is not observed:  marginalize.
		numerator = sum(joint[ , X2, DiseaseIsPresent])
		denominator = sum(joint[ ,X2, ])
	} else if(is.na(X2))	 {		## X2 is not observed:  marginalize.
		numerator = sum(joint[X1 , , DiseaseIsPresent])
		denominator = sum(joint[X1 , , ])
	} else {
		numerator = joint[X1,X2, DiseaseIsPresent]
		denominator = sum(joint[X1,X2,])
	}
	return(numerator/denominator)
}
posterior.from.2.tests()   ## 0.627
posterior.from.2.tests(prevalence=0.58)   ## 0.625
posterior.from.2.tests(X1="EKG-")   ### 0.834
posterior.from.2.tests(prevalence=0.58, X1=NA)    ####   Compare with class notes:  OK.  0.768.
posterior.from.2.tests(X1=NA)  ###  0.7701    Because prior[1] is not exactly 0.58.
posterior.from.2.tests(X2=NA)  ###   0.523
posterior.from.2.tests(prevalence=0.58, X2=NA)  ###   0.520
posterior.from.2.tests(X2="EST-", prevalence=seq(0.1,0.9,0.1), verbose=FALSE)
plot(prevalence<-seq(0.1,0.9,0.1),   .Last.value,  xlab="prevalence", ylab="posterior probability of disease")
abline(a=0,b=1)

####  Exercise answers:
posterior.from.2.tests(prevalence=0.05, X1=NA)		###   0.112
posterior.from.2.tests(prevalence=0.05, X1=NA, verbose=T)   #### 0.112
posterior.from.2.tests(prevalence=0.58, X1=NA, X2="EST-")    #### 0.3046
posterior.from.2.tests(prevalence=0.05, X1=NA, X2="EST-")   #### 0.0164

### EXERCISE:   complete the body of this function
#'	@param posterior.  Prior or current (or imagined future) posterior probability of disease.   
#'	@param lossRatio.  Loss("treat", "no disease")	/ Loss("wait", "disease")	. Default=1.
#' @return  Optimal action.

####### ANSWER ########
expectedlossRatio  = function(
	posterior,
	lossRatio=1
) {
	posterior/(1-posterior) * lossRatio 
}
BayesRule.binaryAction.binaryState = function(
	posterior,
	lossRatio=1,
	...
) {
	ifelse(expectedlossRatio(posterior,lossRatio) > 1, "TREAT", "wait")
}

BayesRule.binaryAction.binaryState(.5)

rbind(
	seq(0,1,.1),
	BayesRule.binaryAction.binaryState(seq(0,1,.1)),
	BayesRule.binaryAction.binaryState(posterior.from.2.tests(X1="EKG-", X2="EST-", prevalence=seq(0,1,.1), verbose=F))
)

cutoffPrior = function(lossRatio=1, X1="EKG-", X2="EST-") {
	uniroot(function(Pr) 
		1 - expectedlossRatio(lossRatio, 
			posterior=posterior.from.2.tests(X1=X1, X2=X2, prevalence=Pr, verbose=F)),
		interval=0:1
	)
}	
cutoffPrior(X1=NA,X2="EST+")    #####   0.29  
cutoffPrior(X1=NA,X2="EST-")    #####   0.76  

#### EXERCISE C:   write a function with inputs of the form	
##    the disease prevalence, 
##   the model for [EKG | disease],  and the model for [EST | EKG, disease state],   
##   and returns the probability of disease conditional on EST result, marginal over EKG.

##  Use this:   [CAD | EST] = [EST | CAD] [CAD] /[EST]
#												 = SUM ( [EST | EKG, CAD] [EKG | CAD] )   [CAD] /[EST]

EKG_EST_Disease_joint = priorArray*CADmodel
EKG_Disease_margin = apply(EKG_EST_Disease_joint, c(1,3), sum)
EKG_Disease_margin_repeated = EKG_EST_Disease_joint  ### just for structure
EKG_Disease_margin_repeated = 
	abind(EKG_Disease_margin, EKG_Disease_margin, along=3)
EKG_Disease_margin_repeated = aperm(EKG_Disease_margin_repeated, c(1,3,2))
dimnames(EKG_Disease_margin_repeated)[[2]] = dimnames(EKG_EST_Disease_joint)[[2]]
EKGmodel = apply(CADmodel, c(1,3), sum)    ###OK:  apply(EKGmodel, 2, sum)
ESTmodelGivenEKG = EKG_EST_Disease_joint / EKG_Disease_margin_repeated
	###OK:  apply(ESTmodelGivenEKG, c(1,3), sum)   ==>  All ones.

ESTmodelGivenEKG[,,"CAD"] * (EKGmodel[,"CAD"])		## OK , same as  CADmodel[, , "CAD"]
ESTmodelGivenEKG[,,"noCAD"] * (EKGmodel[,"noCAD"])  ## OK, same as   CADmodel[, , "noCAD"]
###  So we can swap in a new EKGmodel.

posterior.from.2.tests.VERSION2 = function(
	# X1="EKG+",
	# X2="EST+",
	prevalence=prior[1],	
	X1model= EKGmodel,
	X2model= ESTmodelGivenEKG,
	verbose=TRUE
) {
##  Use this:   [CAD | EST] = [EST | CAD] [CAD] /[EST]
#												 = SUM ( [EST | EKG, CAD] [EKG | CAD] )   [CAD] /[EST]
#												 = SUM ( X2model * X1model )   [CAD] /[EST]
### TO BE CHECKED  ####
	thisPrior = c(prevalence, 1-prevalence)
	jointModel = abind(
					X2model[,,"CAD"] * X1model[,"CAD"],
					X2model[,,"noCAD"] * X1model[,"noCAD"]	, along=3)
	jointDistribution = jointModel * rep(thisPrior, each=4)  ### OK: Same as normalize(twoTestData)
	X2marginGivenDisease = apply(jointModel,  c(2,3), sum)  ## marginalize over EKG. Keep the EST columns.
	X2diseaseMargin = apply(jointDistribution,  c(2,3), sum)
	X2margin = apply(jointDistribution,  2, sum)  ## marginalize over EKG. Keep the EST columns.
	diseaseGivenX2 = X2diseaseMargin /  rep(X2margin, times=2)  ### Check:  apply(diseaseGivenX2, 1, sum)
	return( diseaseGivenX2	)
}
posterior.from.2.tests.VERSION2()
####  should be the same as before:      
posterior.from.2.tests(X1=NA, X2="EST+")   #### OK
posterior.from.2.tests(X1=NA, X2="EST-")   #### OK

###	Now we can use this function to carry the knowledge about [EST | EKG, disease] to other settings,
### 	by replacing both prevalence AND X1model.

###############  STOP HERE FOR NOW  ####################################

####  Investigation into conditional independence  #####

loglin(twoTestData, list(1:3), fit=TRUE)
loglin(twoTestData, 1:3, fit=TRUE)
loglin(twoTestData, list(1:2, 2:3, c(1,3)), fit=TRUE)

twoTestData.IndependentFit =
	loglin(twoTestData, list(1:2, 2:3, c(1,3)), fit=TRUE)$fit
CADmodelIndependent = 
   sweep(twoTestData.IndependentFit, 3, apply(twoTestData.IndependentFit, 3, sum), "/")

posterior.from.2.tests(verbose=FALSE, model=CADmodelIndependent)
posterior.from.2.tests(verbose=FALSE, model=CADmodel)


cutoff = uniroot(
	f=function(prev)
		1 - expectedlossRatio(
			posterior.from.2.tests(
				prevalence=prev, 	X1=NA),
	interval=c(0.05, 0.58)
	)$root

