####  NOTE three exercises A, B, and C, below.

####  See class document 2011-09-06 Diagnosis of coronary artery disease with nuisance parameter.
####  Data from Weiner et al 1979.
####  EST (exercise stress test)  and EKG (electrocardiogram)
####     as predictors of CAD (coronary artery disease)
####  We treat the probabilities estimated from the data as if there were no noise
#### .... for now.

twoTestData.original = array(
	c(	   
	700,221,139,467,		 
	238,33,141,106, 
	938,254,280,573		 
	)
	,
	dim=c(2,2,2),
	dimnames=list(
	EST=c("EST+", "EST-" ),	
	disease=c("CAD","noCAD"),
	EKG=c("EKG-","EKG+" )
	)	
)

twoTestData = aperm(twoTestData.original, c(3,1,2))  ### Permute array dimensions. 
###  Rearranged now to match the table in the notes.

prior = apply(twoTestData, 3, sum)/sum(twoTestData)  ### The prevalence marginal distribution.

prevalenceByTestResults =  
  sweep(twoTestData, 1:2, apply(twoTestData, 1:2, sum), "/")
###  Compare with "Prevalence" table in our class document

CADmodel =    ###  Condition on disease status to get the model(s).
  sweep(twoTestData, 3, apply(twoTestData, 3, sum), "/")
Sensitivities = CADmodel[ , "EST+", "CAD"]/
 	(CADmodel[ , "EST+", "CAD"]+CADmodel[ , "EST-", "CAD"])
Specificities = CADmodel[ , "EST-", "noCAD"]/
 	(CADmodel[ , "EST+", "noCAD"]+CADmodel[ , "EST-", "noCAD"])
### Compare with "Model" table in our class document
 ##  Probability of test results given CAD status.
apply(CADmodel, 3, sum)  ### confirm, each model adds to 1.

priorArray = CADmodel  ### Just to get ther right dimensions and dimnames.
priorArray[ , , ] = rep(prior, each=4)  ### to multiply prior times models.

priorArray * CADmodel  #### prior times models.  Same as twoTestData after normalization.
sum(priorArray * CADmodel)  ### = 1.

###  The data marginal probability,  [X] = [X1,X2]  (for each combination outcome)
dataMargin = apply(twoTestData, 1:2, sum)/sum(twoTestData.original)
###  To get the posterior prob(CAD), divide by [X]
(priorArray * CADmodel)[ , , "CAD"] / dataMargin    #### 4 posteriors, one for each data outcome.


#### EXERCISE A:   COMPLETE THIS FUNCTION

### DOCUMENTATION:   in the style of Roxygen
#'	@param X1 The first test result. May be NA. Default="EKG+".
#'	@param X2 The second test result. May be NA. Default="EST+".
#'	@param prevalence  The prevalence of disease. Default=0.58.
#' @param model This is condit'l distr of all test results, as an array. The first two dimensions are tests with binary outcomes.  The third dimension is the disease state dimension; the first disease state corresponds to "disease present", the second to "disease absent". Default=CADmodel.
#' @return The posterior probability of disease given the data.
posterior.from.2.tests = function(
	X1="EKG+",
	X2="EST+",	prevalence=0.58,		
	model=CADmodel
) {
	### FILL IN THE BODY OF THIS FUNCTION  ##########################
	### TEST YOUR FUNCTION !! 								 ##########################
}

posterior.from.2.tests()
posterior.from.2.tests(X1="EKG-")
posterior.from.2.tests(X1=NA)
posterior.from.2.tests(X2=NA)
posterior.from.2.tests(X2="EST-", prevalence=seq(0.1,0.9,0.1))



#### EXERCISE B:   complete the body of this function
#'	@param posterior.  Prior or current (or imagined future) posterior probability of disease.   
#'	@param lossRatio.  Loss("treat", "no disease")	/ Loss("wait", "disease")	. Default=1.
#' @return  Optimal action.
BayesRule.binaryAction.binaryState = function(
	posterior,
	lossRatio=1
) {
	### FILL IN THE BODY OF THIS FUNCTION  ##########################
	### TEST YOUR FUNCTION !! 								 ##########################
}
BayesRule.binaryAction.binaryState(.5)
BayesRule.binaryAction.binaryState(seq(0,1,.1))
BayesRule.binaryAction.binaryState(posterior.from.2.tests(X1="EKG-", X2="EST-", prevalence=seq(0,1,.1)))

#### EXERCISE C:   write a function with inputs of the form	
##    the disease prevalence, 
##   the model for [EKG | disease],  and the model for [EST | EKG, disease state],   
##   and returns the probability of disease conditional on EST result, marginal over EKG.
##  Provide appropriate default values.
## Test your function.

###  As part of this exercise, run this command:
apply(CADmodel, c(1,3), sum)
###   and discuss what it says about the design of the study.
########  Note that a real-life [EKG | disease] may differ considerably from the data set.



