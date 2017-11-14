###  MULTIPLE IMPUTATION ######	
###   The dataset is generated with the model logit(E(Z)) = a + Xb 
###    with missing values for some X's  (nonresponse).

N = 1000
P = 5  ### number of predictors
beta = rep(0.5,P)  ### regression parameter
alpha = 0 			### intercept parameter

correlation = 0.95  ### correlations of predictors
predictorCovar = 
	correlation*matrix(1, nrow=P, ncol=P) +
	(1-correlation)*diag(1, P)

### Create the design matrix.
require(mvtnorm)
predictorMatrix <- rmvnorm(N, sigma=predictorCovar)
dim(predictorMatrix)
linearPredictor = alpha + predictorMatrix %*% beta
### Generate binary the outcomes via logistic model.
outcomes <- rbinom(N, 1, antilogit(linearPredictor))
### Complete dataset.
datasetComplete = data.frame(
	Z=outcomes,
	as.data.frame(predictorMatrix)
)
####  A function for doing logistic regression ###
analyzeAdataset = function(datasetName, P, print=FALSE, ...) {
	###  construct formula
	formulaString = paste("V", 1:P, sep="",collapse="+")
	glmString = paste(sep="",  "glm(Z ~ ", formulaString,
		 ", data=", datasetName, ",
			family=binomial, ...)")
	if(print) catn("glmString is " %&% glmString)
	glmExpression <- parse(text=glmString)
	if(print) print(glmExpression)
	glm.out = eval(glmExpression)
	return(glm.out)
}
catn = function(...) cat(..., "\n")
`%&%` = function(a, b) paste(a, b, sep="")
glm.out.complete = analyzeAdataset("datasetComplete", P, print=TRUE)
summary(glm.out.complete)
anova(glm.out.complete)
step(glm.out.complete)

#######  But now many observations are missing.
datasetObserved = datasetComplete
probabilityMissing = 0.8
isItMissing = rbinom(N, 1, probabilityMissing) ## At Random!
whichToCensor = sample(1:P, N, replace=TRUE)
for(i in 1:N) 
	if(isItMissing[i])
		datasetObserved[i, -1] [ , whichToCensor[i] ] = NA
datasetObserved[1:6, ]
colMeans(is.na(datasetObserved))

glm.out.observed.omit = analyzeAdataset("datasetObserved", P, T, na.action=na.omit)
summary(glm.out.complete)$coef
summary(glm.out.observed.omit)$coef

#####   Imputation  ########
mean.estimated = apply(datasetObserved[,-1], 2, mean, na.rm=TRUE)
cov.estimated = cov(datasetObserved[ , -1],
	 use="pairwise")
cov.estimated;  cov.estimated - predictorCovar ;   summary(c(cov.estimated - predictorCovar))

### Impute the missing values via modeling ###
imputeAdataset = function() {
	#### NOTE: this imputation IGNORES Z.  Should it???
	####  If not what should be done?
	datasetCompleteImputed = datasetObserved
	predictorsObserved = datasetObserved[,-1]
	for(i_ in 1:P) {
		#This refers to formulas for the mean and variance of
		#   [Vi_ |  (..., V[i_-1], V[i_+1], ...)]
		###   "1"  means Vi_.		"2"  means all except Vi_.
		Sigma12 = cov.estimated[i_, -i_]  # row vector
		Sigma21 = cov.estimated[-i_, i_]  # column vector
		Sigma22 = cov.estimated[-i_, -i_]
		i_IsMissing = which(is.na(predictorsObserved[,i_]))
		Z2 = predictorsObserved[i_IsMissing, -i_] 
		Z2 = as.matrix(Z2)
		meanZ2 = mean.estimated[-i_]
		conditionalMean = mean.estimated[i_] + 
			Sigma12 %*% solve(Sigma22) %*% 
			t(as.matrix(sweep(Z2, 2, meanZ2, "-")))
		conditionalVariance = cov.estimated[i_, i_] -
			Sigma12 %*% solve(Sigma22) %*% Sigma21
		datasetCompleteImputed[,-1][i_IsMissing, i_] = 
			rnorm(length(i_IsMissing),
				mean=conditionalMean,
				sd=sqrt(conditionalVariance))
	}
	return(datasetCompleteImputed)
}

###  So let's compare the imputed values to the true values.  ####
datasetCompleteImputed = imputeAdataset()
head(is.na(datasetObserved))	
c(head(is.na(datasetObserved)))	
head(unlist(datasetComplete))
head(unlist(datasetComplete)[c(is.na(datasetObserved))])
head(unlist(datasetObserved)[c(is.na(datasetObserved))])
head(unlist(datasetCompleteImputed)[c(is.na(datasetObserved))])
plot(x=unlist(datasetComplete)[c(is.na(datasetObserved))],
	xlab="true values, now missing",
	y=unlist(datasetCompleteImputed)[c(is.na(datasetObserved))],
	ylab="imputed values")
diagonal = function() abline(a=0, b=1)
diagonal()

#######  MULTIPLE IMPUTATION  happens here.              ########
mImpute = 10
imputations = lapply(1:mImpute, function(ignoreMe) {
		datasetCompleteImputed <<- imputeAdataset()
		analyzeAdataset("datasetCompleteImputed", P)
	}
)

###  Estimating the parameters from the imputations.
estimatesImputedAll = sapply(imputations, function(imp){
	summary(imp)$coef[,"Estimate"]
})
estimatesImputed = apply(estimatesImputedAll, 1, mean)
print(estimatesImputed)
print(summary(glm.out.complete)$coef[, "Estimate"])
###  Multiple imputation variance estimate  ######################
extractVariances = function(imp)
  summary(imp)$coef[,"Std. Error"]^2
variancesImputedAll = sapply(imputations, extractVariances)
variancesWithinImputed = apply(variancesImputedAll, 1, mean)
variancesBetweenImputed = apply(estimatesImputedAll, 1, var)
variancesComplete = extractVariances(glm.out.complete)
variancesTotal = variancesWithinImputed + variancesBetweenImputed
rbind(within= variancesWithinImputed, 
      between= variancesBetweenImputed,
      total= variancesTotal,
      complete= variancesComplete)

plot(variancesComplete, 
     variancesTotal/variancesComplete,
        pch=c("I", as.character(1:5)))
#	  summary(glm.out.complete)$coef[ , 'z value']^2

