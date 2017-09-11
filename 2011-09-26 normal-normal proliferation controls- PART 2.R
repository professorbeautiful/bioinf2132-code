
##############   TODO:  Formulas for conditional mean and variance.  ###############

##############   TODO:  Use as intro to empirical Bayes.  ###############

##############   TODO:  what if the long-term means of the control subjects are NOT known beforehand?

###############   NOT FINISHED:  2011-10-02  ############
meanAlpha = 2
varAlpha = 10
meanAlphaBetaY = c(rep(meanAlpha,3), rep(0,nDays), rep(meanAlpha,3*nDays))
varBeta = priorVariance
varEpsilon = dataVariance
aNames = paste("a", 1:3, sep=".")
bNames = paste("b", 1:nDays, sep=".")
YNames = paste("Y", outer(1:3,1:nDays,paste,sep="."), sep=".")
eachdimname = c(aNames, bNames, YNames)

SigmaAlphaAlpha =  diag(varAlpha, 3)
dimnames(SigmaAlphaAlpha) = list(aNames, aNames)
SigmaAlphaBeta = matrix(0, nrow=3, ncol=nDays)
dimnames(SigmaAlphaBeta) = list(aNames, bNames)
SigmaAlphaY = kronecker( matrix(1, 1, nDays), diag(varAlpha,3))
dimnames(SigmaAlphaY) = list(aNames, YNames)
SigmaBetaBeta = diag(varBeta, nDays)
dimnames(SigmaBetaBeta) = list(bNames, bNames)
SigmaBetaY = kronecker(diag(varBeta, nDays), matrix(1, 1, 3))
dimnames(SigmaBetaY) = list(bNames, YNames)
SigmaYY = (
  kronecker(matrix(1, nDays, nDays), SigmaAlphaAlpha) 
  + kronecker(SigmaBetaBeta, matrix(1, 3, 3)) 
  + varEpsilon
)
dimnames(SigmaYY) = list(YNames, YNames)

varAlphaBetaY= rbind(
  cbind(
    SigmaAlphaAlpha,
    SigmaAlphaBeta,
    SigmaAlphaY),
  cbind(
    t(SigmaAlphaBeta),
    SigmaBetaBeta,
    SigmaBetaY),
  cbind(
    t(SigmaAlphaY),
    t(SigmaBetaY),
    SigmaYY)
)
dim(varAlphaBetaY)   ###   3 + 30 + 90
dimnames(varAlphaBetaY) = list(eachdimname, eachdimname)

###    what is  [Alpha, Beta | Y ]   ??  #################
###   Apply formula, with Y1=Y, and Y2=(Alpha,Beta)
### So...
..Y1 = controlData$Y  ###  check:  correct order?  yes.
indicesAlphaBeta = 1:(3+nDays)
indicesY = (-1)* indicesAlphaBeta
..mu1 = meanAlphaBetaY[indicesY]
..mu2 = meanAlphaBetaY[indicesAlphaBeta]
..SIG12 = varAlphaBetaY[indicesY, indicesAlphaBeta]
..SIG22 = varAlphaBetaY[indicesAlphaBeta, indicesAlphaBeta]
..SIG11 = varAlphaBetaY[indicesY, indicesY]
posteriorMeanAlphaBeta = 
  ..mu2 + t(..SIG12) %*% solve(..SIG11) %*% (..Y1 - ..mu1)
########   nearly SINGULAR, regardless of varAlpha --  
posteriorMeanAlphaBeta = 
  ..mu2 + t(..SIG12) %*% solve(..SIG11 + diag(1e-10, 3*nDays)) %*% (..Y1 - ..mu1)
posteriorMeanAlpha = posteriorMeanAlphaBeta  [1:3]
posteriorMeanBeta = posteriorMeanAlphaBeta  [4:(3+nDays)]
plot(trueDailyEffects, posteriorMeanBeta)
diagonal()

###### posterior  covariance matrix  ##########
posteriorVarianceAlphaBeta = 
  ..SIG22  -
  t(..SIG12) %*% solve(..SIG11+ diag(1e-10, 3*nDays)) %*% ..SIG12
posteriorSDBeta = sqrt(diag(posteriorVarianceAlphaBeta))   [-(1:3)]
for(day in 1:nDays) lines(rep(trueDailyEffects[day], 2), 
                          posteriorMeanBeta[day]  + c(-1,1)* posteriorSDBeta[day])
######     These variances don't seem right.
####     lme4()