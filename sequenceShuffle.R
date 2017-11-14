#### MCMC example
###   like ClustalW, for aligning sequences

seqData = lapply(6:13, function(L) 
{
  selections = c(rmultinom(1, L, rep(1/8,8)))
  sequence = sample(rep(letters[1:8], selections))
  sequence[(L-6+1):(L-6+3)] = sample(c("X", "Y", "Y"))
  sequence
}
)
seqData
class(seqData) = "seqData"

W = 3

maxLength = max(sapply(seqData, length))
print.seqData = function(s) {
  cat(c(rep("_", maxLength - W ),
        rep("*", W), "\n"))
  for(i in 1:length(s)) {
    cat(rep("_", maxLength - W + 1 - currentPos[i]),
        s[[i]], "\n")
  }
}
seqData


library("rBeta2009")  ## for random dirichlet

iterOnce = function(seqData=seqData, iterNum=1, printIt = FALSE) {
  zBlock <<- t( sapply(1:length(seqData), function(eachSeq)
    seqData[[eachSeq]]
    [currentPos[eachSeq]: 
        (currentPos[eachSeq]+W-1)])
  )
  #print(zBlock)
  ###  zBlock is z!
  ### currentPos is theta!
  for(nSeq in 1:length(seqData)) {
    blockWithRowOmitted = zBlock[-nSeq, ]
    newProb = matrix(NA, ncol=W, nrow=nalphabet)
    for(iCol in 1:W) {
      frequencies = table(blockWithRowOmitted[, iCol])
      parameters = background
      parameters[names(frequencies)] = parameters[names(frequencies)] + 
        frequencies[names(frequencies)]
      newProbThisCol = c(rdirichlet(1, parameters))
      names(newProbThisCol) =  names(background)
      newProb[ , iCol] = newProbThisCol
    }
    rownames(newProb) = alphabet
    positionProbs = sapply( 1:(length(seqData[[nSeq]]) - W + 1),
                            function(theta_l)
                              prod(sapply(1:W, function(w)
                                newProb[ seqData[[nSeq]][w + theta_l - 1],
                                         w])
                              )
    )
    positionProbs = positionProbs/sum(positionProbs)
    currentPos[nSeq] <<- which(rmultinom(1, size=1, positionProbs)[,1]==1)
    
  }
  if(printIt) {cat("\n"); print(seqData)}
  currentPos
}

###   Uncomment this to see bimodality.
seqData[[8]][1:3]= c('Y','Y','Y');  currentPos = rep(1, 8)

doClustalWchain = function(seqData=seqData, nIters=1000) {
  startingPos = rep(1, length(seqData))
  currentPos <<- startingPos
  alphabet = unique(c(unlist(seqData)))
  alphabet = sort(alphabet)
  nalphabet = length(alphabet)
  background = rep(1, nalphabet)
  names(background) = alphabet
  for( iIter in 1:nIters) finalPos <<- iterOnce(seqData)
  par(mfrow=c(2,4))
  for(seqNum in 1:length(seqData)) 
    plot(1:nIters , result[seqNum, ], 
         xlab="iteration", ylab="currentPos",
         main=paste("sequence ", seqNum),
         col=1+seqNum, type="l", ylim=c(1,maxLength-W+1))
  distributions = matrix(0, nrow=length(seqData), ncol=maxLength-W+1)
  sapply(1:length(seqData), function(seqNum) {
    theTable = table(result[seqNum, ])
    distributions[seqNum, as.numeric(names(theTable)) ] <<- 
      theTable
  }
  )
  print(distributions)
  currentPos <<- finalPos
  print(currentPos)
  cat('Aligned to the final position\n')
  print(seqData)
  modePos = sapply(1:nrow(distributions), 
                   function(n) { 
                     thisRow = distributions[ n, ]; 
                     which(max(thisRow)==thisRow)  
                   })
  currentPos <<- modePos
  print(currentPos)
  cat('Aligned to the modes\n')
  print(seqData)
  invisible(seqData)
}

doClustalWchain(seqData)
