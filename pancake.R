pancake = 
function(answer) {
	### Purpose:  given a list of arrays, all of the same shape, construct an array of one greater dimension by piling the list elements together.  Set the dimnames appropriately.
	### Based on sapply.
    if(length(common.dim <- unique(lapply(answer, dim))) != 1L) 
    		stop("Elements should be arrays of the same shape")
    else {
        	theDim = c(common.dim[[1]], length(answer))  ### dim of new array.
        	answer = array(unlist(answer), dim=theDim)
        	###  Create dimnames for answer.
        	eachNames = dimnames(answer[[1]])
        	### Assumes all components have the same dimnames. TODO:  test this, throw error.
        	if(is.null(eachNames))
        		eachNames = lapply(dim(answer[[1]]), seq)
			answerNames = c(eachNames, list(names(answer)))
        	dimnames(answer) = answerNames
    }
    return(answer)
}
