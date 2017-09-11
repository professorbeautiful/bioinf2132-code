simple.bootstrap <- function(X, FUN, ..., B = 1000, seed = 0, block = 50){
	# demonstration program for nonparametric bootstrapping
	# X is a matrix or data frame, rows are observations
	# FUN(X, ...)
	# This version of bootstrap assumes that FUN() returns a scalar and
	# that B is a multiple of block
	# B bootstrap replications
	# seed is an integer between 0 and 1000
	# block is the block size, 
	#   the number of bootstrap values computed simultaneously

	set.seed(seed) # So results are reproducible

	if(is.null(dim(X))) X <- as.matrix(X)
	n <- nrow(X)

	call.stat <- function(i, X, FUN, indices, ...)
		FUN(X[indices[,i], ], ...)
	# The call.stat() function will be called by lapply() to
	# do the actual bootstrapping

	nblocks <- ceiling(B/block) # number of blocks
	result <- numeric(B) # Create space for results
	indices <- matrix( integer(n*block), nrow=n)
	temp <- 1:block
	names(temp) <- rep("",block) # This prevents lapply from coercing temp
	# to a list internally, which uses more memory.

	on.exit({
		cat("Saving replications 1:",
			(i-1)*block," to .bootstrap.results\n")
		assign(".bootstrap.results",replicates,where=1,immediate=T)
	}) # In case function is interrupted

	for(i in 1:nblocks){
		indices[] <- sample(1:n, n*block, T) # Sample the indices
		result[temp+block*(i-1)] <-
		unlist(lapply(temp, call.stat, X, FUN, indices, ...))
	}
	on.exit()
	result
}
