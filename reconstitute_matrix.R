reconstitute = function(m) {
	if(! identical(m, t(m)))
		stop("Fails unless matrix is self-adjunct (symmetric if real)")
	N = dim(m)[1]
	eig.l = eigen(t(m))
	lv = eig.l$vectors
	eig.r = eigen(m)
	rv = eig.r$vectors
	eigenvalues = rbind(eig.l=eig.l$values,eig.r=eig.r$values)
	print(eigenvalues)
	components <- lapply(as.list(1:N), 
		function(i, eig.r, rv, lv) {  ### Extra args required by Splus:  different scoping rules.
				eig.r$values[i] * outer(rv[,i], lv[,i])}, 
			eig.r=eig.r, rv=rv, lv=lv)
	print(components)
	reconst.expression = paste ( paste("components[[", 1:N, "]]"), collapse=" + ")
	catn(reconst.expression)
	reconst = eval(parse(text=reconst.expression))
	catn("Input matrix")
	print(m)
	catn("Output matrix")
	reconst
}
	
m = matrix(c(.8,.2,.2,.8), nrow=2); reconstitute(m)   #OK
m = matrix(c(0,1,1,0), nrow=2); reconstitute(m)  #OK
m = matrix(c(0,1,-1,0), nrow=2); reconstitute(m)  #Wrong in Splus.  Correct in R.
m = matrix(c(	0,1,0,		
				1,0,0,
				0,0,1), nrow=3); reconstitute(m)  #OK
m = random.walk.matrix(n1=n1<-3);  reconstitute(m)  ### Not OK.
m = random.walk.matrix(n1=n1<-9);  reconstitute(m)
test.eigen.accuracy <- function(m) 
		max(abs(m%*%eigen(m)$vectors - eigen(m)$vectors * 				matrix(eigen(m)$values,nrow=n1+1,ncol=n1+1, byrow=T) ))
test.eigen.accuracy(m)
