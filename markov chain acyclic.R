eps = 0.1; w = 1-eps
P = matrix(
	c(eps, w, 0, 
		0, eps, w, 
		w, 0, eps),
	byrow=TRUE, nrow=3
)
