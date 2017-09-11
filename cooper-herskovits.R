table(rbinom(1e5, 9, rbeta(1e5,1,1)))
## Uniform on 0,...,9.  good.

#install.packages("rBeta2009")
library("rBeta2009")

myRdir = function(r, shape){
  if(length(shape)==2) 
    return(c(p<-rbeta(r, shape[1], shape[2]), 1-p))
  return(rdirichlet(r, shape))
} # Necessary because of flaw in design of rdirichlet
K = 10  ## ri
N = 27  ## Nij
alphas = rep(1,K)
Nrandom = 1e5
D <- table(sapply(1:Nrandom, function(ignoreMe) paste(collapse=",",
  rmultinom(1, N, myRdir(1, alphas))
)))
#D
#summary(as.vector(normalize(D)))

# Again, uniform.
1/length(D)
factorial(N)/factorial(N+K-1)*factorial(K-1)
## The ratio is the coverage. 
## check!  (for small sets)

## Now let's look at C-H.
sortedNames=lapply(strsplit(names(D), split=","), sort)
valuesChar = pancake(  unique(sortedNames))
values = as.numeric.matrix(valuesChar)
valuesChar = apply(valuesChar, 2, paste, collapse=",")
valuesUngrouped = as.numeric.matrix(pancake(
  (sortedNames)
))
valueMap = sapply(1:length(D), function(col)
  match(paste(valuesUngrouped[,col], collapse=","), valuesChar))
Dunorder = sapply(1:length(valuesChar), function(v)
  sum(D[valueMap==v]))
normalize(Dunorder)
### C-H formula:
factorial(K-1)/factorial(N+K-1) * 
  apply(values, 2, function(x) prod(sapply(x, factorial)))
###  Does not add to 1.


outcomes=apply(FUN=as.numeric, MARGIN=1:2,
               as.matrix(pancake(strsplit(names(D), split=",")))
)

impurity = data.frame(outcomes = names(D),
          multinom=factorial(N) / apply(outcomes, 2, 
           function(nn) prod(factorial(nn))), 
     gini=apply(outcomes, 2, 
           function(nn) sum(nn/sum(nn) * (1-nn/sum(nn)))))
with(impurity, plot(multinom, gini,
     xlab="multinomial coef",
     ylab="gini index", log="x") )
title(paste0("K = ", K, " (ri), N = ", N, "  (Nij)"))

impurity.sorted = impurity[order(impurity$gini), ]
head(impurity.sorted)
tail(impurity.sorted)
##  with K=4, N=7, not quite monotone.
# Gini = how often a randomly chosen element 
# from the set would be incorrectly labeled 
# if it were randomly labeled according to the 
# distribution of labels in the subset =  sum( f * (1-f)).
