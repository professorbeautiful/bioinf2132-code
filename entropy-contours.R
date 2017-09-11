entropy = function(pi1, pi2) (-pi1*log(pi1) + -pi2*log(pi2))

p = seq(0,1,length=100)

pi1pi2 = expand.grid(p, p)
pi1pi2m = as.matrix(pi1pi2)
dim(pi1pi2m)
pi1pi2[sample(size=40, 1:nrow(pi1pi2)), ]
z = matrix(nrow=100, 
           apply(pi1pi2m, 1, 
                 function(row) entropy(row[1], row[2])))
contour(x=p, y=p, z=z)

which(z==max(z,na.rm=T))
pi1pi2[which(z==max(z,na.rm=T)), ]  ## 36/99. actual argmax is exp(-1)
points(pi1pi2[which(z==max(z,na.rm=T)), ],
       col="red", pch="M")
abline(a=1, b=-1)
points(0.5, 0.5,, col="green", pch=19, cex=2)
title("Contour plot of entropy", sub = "red: unconstrained, green: constrained")
