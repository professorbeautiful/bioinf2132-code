logit = function(p) log(p/(1-p))
antilogit = function(x)  1-1/(1+exp(x))

DLdata = matrix(c(3,5,2,90),nrow=2)
dimnames(DLdata) = list(c("dark","light"),c("R","N"))
logit.hat = logit(apply(DLdata, 1, function(r)r[1]/sum(r)))


### Split:

logML.Split = lgamma(sum(DLdata)+1)-sum(lgamma(DLdata+1)) +
	sum(DLdata*(log(DLdata/sum(DLdata))))
#AIC

AIC.Split = 2*logML.Split - 2*3
BIC.Split = 2*logML.Split - log(sum(DLdata))*3
BF.Split = exp(BIC.Split/2)

### Lump:

collapsedData = apply(DLdata,2,sum)
logML.Lump = log(dbinom(sum(DLdata[1,]), size=sum(DLdata), p=sum(DLdata[1,])/sum(DLdata))) + 
			log(dbinom(sum(DLdata[,1]), size=sum(DLdata), p=sum(DLdata[,1])/sum(DLdata))) + 
			log(phyper(DLdata[1,1], sum(DLdata[1,]), sum(DLdata[2,]), sum(DLdata[,1])))
	
	
#AIC

AIC.Lump = 2*logML.Lump - 2*2
BIC.Lump = 2*logML.Lump - log(sum(DLdata))*2
BF.Lump = exp(BIC.Lump/2)

BF = BF.Lump/BF.Split


results = matrix(nrow=2, byrow=T,
	dimnames=list(c("split", "lump"),  c("LL", "AIC", "BIC", "BF")),
	c(
	logML.Split, AIC.Split, BIC.Split, BF.Split,
	logML.Lump, AIC.Lump, BIC.Lump, BF.Lump)
)

comparison = c(results[1,1:3] - results[2,1:3])
comparison = c(comparison, results[1,4] / results[2,4])
names(comparison) =  words("LLdiff AICdiff BICdiff BCratio")
print(results)
print(comparison)
