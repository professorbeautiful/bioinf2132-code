logit = function(p) log(p/(1-p))
antilogit = function(x)  1-1/(1+exp(x))

DLdata = matrix(c(3,5,2,90),nrow=2)
dimnames(DLdata) = list(c("dark","light"),c("R","N"))
logit.hat = logit(apply(DLdata, 1, function(r)r[1]/sum(r)))


### Split: independent poissons

ML.Split = prod(dpois(DLdata, lambda=DLdata))

logML.Split = log(ML.Split)

df.Split = 4


##Lump: independent poissons for row totals, binomial for columns

rowTotals = apply(DLdata, 1, sum)
pHat = sum(DLdata[, "R"]) / sum(DLdata)
colTotals = apply(DLdata, 2, sum)
ML.Lump = prod(dpois(rowTotals, lambda=rowTotals)) * 
	prod(dbinom(DLdata[, "R"], colTotals, pHat))

logML.Lump = log(ML.Lump)

df.Lump = 3


results = data.frame(
	logML=c(logML.Split, logML.Lump), 
	df=c(df.Split,df.Lump))
rownames(results) = c("split", "lump")

#AIC

results$AIC = 
	2*results$logML - 2*results$df

results$BIC = 
	2*results$logML - log(sum(DLdata))* results$df


comparison = results[1, ] - results[2, ]
names(comparison) =  words("logML.diff df.diff AIC.diff BIC.diff")
print(results)
print(comparison)


BayesFactor.estimate = exp(comparison$BIC.diff/2)
## From class: exact Bayes factor is...

cat("Bayes factor estimate = ", BayesFactor.estimate,
	"   \nFrom class: exact Bayes factor = ",
	 1/0.05121
)
