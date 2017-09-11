###     2011-09-26   normal-normal proliferation controls example
###  Presented  2014-01-22

#png("dailyControls-raw-data.png", width=960)
###  Regarding graphics:   I find the best graphics into a document is to plot into the native plot window. adjust the way you want, then use a screen grabber.  In windows, Faststone capture is great.  In Mac use cmd-control-shift-4 to put selected area into pastebuffer.

nDays = 30
nControls = 3
controlMeans = 1:nControls   ###  alpha 
dataVariance = 1  					###tau
priorMean = 0							### m
priorVariance = 0.2					### V
controlData = expand.grid(control=1:nControls, day=1:nDays)
trueDailyEffects = rnorm(nDays, priorMean, sqrt(priorVariance))
controlData$trueDailyEffect = trueDailyEffects[controlData$day]
controlData$sampleEffect = controlMeans[controlData$control]
controlData$Y = controlMeans[controlData$control] + controlData$trueDailyEffect + rnorm(1:nrow(controlData), 0, sqrt(dataVariance))

with(controlData,	plot(day, Y))
for(control in 1:nControls) with(controlData[controlData$control == control, ],
			lines(day, Y, lty=2))
aggregate(data=controlData, Y ~ control, FUN=mean)
dailyMeans = print(aggregate(data=controlData, Y ~ day, FUN=mean))$Y
lines(1:nDays, dailyMeans, lwd=2, col="red")
title("Raw data and raw daily means")
#dev.off()

#jpeg("dailyControls-dailyEffects.jpeg", width=960)
###  We suppose that we know controlMeans, so we subtract that off to get raw daily effects.
estimatedDailyEffects = dailyMeans - mean(controlMeans)
plot(1:nDays, estimatedDailyEffects, lwd=1, ty="b", pch="R", col="red")
###  Compare to the truth:
points(1:nDays, trueDailyEffects, lwd=1, pch="T")

###  Prior means for dailyEffects all equal zero:
lines(1:nDays, rep(priorMean,nDays), lwd=3, col="blue", lty=2)

###  Finally, the Bayesian posterior means:
shrinkerB = priorVariance/(priorVariance + dataVariance/nControls)
posteriorMeans = estimatedDailyEffects*shrinkerB + priorMean*(1-shrinkerB)
points(1:nDays, posteriorMeans, pch="B", ty="b", col="blue")

title("R=raw means\n T= true day effects\n B=Bayes estimates (dotted=prior)")

summary((estimatedDailyEffects - trueDailyEffects)^2)
summary((posteriorMeans - trueDailyEffects)^2)
mean(trueDailyEffects^2)
priorVariance*(1-shrinkerB)  ##  posterior (conditional) variance

###  With 3000 days, you can see that the means of the squared errors converge to theoretical.

plot(estimatedDailyEffects - trueDailyEffects, posteriorMeans - trueDailyEffects)
abline(a=0, b=1, lwd=3)

shrinkerVector = seq(0, 1, 0.02)
mseFunction = function(w) 
  mean(
    (estimatedDailyEffects*w + priorMean*(1-w) - trueDailyEffects)^2
  ) 
mseVector = sapply(shrinkerVector, mseFunction)

plot(shrinkerVector, mseVector, type="l", xlab="shrinker", ylab="mean(squared error)", ylim=c(0, 0.35))

points(x=0, y=mseVector[1], pch=".", cex=10)
text(0, mseVector[1] + 0.03, expression(hat(beta) == zero), pos=4)
text(0, mseVector[1] , "(prior)", pos=4)

points(x=1, y=mseVector[length(mseVector)], pch=".", cex=10)
text(1  - 0.15, y=mseVector[length(mseVector)]+0.05, 
     expression(hat(beta) == X(MLE)), pos=1)

points(x=shrinkerB, y=mseFunction(shrinkerB), pch=".", cex=10,
       col="red")
text(col="red", adj=c(0,0),
     shrinkerB, mseFunction(shrinkerB)+0.05, 
     expression(hat(beta)== a[B](X)))
text(col="red", adj=c(0,0),
     shrinkerB+0.05, mseFunction(shrinkerB)-0.02,  
     "Bayes rule")
# text(col="red", adj=c(1,1),
#      shrinkerB, mseFunction(shrinkerB)-0.09, "(a[B])"))  ### roughly equivalent to expression()
#abline(col="red", v=shrinkerB)

whichAppearsBest = which(mseVector == min(mseVector))[1]
points(x=shrinkerVector[whichAppearsBest], 
       y=mseVector[whichAppearsBest], pch=".", cex=10,
       col="darkgreen")
text(shrinkerVector[whichAppearsBest], 
     mseVector[whichAppearsBest],
     expression(hat(beta) == delta[hat(B)] ), adj=c(1,1),
     col="darkgreen")
text(col="darkgreen", adj=c(1/2,1),
     shrinkerVector[whichAppearsBest],
     mseVector[whichAppearsBest] - 0.05,
     "estimated Bayes rule")


#### For the "plotmath" stuff to work, 3rd arg  has to be parsable as an R expression.
library(mvbutils)
title("mean over " %&%  nDays %&% " days")
