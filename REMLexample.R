###   https://www.r-bloggers.com/exercise-in-remlmixed-model/

r1 <- read.table(textConnection('
                                1 1 1 56 1 1 2 41
                                1 2 1 50 1 2 2 36
                                1 3 1 39 1 3 2 35
                                2 1 1 30 2 1 2 25
                                2 2 1 36 2 2 2 28
                                2 3 1 33 2 3 2 30
                                3 1 1 32 3 1 2 24
                                3 2 1 31 3 2 2 27
                                3 3 1 15 3 3 2 19
                                4 1 1 30 4 1 2 25
                                4 2 1 35 4 2 2 30
                                4 3 1 17 4 3 2 18'),
                 
                 col.names=c('Block1','A1','B1','Y1','Block2','A2','B2','Y2'))

sp <- with(r1,data.frame(
  Block=factor(c(Block1,Block2)),
  A=factor(c(A1,A2)),
  B=factor(c(B1,B2)),
  Y=c(Y1,Y2)))

require(lme4)
l1 <- lmer(
  Y ~ A + B + A*B + (1 |Block) + (   1| A : Block) ,
  data=sp)
summary(l1)
l2 <- lmer(
  Y ~ A + B + (1 |Block) + (   1| A : Block) ,
  data=sp)
anova(l2,l1)

##  We check the results.
##  This will be a PARAMETRIC 'bootstrap',
##   not the classic resampling bootstrap.
##  Data are generated out of the null model.
library(boot)
bootFun <- function(m0,m1) {
  s <- simulate(m0)
  L0 <- logLik(refit(m0,s), REML=TRUE)
  L1 <- logLik(refit(m1,s), REML=TRUE)
  2*(L1-L0)
}
# Observed deviance:
obsdev <- 2*( as.numeric(logLik(l1, REML=TRUE))-as.numeric(logLik(l2, REML=TRUE)))
set.seed(1001)
reps <- replicate(1000, bootFun(l2,l1))
cat("Observed deviance is ", obsdev, "\n")
cat("P value relative to the bootstrap distribution is",
    mean(reps>obsdev), "\n")
hist(reps)
abline(v=obsdev, col="red") 
p.vec=seq(0,1, length=1000+2)[-c(1,1002)]
plot(qchisq(p.vec, df = 2), sort(reps))
abline(a=0, b=1)
abline(a=-5, b=1)

nCopies = 2
spBig = Reduce(rbind, x=
                 lapply(
  1:nCopies, function(ignore){data.frame(sp[,-4], simulate(l2))} ) )
names(spBig)[4] = "Y"
dim(spBig); names(spBig)
l1Big <- lmer(
  Y ~ A + B + A*B + (1 |Block) + (   1| A : Block) ,
  data=spBig)
l2Big <- lmer(
  Y ~ A + B + (1 |Block) + (   1| A : Block) ,
  data=spBig)
reps <- replicate(1000, bootFun(l2Big,l1Big))
plot(qchisq(p.vec, df = 2), sort(reps))
title(paste("sample size = ", nrow(spBig)))
abline(a=0, b=1)
abline(a=-5, b=1)


####
mcs <- mcmcsamp(l1,10000)
mcsdf <- as.data.frame(mcs)
c(mean=mean(mcsdf[,1]+.5*mcsdf[,'B2f']),sd=sd(mcsdf[,1]+.5*mcsdf[,'B2']))
library(MCMCglmm)
m1 <- MCMCglmm(Y ~ A + B + A*B, random= ~Block + A : Block ,
               data=sp,family='gaussian',nitt=500000,thin=20,
               burnin=50000,verbose=FALSE)
summary(m1)

#Test for fixed effects
table(sign(m1$Sol[,'A2:B2']),sign(m1$Sol[,'A3:B2']))/nrow(m1$Sol)

# mean for a level 1
c(mean=mean(m1$Sol[,1]+.5*m1$Sol[,'B2']),
  sd=sd(m1$Sol[,1]+.5*m1$Sol[,'B2']))

## MCMC variances
lapply(1:3,function(i) quantile(m1$VCV[,i],seq(.1,.9,.2)))

l3 <- lmer(Y ~ A + B + A*B + (1 |Block)  ,data=sp)
summary(l3)
obsdev13 <- 2*( as.numeric(logLik(l1))-as.numeric(logLik(l3)))
reps13 <- replicate(1000,pboot(l3,l1))
mean(reps13>obsdev13)



