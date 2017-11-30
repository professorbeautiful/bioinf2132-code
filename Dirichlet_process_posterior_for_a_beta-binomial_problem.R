### Dirichlet process posterior for a beta-binomial problem.

library(DPpackage)

data(rolling)
y <- cbind(rolling$y1,rolling$y2)

table(y)

# Prior information

prior<-list(alpha=1,
            a1=1,
            b1=1)

# Initial state
state <- NULL

# MCMC parameters

mcmc <- list(nburn=5000,
             nsave=10000,
             nskip=3,
             ndisplay=100)

# Fitting the model

fit <- DPbetabinom(y=y,ngrid=100, 
                   prior=prior, 
                   mcmc=mcmc, 
                   state=state, 
                   status=TRUE)
str(fit$state)
length(unique(fit$state$p))
table(fit$state$ss)
table(fit$state$ss, y[ , 1])

fit
summary(fit)
names(fit)
str(fit$save.state)


# density estimate
plot(fit,output="density")
points((1:9)/9, table(y)/100, col='red', type='h')

# parameters
plot(fit,output="param")

