### To use this file, 
### you need to install.packages('R2OpenBUGS'),
### and install the OpenBUGS program.
### That is not hard on Windows and Linux,
### but challenging on OSX.
### If you need it, I can give you help.

library(R2OpenBUGS)
library(coda)
data(schools)
schools.model = "
model {
       for (j in 1:J)
{
  y[j] ~ dnorm (theta[j], tau.y[j])
  theta[j] ~ dnorm (mu.theta, tau.theta)
  tau.y[j] <- pow(sigma.y[j], -2)
}
mu.theta ~ dnorm (0.0, 1.0E-6)
tau.theta <- pow(sigma.theta, -2)
sigma.theta ~ dunif (0, 1000)
}
"
write(schools.model, file = 'schools.model')
file.show("schools.model")

J <- nrow(schools)
y <- schools$estimate
sigma.y <- schools$sd
data <- list ("J", "y", "sigma.y")
inits <- function(){
  list(theta = rnorm(J, 0, 100), mu.theta = rnorm(1, 0, 100), 
       sigma.theta = runif(1, 0, 100))
}
OpenBUGS.dir=
  "/Users/Roger/.wine/drive_c/Program Files (x86)/OpenBUGS/OpenBUGS323"
OpenBUGS.pgm=paste0(OpenBUGS.dir, "/OpenBUGS.exe")
#OpenBUGS.pgm="/Users/Roger/Applications/Wineskin/OpenBUGS.app"

schools.sim <- bugs(data, inits, 
                    model.file = "schools.model",
                    parameters = c("theta", "mu.theta", "sigma.theta"),
                    n.chains = 3, 
                    n.iter = 1000,
                    OpenBUGS.pgm = OpenBUGS.pgm,
                    useWINE = TRUE,
                    newWINE = TRUE
                    #,debug=TRUE   #   Leads to error
                    #arguments 'show.output.on.console', 'minimized' and 'invisible' are for Windows only
                    )
## You can ignore the error err:ole:CoReleaseMarshalData
schools.sim
plot(schools.sim)
names(schools.sim)
str(schools.sim)

attach.bugs(schools.sim)
# > ls(pos=2)
# [1] "deviance"    "mu.theta"    "n.sims"      "sigma.theta" "theta" 
print(mean(theta[,1] > theta[,3]))
detach.bugs()
print(mean(schools.sim$sims.list$theta[,1] > schools.sim$sims.list$theta[,3]))



