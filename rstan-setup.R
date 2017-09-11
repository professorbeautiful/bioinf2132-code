### https://github.com/stan-dev/rstan/wiki/RStan-Mac-OS-X-Prerequisite-Installation-Instructions
### https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started#how-to-use-rstan

# dotR <- file.path(Sys.getenv("HOME"), ".R")
# if (!file.exists(dotR)) dir.create(dotR)
# M <- file.path(dotR, "Makevars")
# if (!file.exists(M)) file.create(M)
# cat("\nCXXFLAGS=-O3", file = M, sep = "\n", append = TRUE)
# cat("\nCC=clang", "CXX=clang++ -arch x86_64 -ftemplate-depth-256", 
#     file = M, sep = "\n", append = TRUE)
# Sys.setenv(MAKEFLAGS = "-j4") 
# install.packages("rstan", dependencies = TRUE)

schools_dat <- list(J = 8, 
                    y = c(28,  8, -3,  7, -1,  1, 18, 12),
                    sigma = c(15, 10, 16, 11,  9, 11, 10, 18))
library(rstan)
fit <- stan(file = '8schools.stan', data = schools_dat, 
            iter = 1000, chains = 4)
print(fit)
plot(fit, pars=c("mu", "tau"))
plot(fit)
pairs(fit, pars = c("mu", "tau", "lp__"))
### lp__  is the log posterior

la <- extract(fit, permuted = TRUE) # return a list of arrays 
mu <- la$mu 
a <- extract(fit, permuted = FALSE) 
str(a)
a2 <- as.array(fit)  ### same as a
m <- as.matrix(fit)

print(fit, digits = 3)
