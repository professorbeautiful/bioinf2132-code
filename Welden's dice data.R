#  Welden's dice data


n = 26306
xbar = 4.0524
x = round(n * 4.0524)

(xbar - 4)/sqrt(24/9/n)
### 5.20,  not 5.17.
z = (xbar - 4)/sqrt(2.6983/n)

sigma = sqrt(2.6983)

m0 = dnorm((xbar - 4)/sqrt(2.6983/n))

###    B ~ sqrt(n)* tau/sigma * exp(-z^2/2)

tau = sigma

sqrt(n)* tau/sigma * exp(-z^2/2)