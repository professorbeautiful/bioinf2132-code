Q1 = 10
Q3 = 50

f = function(p) {
  a = p[1]
  b = p[2]
  value = (pgamma(Q1, a, b) - 0.25)^2 +
    (pgamma(Q3, a, b) - 0.75)^2
  value
}
nlm(f, p=c(a=4, b=10), fscale=0, print.level=1)
