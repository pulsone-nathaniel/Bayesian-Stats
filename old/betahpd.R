# betahpd.R by Nathaniel Pulsone
# returns the 100×(1−α)% highest HPD region for a beta(a, b) distribution
library(coda)
set.seed(123) # for reproduceability with the interval check
betahpd = function(alpha = 0.05, a = 1, b = 1){
  # to check accuracy of the interval
  samples = as.mcmc(rbeta(10000, a, b))
  result1 = HPDinterval(samples, 1-alpha)
  
  lower = (a-1)/(a+b-2)
  upper = lower
  d = dbeta(lower, a, b)
  r = d
  coverage = 0
  while(coverage < 1-alpha){
    lower = lower - .0001
    d = dbeta(lower, a, b)
    r = dbeta(upper, a, b)
    while(r > d){
      upper = upper + .00001
      r = dbeta(upper, a, b)
    }
    coverage = pbeta(upper, a, b) - pbeta(lower, a, b)
  }
  result2 = cbind(lower, upper)
  
  return(list("MC" = result1,
              "NU" = result2))
}
