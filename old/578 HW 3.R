# MA 578 Bayesian Stats HW 3 
# Nathaniel Pulsone

set.seed(123) # For Monte Carlo Reproduceability

# Question 1:

s = 10000
x = runif(s, 0, 1) # Draw iid U(0, 1) samples
z = log(x / (1 - x)) # Apply inverse CDF transform

# Implement Integral Bounds
for(i in seq(1:length(z))){
  if(z[i] < -1){
    z[i] = 0
  }

}

fz = z^2/(exp(z))
est = sum(fz) / s # Monte Carlo Estimator for integral

sf2 = sum((fz - est)^2) / (s - 1) # Monte Carlo Standard Error
lower = est - qnorm(.025, lower.tail = F) * sqrt(sf2/s)
upper = est + qnorm(.025, lower.tail = F) * sqrt(sf2/s)
CI = cbind(lower, upper) # 95% Confidence Interval
print("Question 1 Answer:")
print(CI)

# Question 2

# (b)
rz = function(n){
  z = c()
  y = c()
  status = c()
  while(sum(status == "y") < n){
    u = runif(1)
    zst = rexp(1)
    b = exp(2 * zst) / (1 + exp(zst))^2
    
    z = c(z, zst)
    y = c(y, u * 2 * exp(-zst))
    
    if(u <= b){
      status = c(status, "y") # sample is accepted
    }
    else{
      status = c(status, "n") # sample is rejected
    }
  }
  results = list("z" = z, # samples
                 "y" = y, # U * Cq(z),
                 "status" = status) # rejected or accepted
  return(results)
}

rs30 = rz(30)
plot(rs30$z, rs30$y, col = as.factor(rs30$status))

# (c)
rs10k = rz(10000)
prob = mean(rs10k$z[which(rs10k$status == "y")] > 1)

sf = sum((rs10k$z[which(rs10k$status == "y")] - prob)^2) / 9999
lower = prob - qnorm(0.025, lower.tail = F) * sqrt(sf / 10000)
upper = prob + qnorm(0.025, lower.tail = F) * sqrt(sf / 10000)

print("Question 2(c) Answer:")
print(prob)
print(cbind(lower, upper))
plot(rs10k$z, rs10k$y, col = as.factor(rs10k$status))

# Question 3
# (a)

theta_A = rgamma(10000, 237, 20)
theta_B = rgamma(10000, 125, 14)

prob = mean(theta_B < theta_A)
prob

# (b)

n0 = seq(1:30)

# Summary statistics from HW 2 Data
sa = sum(c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6))
na = 10
sb = sum(c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7))
nb = 13


probs = c()
for(i in n0){
  theta_A = rgamma(10000, 120 + sa, 10 + na)
  theta_B = rgamma(10000, 12*i + sb, i + nb)
  probs = c(probs, mean(theta_B < theta_A))
}
plot(n0, probs)

# (c)
Ytil_A = rnbinom(10000, 120 + sa, (10 + na)/(11 + na))
Ytil_B = rnbinom(10000, 12 + sb, (1 + nb)/(2 + nb))
prob = mean(Ytil_B < Ytil_A)
prob


probs_pred = c()
for(i in n0){
  Ytil_A = rnbinom(10000, 120 + sa, (10 + na)/(11 + na))
  Ytil_B = rnbinom(10000, 12 + sb, (i + nb)/(1 + i + nb))
  probs_pred = c(probs_pred, mean(Ytil_B < Ytil_A))
}
plot(n0, probs_pred)