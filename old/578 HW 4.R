# 578 HW 4.R by Nathaniel Pulsone

set.seed(123) # for reproduceability

# Question 1
# (a)

dat1 = read.table("school1.dat.txt")[,1]
dat2 = read.table("school2.dat.txt")[,1]
dat3 = read.table("school3.dat.txt")[,1]

# Function returns posterior means and CIs for an input data set
# s = Monte Carlo samples
Estimates = function(dat, s = 10000){
  # prior
  mu0 = 5
  sig20 = 4
  kap0 = 1
  nu0 = 2
  
  # sample data
  n = length(dat)
  ybar = sum(dat)/n
  s2 = sum((dat - ybar)^2)/(n-1)
  
  # posterior update
  mun = (kap0 / (kap0 + n)) * mu0 + (n/(kap0 + n))*ybar
  nun = nu0 + n 
  sig2n = (1/nun) * ((nu0*sig20) + (n - 1)*s2 + (kap0 * n / (kap0 + n) * (ybar-mu0)^2))
  
  # Monte Carlo estimation
  sig2_MC = 1 / (rgamma(s, nun/2, nun/2 * sig2n))
  theta_MC = c()
  for(i in seq(1:s)){
    theta_MC = c(theta_MC, rnorm(1, mun, sqrt(sig2_MC[i]/(kap0 + n))))
  }
  
  theta_hat = mean(theta_MC)
  sig2_hat = mean(sig2_MC)
  
  theta_CI = quantile(theta_MC, probs = c(0.025, 0.975))
  sig2_CI = quantile(sig2_MC, probs = c(0.025, 0.975))
  
  return(list("output" = rbind(c(theta_hat, theta_CI),
               c(sqrt(sig2_hat), sqrt(sig2_CI))), # part a estimates
              "theta_MC" = theta_MC, # MC samples for part b
              "sig2_MC" = sig2_MC)
         )
}

school1 = Estimates(dat1)
school2 = Estimates(dat2)
school3 = Estimates(dat3)

school1$output
school2$output
school3$output

# (b)

th1 = school1$theta_MC
th2 = school2$theta_MC
th3 = school3$theta_MC

mean((th2<th1)*(th3<th1))


# (c)
# function to generate predictive samples for each school
rpred = function(school){
  th = school$theta_MC
  sg = sqrt(school$sig2_MC)
  
  yt = c()
  for(i in seq(1:length(sg))){
    yt = c(yt, rnorm(1, th[i], sg[i]))
  }
  return(yt)
}

yt1 = rpred(school1)
yt2 = rpred(school2)
yt3 = rpred(school3)

mean((yt3<yt1)*(yt2<yt1))

# Question 2
sizes = c(1, 2, 4, 8, 16, 32)
probs = c()
for(j in seq(1,length(sizes))){
  # prior
  mu0 = 75
  sig20 = 100
  kap0 = sizes[j]
  nu0 = sizes[j]
  
  # data
  n = 16
  yba = 75.2
  s2a = 7.3^2
  ybb = 77.5
  s2b = 8.1^2
  
  # updates
  kapn = kap0 + n
  nun = nu0 + n
  
  muna = kap0/kapn * nu0 + n/kapn * yba
  munb = kap0/kapn * nu0 + n/kapn * ybb
  
  sig2na = (1/nun) * ((nu0 * sig20) + (n-1)*s2a + (kap0 * n)/kapn * (yba - mu0)^2)
  sig2nb = (1/nun) * ((nu0 * sig20) + (n-1)*s2b + (kap0 * n)/kapn * (ybb - mu0)^2)
  
  # Monte Carlo Samples
  s = 10000
  
  sig2a_MC = 1 / rgamma(s, nun/2, nun/2 * sig2na)
  sig2b_MC = 1 / rgamma(s, nun/2, nun/2 * sig2nb)
  
  thetaa_MC = c()
  thetab_MC = c()
  for(i in seq(1,s)){
    thetaa_MC = c(thetaa_MC, rnorm(1, muna, sqrt(sig2a_MC[i]/kapn)))
    thetab_MC = c(thetab_MC, rnorm(1, munb, sqrt(sig2b_MC[i]/kapn)))
  }
  
  # Probability
  probs = c(probs, mean(thetaa_MC < thetab_MC))
}
plot(sizes, probs, main = "Probability vs. Prior Weights")

# Question 4
# (b)
s = 10000
y = c(2.16, 0.74, 1.87, 3.03, 3.11, 2.74, 1.23, 3.64, 1.57, 2.12)

ptl = function(z){
  result = zs * exp(-0.5 * (sum((y-zs)^2) + zs))
  return(result)
}

zs = rgamma(s, 0.5, 2)
qzs = dgamma(zs, 0.5, 2)
wtls = c()
for(i in seq(1,s)){
  wtls = c(wtls, (zs[i] * exp(-0.5 * (sum((y-zs[i])^2) + zs[i])))/qzs[i])
}
ws = wtls/(sum(wtls))

psih = sum(ws * zs)
psih

# (c)
ybar = mean(y)
C = (1/sqrt(pi))*exp(-0.5*sum((y-ybar)^2))

zs = c()
while(length(zs) < s){
  u = runif(1)
  zst = rgamma(1, 0.5, 2)
  b = (zst * exp(-0.5*(sum((y-zst)^2)) + zst))/(C*dgamma(zst, 0.5, 2))
  
  if(u < b){
    zs = c(zs, zst)
  }
}

mean(zs)


