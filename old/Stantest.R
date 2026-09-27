data(iris)
head(iris)
iris_v <- subset(iris, Species == "versicolor")
x <- iris_v$Sepal.Length # predictor
y <- iris_v$Petal.Length # response
N <- length(y)
data_list <- list(N = N, x = as.vector(x), y = as.vector(y))

library(rstan)
options(mc.cores = parallel::detectCores())
fit <- stan(file = "linreg.stan",
            data = data_list,
            chains = 4, iter = 4000, warmup = 2000, seed = 2025,
            control = list(adapt_delta = 0.9, max_treedepth = 12))
print(fit, pars = c("alpha","beta","sigma"),
      probs = c(0.025, 0.5, 0.975))



library(rstan)
library(coda)
library(bayesplot)
glucose.dat = read.table("glucose.dat.txt")[,1]

# initialize priors
a = 1
b = 1
mu_0 = 120
tao2_0 = 200
sig2_0 = 1000
nu_0 = 10
k_0 = 10

# initialize algorithm parameters and samples
Ytild_s = c()
theta_mins = c()
theta_maxs = c()

# initialize sample data metrics
Y = glucose.dat
n = length(Y)

# Gibbs sampler
iters = 10000
for(i in seq(1,iters)){
  p = rbeta(1, a, b)
  X = rbinom(n, 1, p)
  
  # here, 0s are assigned to group 1, and 1s assigned to group 2
  
  n_1s = sum(1 - X)
  n_2s = n - n_1s
  
  a = a + n_2s
  b = b + n_1s
  
  Y1s = Y[which(X == 0)] # observations assigned to group 1
  Y2s = Y[which(X == 1)] # to group 2
  
  Ybar_1s = mean(Y1s)
  Ybar_2s = mean(Y2s)
  
  # update group 1 parameters
  
  sig2_n1s = (1 / (nu_0 + n_1s)) * (nu_0*sig2_0 + sum((Y1s - Ybar_1s)^2) + k_0*n_1s/(k_0 + n_1s) *(Ybar_1s - mu_0)^2)
  nu_1s = nu_0 + n_1s
  
  sig2_1s = 1 / (rgamma(1, nu_1s/2, nu_1s/2 * sig2_n1s)) # sample sig2_1
  
  lambda_0 = 1/tao2_0
  lambda_n = lambda_0 + n_1s / sig2_1s
  
  mu_1s = (lambda_0/lambda_n)*mu_0 + ((n_1s / sig2_1s)/lambda_n)*Ybar_1s
  tao_1s = 1/lambda_n
  
  theta_1s = rnorm(1, mu_1s, sqrt(tao_1s)) # sample theta_1
  
  
  # update group 2 parameters
  
  sig2_n2s = (1 / (nu_0 + n_2s)) * (nu_0*sig2_0 + sum((Y2s - Ybar_2s)^2) + k_0*n_2s/(k_0*n_2s) *(Ybar_2s - mu_0)^2)
  nu_2s = nu_0 + n_2s
  
  sig2_2s = 1 / (rgamma(1, nu_2s/2, nu_2s/2 * sig2_n2s)) # sample sig2_2
  
  lambda_n = lambda_0 + n_2s / sig2_2s
  
  mu_2s = (lambda_0/lambda_n)*mu_0 + ((n_2s / sig2_2s)/lambda_n)*Ybar_2s
  tao_2s = 1/lambda_n
  
  theta_2s = rnorm(1, mu_2s, sqrt(tao_2s)) # sample theta_2
  
  # update theta samples for part (c)
  theta_mins = c(theta_mins, min(c(theta_1s, theta_2s)))
  theta_maxs = c(theta_maxs, max(c(theta_1s, theta_2s)))
  
  # update Y samples for part (d)
  xs = rbinom(1, 1, p)
  if(xs == 0){
    Ytild_s = c(Ytild_s, rnorm(1, theta_1s, sqrt(sig2_1s)))
  }
  if(xs == 1){
    Ytild_s = c(Ytild_s, rnorm(1, theta_2s, sqrt(sig2_2s)))
  }
  
}
acf(theta_maxs)
acf(theta_mins)
effectiveSize(theta_maxs)
effectiveSize(theta_mins)

hist(Ytild_s, breaks = 30)
