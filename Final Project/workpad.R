# Final Project Raw code

# init libraries
library(rstan)
library(coda)
library(bayesplot)

# init data
df = read.csv("Final Project/diabetes.csv")

Y = df$Outcome

pregnancies = df$Pregnancies
glucose = df$Glucose
age = df$Age

X = cbind(1, pregnancies, glucose, age)

rm(pregnancies, glucose, age)


# define stan model
stan_data <- list(
  N = length(Y),
  K = ncol(X),
  Y = Y,
  X = X
)

fit <- stan(
  file = "Final Project/log_reg.stan",
  data = stan_data,
  iter = 2000,
  warmup = 1000,
  chains = 1,
  seed = 123
)

mcmc_trace(fit)
mcmc_dens(as.array(fit), pars = "beta[1]")

p <- mcmc_trace(fit, pars = c("beta[1]", "beta[2]", "beta[3]"))
print(p)


par(mfrow = c(2,2))
plot(n_prior, b_hats[1,], ylab = "Inctercept", type = "l")
plot(n_prior, b_hats[2,], ylab = "Glucose", type = "l")
plot(n_prior, b_hats[3,], ylab = "Skin Thickness", type = "l")
plot(n_prior, b_hats[4,], ylab = "Insulin", type = "l")

par(mfrow = c(1,4))
hist(y, main = "Histogram of Diabetes Pedigree")
hist(gluc, main = "Histogram of Glucose Levels")
hist(bmi, main = "Histogram of Body Mass Index")
hist(bp, main = "Histogram of Blood Pressure")

par(mfrow = c(1, 2))
mcmc_acf(sim)
mcmc_trace(sim)
mcmc_dens(sim)
pp_check(fit)

test_ind = sample(1:length(y), 100)

y_test = y[test_ind]
X_test = X[test_ind,]

beta_bma = c(mean(beta[,1]),
              mean(beta[,2]),
              mean(beta[,3]),
              mean(beta[,4]),
              mean(beta[,5]),
              mean(beta[,6]),
              mean(beta[,7]),
              mean(beta[,8]))

beta_map = as.numeric(b_hats[,2])

sum((y_test - X_test%*%beta_bma)^2) / 100
sum((y_test - X_test[,c(1, 3, 5, 6)]%*%beta_map)^2) / 100
