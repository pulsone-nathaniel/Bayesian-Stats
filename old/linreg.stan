data {
  int<lower=1> N; // number of observations
  vector[N] x; // predictor (sepal length)
  vector[N] y; // response (petal length)
}

parameters {
  real alpha; // intercept
  real beta; // slope
  real<lower=0> sigma; // residual SD
}

model {
  // weakly informative priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(0, 1); // half-N(0,1) due to lower=0
  // sampling model (vectorized)
  y ~ normal(alpha + beta * x, sigma);
}

generated quantities {
  // posterior predictive for each observed x (useful for checks/plots)
  vector[N] y_rep;
  for (n in 1:N) y_rep[n] = normal_rng(alpha + beta * x[n], sigma);
}
  