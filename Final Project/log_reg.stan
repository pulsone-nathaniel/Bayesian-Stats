data {
  int<lower=1> N;               // number of observations
  int<lower=1> K;               // number of predictors
  int<lower=0,upper=1> Y[N];    // binary outcome
  matrix[N, K] X;               // design matrix
}

parameters {
  vector[K] beta;               // regression coefficients
}

model {
  // Priors
  beta ~ normal(0, 5);          // weakly-informative prior

  // Likelihood
  Y ~ bernoulli_logit(X * beta);
}