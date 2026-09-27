data {
  int<lower = 0> N;
  int<lower = 0> y[N];
}
parameters {
  real<lower = 0> lambda;
}
model {
  lambda ~ gamma(3, 1);
  y ~ poisson(lambda);
}
generated quantities {
  int<lower = 0> y_tilde = poisson_rng(lambda);
}