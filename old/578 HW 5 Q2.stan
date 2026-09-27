data {
  real<lower=0> alpha;
  real<lower=0> beta;
  int<lower=1> n;
  int<lower=0, upper=n> Y;
}
parameters {
  real<lower=0, upper=1> theta;
}
model {
  Y ~ binomial(n, theta);
  theta ~ beta(alpha, beta);
}