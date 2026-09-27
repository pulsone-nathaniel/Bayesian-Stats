data {
int<lower=0> N;
int<lower=0> n[N];
int<lower=0> y[N];
vector[N] x;
}

parameters {
real alpha;
real beta;
}

model {
y ~ binomial_logit(n, alpha + beta * x);
}