data {
  int<lower=0> N;
  int<lower=1> J; //dimension of observation
  vector[J] X[N]; //observations
  vector[N] y;
  
}

parameters {
  corr_matrix[J] Omega;
  vector<lower=0>[J] sigma;
  vector[J] mu;
  vector[J] beta;
  real<lower=0> sig;
  real alpha;
}

transformed parameters {
  cov_matrix[J] Sigma;
  Sigma = quad_form_diag(Omega, sigma);
}

model {
  X ~ multi_normal(mu, Sigma);
  sigma ~ cauchy(0,5);
  Omega ~ lkj_corr(1);
  for(i in 1:N) {
    y[i] ~ normal(alpha + X[i]' * beta, sig);
  }
}

