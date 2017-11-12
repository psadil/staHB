functions {
  real lb(real mu_i, real mu_j, real top) {
    return mu_i - mu_j < 0 ? negative_infinity() : top;

  }
  real ub(real mu_i, real mu_j, real bot) {
    return mu_i - mu_j < 0 ? bot: positive_infinity();

  }
}
data {
  int n;
  matrix[n, 2] y;
  int condition[n];
  int n_condition;
}
parameters{
  cholesky_factor_corr[2] Omega;
  ordered[n_condition] mu_x;
  real mu1_y;
  real<lower = lb(mu_x[1], mu_x[2], mu1_y), upper = ub(mu_x[1], mu_x[2], mu1_y)> mu2_y;
}
model {
  vector[n_condition] mu_y = [mu1_y, mu2_y]';
  matrix[n_condition, 2] mu = append_col(mu_x, mu_y);

  Omega ~ lkj_corr_cholesky(1.5);

  mu_x ~ normal(0,1);
  mu1_y ~ normal(0,1);
  mu2_y ~ normal(0, 1) T[lb(mu_x[1], mu_x[2], mu1_y), ub(mu_x[1], mu_x[2], mu1_y)];

  for(i in 1:n){
    y[i] ~ multi_normal_cholesky(mu[condition[i]], Omega);
  }

}
