functions {
    real binormal_cdf(real z1, real z2, real rho) {
    if (z1 != 0 || z2 != 0) {
      real denom = fabs(rho) < 1.0 ? sqrt((1 + rho) * (1 - rho)) : not_a_number();
      real a1 = (z2 / z1 - rho) / denom;
      real a2 = (z1 / z2 - rho) / denom;
      real product = z1 * z2;
      real delta = product < 0 || (product == 0 && (z1 + z2) < 0);
      return 0.5 * (Phi(z1) + Phi(z2) - delta) - owens_t(z1, a1) - owens_t(z2, a2);
    }
    return 0.25 + asin(rho) / (2 * pi());
  }
  real biprobit_lpdf(row_vector Y, real mu1, real mu2, real rho) {
    real q1;
    real q2;
    real w1;
    real w2;
    real rho1;

    // from Greene's econometrics book
    q1 = 2*Y[1] - 1.0;
    q2 = 2*Y[2] - 1.0;

    w1 = q1*mu1;
    w2 = q2*mu2;

    rho1 = q1*q2*rho;
    return log(binormal_cdf(w1, w2, rho1));
  }
}
data {
  int<lower=1> n_subject; // The number of subjects
  int<lower=1> n_item; // The number of items
  int<lower=0> n; // number of observations
  int<lower=1> n_condition; // number of conditions (4)
  int<lower=1> D; // number of outcomes (2)
  int<lower=1,upper=n_subject> subject[n]; // Index indicating the subject for a current trial
  int<lower=1,upper=n_item> item[n]; // Index indicating the subject for a current trial
  int<lower=1,upper=n_condition> condition[n]; // Index indicating which questions were gotten correct
  matrix<lower = 0, upper=1>[n, D] y;
  vector[7] priors;
  int n_orders;
  int order_X[D, n_orders, n_condition]; // conditions permutated by appropriate order
}
parameters{
  matrix[n_condition,D] condition_mu; // condition effects on mean of latent
  matrix[D, n_subject] subject_mu_raw; // Subject-level coefficients for the bivariate normal means
  matrix[D, n_item] item_mu_raw; // Subject-level coefficients for the bivariate normal means
  vector<lower=0.0>[D] subject_scale; // Variance of subject-level effects
  vector<lower=0.0>[D] item_scale; // Variance of subject-level effects
  corr_matrix[D] condition_omega[n_condition];
  cholesky_factor_corr[D] subject_L;
  cholesky_factor_corr[D] item_L;
  // vector<lower=0>[n_condition] eta_lkj; // condition effects on mean of latent
  vector[n_orders] theta_raw;
}
transformed parameters {
  matrix[D, n] Mu_unordered; // Bivariate nrmal means per trial
  matrix[D, n_subject] subject_mu; // Subject-level coefficients for the bivariate nrmal means
  matrix[D, n_item] item_mu; // Subject-level coefficients for the bivariate nrmal means
  matrix[n, D] Mu_ordered[n_orders];
  vector[n_orders] theta_log;
  matrix[n_condition, D] condition_mu_ordered;

  subject_mu = diag_pre_multiply(subject_scale, subject_L) * subject_mu_raw;
  item_mu = diag_pre_multiply(item_scale, item_L) * item_mu_raw;

  Mu_unordered = subject_mu[,subject] + item_mu[,item];

  for(order in 1:n_orders){
    for(d in 1:D){
      condition_mu_ordered[,d] = cumulative_sum(condition_mu[order_X[d, order,],d]);
    }
    Mu_ordered[order] = Mu_unordered' + condition_mu_ordered[condition];
  }

  theta_log = log_softmax(theta_raw);

}
model {
  vector[n_orders] lps; // temporary variable to store

  // priors
  to_vector(condition_mu) ~ normal(0, 1);
  // eta_lkj ~ gamma(priors[2], priors[3]);
  for(k in 1:n_condition){
    condition_omega[k] ~ lkj_corr(priors[6]);
  }
  subject_scale ~ gamma(priors[4], priors[5]);
  item_scale ~ gamma(priors[4], priors[5]);
  subject_L ~ lkj_corr_cholesky(priors[6]);
  item_L ~ lkj_corr_cholesky(priors[6]);

  to_vector(subject_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, scale_s)
  to_vector(item_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, scale_s)

  theta_raw ~ normal(0, priors[7]);

  // likelihood determined by mixture of possible orders
  for (i in 1:n){
    for(order in 1:n_orders){
      lps[order] = theta_log[order] + biprobit_lpdf(row(y,i) | Mu_ordered[order,i,1], Mu_ordered[order,i,2], condition_omega[condition[i],1,2]);
    }
    target += log_sum_exp(lps);
  }

}
generated quantities {
  vector[n_condition] condition_rho; // correlation to assemble
  real subject_rho; // correlation to assemble
  real item_rho; // correlation to assemble
  vector[n] log_lik;
  vector[n_orders] theta = exp(theta_log);

  {
    matrix[D, D] Sigma;
    for (k in 1:n_condition){
      condition_rho[k] = condition_omega[k,1,2];
    }
    Sigma = multiply_lower_tri_self_transpose(subject_L);
    subject_rho = Sigma[1,2];
    Sigma = multiply_lower_tri_self_transpose(item_L);
    item_rho = Sigma[1,2];
  }


  {
    vector[n_orders] lps; // temporary variable to store
    for (i in 1:n){
      for(order in 1:n_orders){
        lps[order] = theta_log[order] + biprobit_lpdf(row(y,i) | Mu_ordered[order,i,1], Mu_ordered[order,i,2], condition_omega[condition[i],1,2]);
      }
      log_lik[i] = log_sum_exp(lps);
    }
  }


}
