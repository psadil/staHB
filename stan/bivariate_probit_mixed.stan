functions {
  row_vector binormal_cdf_vector(matrix z, vector rho) {
    // NOTE: this funciton has a hard time computing values very close to 0.
    // In those cases, if often incorrectly outputs a negative number.
    // expressions come from https://github.com/stan-dev/stan/issues/2356
    // and https://papers.ssrn.com/sol3/papers.cfm?abstract_id=2924071

    int n = num_elements(rho);
    vector[n] denom = sqrt((1.0 + rho) .* (1.0 - rho));
    vector[n] a1 = (col(z,2) ./ col(z,1) - rho) ./ denom;
    vector[n] a2 = (col(z,1) ./ col(z,2) - rho) ./ denom;
    vector[n] product = col(z,1) .* col(z,2);
    matrix[n,2] z_probit = Phi(z); // too unstable to use Phi_approx safely
    real delta;
    row_vector[n] out;

    for(i in 1:n){
      if( rho[i] == 1){
        out[i] = min( row(z_probit,i) );
      }else if(rho[i] == -1){
        out[i] = sum( row(z_probit,i) ) - 1;
      }else if(rho[i] == 0){
        out[i] = prod( row(z_probit,i) );
      }else if(z[i,1] == 0 && z[i,2] == 0){
        out[i] = 0.25 + asin(rho[i]) / (2.0 * pi());
      }else {
        delta = product[i] < 0 || (product[i] == 0 && (z[i,1] + z[i,2]) < 0);
        out[i] = 0.5 * (z_probit[i,1] + z_probit[i,2] - delta) - owens_t(z[i,1], a1[i]) - owens_t(z[i,2], a2[i]);
      }
    }
    return out;
}
  row_vector biprobit_lpdf_vector(matrix Y, matrix mu, vector rho) {
    int n = num_elements(rho);
    matrix[n, 2] q = 2.0 * Y - 1.0;
    matrix[n, 2] z = q .* mu;
    vector[n] rho1 = col(q,1) .* col(q,2) .* rho;

    return log(binormal_cdf_vector(z, rho1));
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
  int<lower=1,upper=n_condition> condition[n];
  matrix<lower = 0, upper=1>[n, D] y;
  vector[7] priors;
  int n_orders;
  int order_X[D, n_orders, n_condition]; // conditions permutated by appropriate order
}
parameters{
  ordered[n_condition] condition_mu_ordered[D]; // condition effects on mean of latent
  matrix[D, n_subject] subject_mu_raw; // Subject-level coefficients for the bivariate normal means
  matrix[D, n_item] item_mu_raw;
  vector<lower=0.0>[D] subject_scale; // Variance of subject-level effects
  vector<lower=0.0>[D] item_scale;
  corr_matrix[D] condition_omega[n_condition];
  cholesky_factor_corr[D] subject_L;
  cholesky_factor_corr[D] item_L;
  vector[n_orders] theta_raw;
}
transformed parameters {
  matrix[D, n_subject] subject_mu; // Subject-level coefficients for the bivariate nrmal means
  matrix[D, n_item] item_mu; // Subject-level coefficients for the bivariate nrmal means
  matrix[n, D] Mu_ordered[n_orders];
  vector[n_orders] theta_log;

  subject_mu = diag_pre_multiply(subject_scale, subject_L) * subject_mu_raw;
  item_mu = diag_pre_multiply(item_scale, item_L) * item_mu_raw;

  {
    matrix[D, n] Mu_random = subject_mu[,subject] + item_mu[,item];

    for(order in 1:n_orders){
      Mu_ordered[order] = Mu_random' + append_col(condition_mu_ordered[1,order_X[1,order,condition]],
        condition_mu_ordered[2,order_X[2,order,condition]]);
    }
  }

  theta_log = log_softmax(theta_raw);

}
model {
  matrix[n_orders, n] lps;

  // priors
  for(d in 1:D){
    condition_mu_ordered[d] ~ normal(0, priors[1]);
  }
  for(k in 1:n_condition){
    condition_omega[k] ~ lkj_corr(priors[6]);
  }
  subject_scale ~ gamma(priors[4], priors[5]);
  item_scale ~ gamma(priors[4], priors[5]);
  subject_L ~ lkj_corr_cholesky(priors[6]);
  item_L ~ lkj_corr_cholesky(priors[6]);

  to_vector(subject_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, subject_scale)
  to_vector(item_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, item_scale)

  theta_raw ~ normal(0, priors[7]);

  // likelihood determined by mixture of possible orders
  for(order in 1:n_orders){
    lps[order,] = theta_log[order] + biprobit_lpdf_vector(y, Mu_ordered[order,,1:2], to_vector(condition_omega[condition,1,2]));
  }
  for (i in 1:n){
    target += log_sum_exp(col(lps,i));
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
    matrix[n_orders,n] lps;
    for(order in 1:n_orders){
      lps[order,] = theta_log[order] + biprobit_lpdf_vector(y, Mu_ordered[order,,1:2], to_vector(condition_omega[condition,1,2]));
    }
    for (i in 1:n){
      log_lik[i] = log_sum_exp(col(lps,i));
    }
  }
}
