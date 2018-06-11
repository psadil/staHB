functions{
#include /functions/binormal.stan
}
data {
  int<lower=0> n; // number of observations
  int<lower=1> n_condition; // number of conditions (4)
  int<lower=1> n_subject; // The number of subjects
  int<lower=1,upper=n_subject> subject[n]; // Index indicating the subject for a current trial
  int<lower=1> n_item; // The number of subjects
  int<lower=1,upper=n_item> item[n]; // Index indicating the subject for a current trial
  int<lower=1> D; // number of outcomes (2)
  int<lower=1,upper=n_condition> condition[n];
  matrix<lower = 0, upper=1>[n, D] y;
  vector[7] priors;
  int n_orders;
  matrix[n,n_condition] X[n_orders, D];
}
transformed data{
  row_vector[n_orders] zeroes_order = rep_row_vector(0,n_orders);
}
parameters{
  corr_matrix[D] condition_omega[n_condition];
  vector[n_orders-1] theta_raw; // first value must be pinned for identifiability
  matrix[n_condition-2, n_orders] zeta_raw[D]; // -1 for intercept, -1 for simplex identification
  row_vector[D] intercept;
  row_vector[D] condition_mu_raw; // basically forced to be positive by data
  vector<lower=0.0>[D] subject_scale; // Variance of subject-level effects
  matrix[D, n_subject] subject_mu_raw; // Subject-level coefficients for the bivariate normal means
  cholesky_factor_corr[D] subject_L;
  vector<lower=0.0>[D] item_scale; // Variance of subject-level effects
  matrix[D, n_item] item_mu_raw; // Subject-level coefficients for the bivariate normal means
  cholesky_factor_corr[D] item_L;
}
transformed parameters {
  vector[n_orders] theta_log = log_softmax(append_row(0,theta_raw));
  matrix[n_orders, n] lps;
  vector[n] log_lik;
  matrix[n_condition, D] condition_mu_ordered[n_orders]; // condition effects on mean of latent
  matrix[n_condition-1, n_orders] zeta[D]; // proportion of effect in each dimension
  vector[n_condition-1]  cumulative_sum_softmax_zeta[D, n_orders];
  matrix[n_subject, D] subject_mu = (diag_pre_multiply(subject_scale, subject_L) * subject_mu_raw)';
  matrix[n_item, D] item_mu = (diag_pre_multiply(item_scale, item_L) * item_mu_raw)';

  for (d in 1:D){
		 zeta[d] = append_row(zeroes_order, zeta_raw[d]);
	}

  // likelihood determined by mixture of possible orders
  // this would usually go in model block, but I want the log_lik for waic so
  // if put there I might need to calculate twice
  for(order in 1:n_orders){
		for(d in 1:D){
			cumulative_sum_softmax_zeta[d, order] = cumulative_sum( softmax( col( zeta[d], order)));
		}

    condition_mu_ordered[order] = append_row(intercept,
      append_col(condition_mu_raw[1] * cumulative_sum_softmax_zeta[1, order],
        condition_mu_raw[2] * cumulative_sum_softmax_zeta[2, order] ));

    lps[order,] = theta_log[order] + biprobit_lpdf_vector(y,
      subject_mu[subject] + item_mu[item] + append_col(X[order,1] * condition_mu_ordered[order,, 1], X[order,2] * condition_mu_ordered[order,, 2]),
      to_vector(condition_omega[condition,1,2]) );
  }
  for (i in 1:n){
    log_lik[i] = log_sum_exp(col(lps,i));
  }

}
model {

  // priors
  intercept ~ normal(0, priors[1]); // need an intercept to allow mix of positive and negative mo effects
  condition_mu_raw ~ normal(0, priors[2]);
	for(d in 1:D){
		to_vector(zeta_raw[d]) ~ normal(0, priors[3]);
	}

  subject_scale ~ gamma(priors[4], priors[5]);
  subject_L ~ lkj_corr_cholesky(priors[6]);
  to_vector(subject_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, subject_scale)
  item_scale ~ gamma(priors[4], priors[5]);
  item_L ~ lkj_corr_cholesky(priors[6]);
  to_vector(item_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, subject_scale)

  for(cond in 1:n_condition){
		condition_omega[cond] ~ lkj_corr(priors[6]);
	}

  theta_raw ~ normal(0, priors[7]);

  target += sum(log_lik);

}
generated quantities {
  vector[n_orders] theta = exp(theta_log);
  real<lower=-1.0, upper=1.0> subject_rho; // correlation to assemble
  real<lower=-1.0, upper=1.0> item_rho; // correlation to assemble
	vector<lower=-1, upper=1>[n_condition] condition_rho;

	for(cond in 1:n_condition){
		condition_rho[cond] = condition_omega[cond,1,2]; // correlation to assemble
	}
  {
    matrix[D, D] Sigma;
    Sigma = multiply_lower_tri_self_transpose(subject_L);
    subject_rho = Sigma[1,2];
    Sigma = multiply_lower_tri_self_transpose(item_L);
    item_rho = Sigma[1,2];
  }

}
