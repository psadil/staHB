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
  row_vector biprobit_lpdf_vector(matrix Y, matrix mu, real rho) {
    int n = rows(Y);
    matrix[n, 2] q = 2.0 * Y - 1.0;
    matrix[n, 2] z = q .* mu;
    vector[n] rho1 = col(q,1) .* col(q,2) * rho;

    return log(binormal_cdf_vector(z, rho1));
  }
}
data {
  int<lower=0> n; // number of observations
  int<lower=1> n_condition; // number of conditions (4)
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
  corr_matrix[D] condition_omega;
  vector[n_orders-1] theta_raw; // first value must be pinned for identifiability
  matrix[n_condition-2, n_orders] zeta_raw; // -1 for intercept, -1 for simplex identification
  row_vector[D] intercept;
  row_vector[D] condition_mu_raw; // basically forced to be positive by data
}
transformed parameters {
  vector[n_orders] theta_log = log_softmax(append_row(0,theta_raw));
  matrix[n_orders, n] lps;
  vector[n] log_lik;
  matrix[n_condition, D] condition_mu_ordered[n_orders]; // condition effects on mean of latent
  matrix[n_condition-1, n_orders] zeta = append_row(zeroes_order, zeta_raw);
  vector[n_condition-1]  cumulative_sum_softmax_zeta[n_orders];

  // likelihood determined by mixture of possible orders
  // this would usually go in model block, but I want the log_lik for waic so
  // if put there I might need to calculate twice
  for(order in 1:n_orders){
    cumulative_sum_softmax_zeta[order] = cumulative_sum( softmax( col( zeta, order)));

    condition_mu_ordered[order] = append_row(intercept,
      append_col(condition_mu_raw[1] * cumulative_sum_softmax_zeta[order],
        condition_mu_raw[2] * cumulative_sum_softmax_zeta[order] ));

    lps[order,] = theta_log[order] + biprobit_lpdf_vector(y,
      append_col(X[order,1] * condition_mu_ordered[order,, 1], X[order,2] * condition_mu_ordered[order,, 2]),
      condition_omega[1,2] );
  }
  for (i in 1:n){
    log_lik[i] = log_sum_exp(col(lps,i));
  }

}
model {

  // priors
  intercept ~ normal(0, priors[1]); // need an intercept to allow mix of positive and negative mo effects
  condition_mu_raw ~ normal(0, priors[2]);
  to_vector(zeta_raw) ~ normal(0, priors[3]);

  condition_omega ~ lkj_corr(priors[6]);

  theta_raw ~ normal(0, priors[7]);

  target += sum(log_lik);

}
generated quantities {
  real<lower=-1, upper=1> condition_rho = condition_omega[1,2]; // correlation to assemble
  vector[n_orders] theta = exp(theta_log);
}
