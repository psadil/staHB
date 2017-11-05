data {
  int<lower=1> n_subject; // The number of subjects
  int<lower=1> n_item; // The number of items
  int<lower=0> n; // number of observations
  int<lower=1> n_condition; // number of conditions (4)
  int<lower=1> D; // number of outcomes (2)
  int<lower=1,upper=n_subject> subject[n]; // Index indicating the subject for a current trial
  int<lower=1,upper=n_item> item[n]; // Index indicating the subject for a current trial
  int<lower=1,upper=n_condition> condition[n]; // Index indicating which questions were gotten correct
  matrix<lower=0,upper=1>[n, n_condition] X; // design matrix
  int<lower = 0, upper=1> y[n, D];
  int<lower = 1> sum_y;
  vector[6] priors; // lkj + (3*mean of latent) + eta
  int idn[n_condition];
  int id1[idn[1]];
  int id2[idn[2]];
  int id3[idn[3]];
  int id4[idn[4]];
}
transformed data {
  int<lower=0> N_pos;
  int<lower=1,upper=n> n_pos[sum_y];
  int<lower=1,upper=D> d_pos[size(n_pos)];
  int<lower=0> N_neg;
  int<lower=1,upper=n> n_neg[(n * D) - size(n_pos)];
  int<lower=1,upper=D> d_neg[size(n_neg)];

  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i = 1;
    int j = 1;

    for (o in 1:n) {
      for (d in 1:D) {
        if (y[o,d] == 1) {
          n_pos[i] = o;
          d_pos[i] = d;
          i = i + 1;
        } else {
          n_neg[j] = o;
          d_neg[j] = d;
          j = j + 1;
        }
      }
    }
  }
}
parameters{
  matrix[D, n_condition] condition_mu; // condition effects on mean of latent
  matrix[D, n_subject] subject_mu_raw; // Subject-level coefficients for the bivariate normal means
  matrix[D, n_item] item_mu_raw; // Subject-level coefficients for the bivariate normal means
  vector<lower=0.0>[D] subject_scale; // Variance of subject-level effects
  vector<lower=0.0>[D] item_scale; // Variance of subject-level effects
  cholesky_factor_corr[D] condition_L[n_condition];
  cholesky_factor_corr[D] subject_L;
  cholesky_factor_corr[D] item_L;
  vector<lower=0>[n_condition] eta_lkj; // condition effects on mean of latent
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
}
transformed parameters {
  // We reconstruct the latent variable z for each subject
  // from its positive and negative components.
  vector[D] Mu[n]; // Bivariate nrmal means per trial
  vector[D] z[n];
  matrix[D, n_subject] subject_mu; // Subject-level coefficients for the bivariate nrmal means
  matrix[D, n_item] item_mu; // Subject-level coefficients for the bivariate nrmal means

  for (o in 1:N_pos)
    z[n_pos[o], d_pos[o] ] = z_pos[o];
  for (o in 1:N_neg)
    z[n_neg[o], d_neg[o] ] = z_neg[o];

  subject_mu = diag_pre_multiply(subject_scale, subject_L) * subject_mu_raw;
  item_mu = diag_pre_multiply(item_scale, item_L) * item_mu_raw;

  for (o in 1:n)
    Mu[o] = col(condition_mu, condition[o]) + col(subject_mu, subject[o]) + col(item_mu, item[o]);

}
model {

  // priors
  to_vector(condition_mu) ~ normal(0, priors[1]);
  eta_lkj ~ gamma(priors[2], priors[3]);
  for(k in 1:n_condition){
    condition_L[k] ~ lkj_corr_cholesky(eta_lkj[k]);
  }
  subject_scale ~ gamma(priors[4], priors[5]);
  item_scale ~ gamma(priors[4], priors[5]);
  subject_L ~ lkj_corr_cholesky(priors[6]);
  item_L ~ lkj_corr_cholesky(priors[6]);

  to_vector(subject_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, scale_s)
  to_vector(item_mu_raw) ~ normal(0, 1);  // implies ~ normal(0, scale_s)

  // likelihood
  for (k in 1:n_condition){
    if(k == 1){
      z[id1] ~ multi_normal_cholesky(Mu[id1], condition_L[k]);
    }else if(k == 2){
      z[id2] ~ multi_normal_cholesky(Mu[id2], condition_L[k]);
    }else if(k == 3){
      z[id3] ~ multi_normal_cholesky(Mu[id3], condition_L[k]);
    }else{
      z[id4] ~ multi_normal_cholesky(Mu[id4], condition_L[k]);
    }
  }

}
generated quantities {
  vector[n_condition] condition_rho; // correlation to assemble
  real subject_rho; // correlation to assemble
  real item_rho; // correlation to assemble

  {
    matrix[D, D] Sigma;
    for (k in 1:n_condition){
      Sigma = multiply_lower_tri_self_transpose(condition_L[k]);
      condition_rho[k] = Sigma[1,2];
    }
    Sigma = multiply_lower_tri_self_transpose(subject_L);
    subject_rho = Sigma[1,2];
    Sigma = multiply_lower_tri_self_transpose(item_L);
    item_rho = Sigma[1,2];
  }

}
