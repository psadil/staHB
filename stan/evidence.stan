// generated with brms 1.10.2
functions { 
} 
data { 
  int<lower=1> N;  // total number of observations 

  int<lower=1> nresp;  // number of responses   
  int nrescor;  // number of residual correlations 
  vector[nresp] Y[N];  // response matrix
  int<lower=1> K_evidencex;  // number of population-level effects 
  matrix[N, K_evidencex] X_evidencex;  // population-level design matrix 
  int<lower=1> K_evidencey;  // number of population-level effects 
  matrix[N, K_evidencey] X_evidencey;  // population-level design matrix 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
  int Kc_evidencex = K_evidencex - 1; 
  matrix[N, K_evidencex - 1] Xc_evidencex;  // centered version of X_evidencex 
  vector[K_evidencex - 1] means_X_evidencex;  // column means of X_evidencex before centering 
  int Kc_evidencey = K_evidencey - 1; 
  matrix[N, K_evidencey - 1] Xc_evidencey;  // centered version of X_evidencey 
  vector[K_evidencey - 1] means_X_evidencey;  // column means of X_evidencey before centering 
  for (i in 2:K_evidencex) { 
    means_X_evidencex[i - 1] = mean(X_evidencex[, i]); 
    Xc_evidencex[, i - 1] = X_evidencex[, i] - means_X_evidencex[i - 1]; 
  } 
  for (i in 2:K_evidencey) { 
    means_X_evidencey[i - 1] = mean(X_evidencey[, i]); 
    Xc_evidencey[, i - 1] = X_evidencey[, i] - means_X_evidencey[i - 1]; 
  } 
} 
parameters { 
  vector[Kc_evidencex] b_evidencex;  // population-level effects 
  real temp_evidencex_Intercept;  // temporary intercept 
  vector[Kc_evidencey] b_evidencey;  // population-level effects 
  real temp_evidencey_Intercept;  // temporary intercept 
  // parameters for multivariate linear models 
  vector<lower=0>[nresp] sigma; 
  cholesky_factor_corr[nresp] Lrescor; 
} 
transformed parameters { 
  // cholesky factor of residual covariance matrix 
  cholesky_factor_cov[nresp] LSigma = diag_pre_multiply(sigma, Lrescor); 
} 
model { 
  vector[N] mu_evidencex = Xc_evidencex * b_evidencex + temp_evidencex_Intercept; 
  vector[N] mu_evidencey = Xc_evidencey * b_evidencey + temp_evidencey_Intercept; 
  // multivariate linear predictor matrix 
  vector[nresp] Mu[N]; 
  for (n in 1:N) { 
    Mu[n, 1] = mu_evidencex[n]; 
    Mu[n, 2] = mu_evidencey[n]; 
  } 
  // priors including all constants 
  target += student_t_lpdf(sigma | 3, 0, 10)
    - 2 * student_t_lccdf(0 | 3, 0, 10); 
  target += lkj_corr_cholesky_lpdf(Lrescor | 1); 
  // likelihood including all constants 
  if (!prior_only) { 
    target += multi_normal_cholesky_lpdf(Y | Mu, LSigma); 
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_evidencex_Intercept = temp_evidencex_Intercept - dot_product(means_X_evidencex, b_evidencex); 
  // actual population-level intercept 
  real b_evidencey_Intercept = temp_evidencey_Intercept - dot_product(means_X_evidencey, b_evidencey); 
  matrix[nresp, nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor); 
  vector<lower=-1,upper=1>[nrescor] rescor; 
  // take only relevant parts of residual correlation matrix 
  rescor[1] = Rescor[1, 2]; 
} 