// generated with brms 1.10.2
functions { 

  /* compute monotonic effects
   * Args:
   *   scale: a simplex parameter
   *   i: index to sum over the simplex
   * Returns:
   *   a scalar between 0 and 1
   */
  real mo(vector scale, int i) {
    if (i == 0) {
      return 0;
    } else {
      return(sum(scale[1:i]));
    }
  }
} 
data { 
  int<lower=1> N;  // total number of observations 

  int<lower=1> nresp;  // number of responses   
  int nrescor;  // number of residual correlations 
  vector[nresp] Y[N];  // response matrix
  int<lower=1> Kmo_evidencex;  // number of monotonic effects
  int<lower=1> Imo_evidencex;  // number of monotonig variables
  int<lower=2> Jmo_evidencex[Imo_evidencex];  // length of simplexes
  // monotonic variables 
  int Xmo_evidencex_1[N];
  // prior concentration of monotonic simplexes
  vector[Jmo_evidencex[1]] con_simo_evidencex_1;
  int<lower=1> Kmo_evidencey;  // number of monotonic effects
  int<lower=1> Imo_evidencey;  // number of monotonig variables
  int<lower=2> Jmo_evidencey[Imo_evidencey];  // length of simplexes
  // monotonic variables 
  int Xmo_evidencey_1[N];
  // prior concentration of monotonic simplexes
  vector[Jmo_evidencey[1]] con_simo_evidencey_1;
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
} 
parameters { 
  real temp_evidencex_Intercept;  // temporary intercept 
  // scale of monotonic effects 
  vector[Kmo_evidencex] bmo_evidencex; 
  // simplexes of monotonic effects 
  simplex[Jmo_evidencex[1]] simo_evidencex_1; 
  real temp_evidencey_Intercept;  // temporary intercept 
  // scale of monotonic effects 
  vector[Kmo_evidencey] bmo_evidencey; 
  // simplexes of monotonic effects 
  simplex[Jmo_evidencey[1]] simo_evidencey_1; 
  // parameters for multivariate linear models 
  vector<lower=0>[nresp] sigma; 
  cholesky_factor_corr[nresp] Lrescor; 
} 
transformed parameters { 
  // cholesky factor of residual covariance matrix 
  cholesky_factor_cov[nresp] LSigma = diag_pre_multiply(sigma, Lrescor); 
} 
model { 
  vector[N] mu_evidencex = rep_vector(0, N) + temp_evidencex_Intercept; 
  vector[N] mu_evidencey = rep_vector(0, N) + temp_evidencey_Intercept; 
  // multivariate linear predictor matrix 
  vector[nresp] Mu[N]; 
  for (n in 1:N) { 
    mu_evidencex[n] = mu_evidencex[n] + (bmo_evidencex[1]) * mo(simo_evidencex_1, Xmo_evidencex_1[n]); 
    mu_evidencey[n] = mu_evidencey[n] + (bmo_evidencey[1]) * mo(simo_evidencey_1, Xmo_evidencey_1[n]); 
    Mu[n, 1] = mu_evidencex[n]; 
    Mu[n, 2] = mu_evidencey[n]; 
  } 
  // priors including all constants 
  target += dirichlet_lpdf(simo_evidencex_1 | con_simo_evidencex_1); 
  target += dirichlet_lpdf(simo_evidencey_1 | con_simo_evidencey_1); 
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
  real b_evidencex_Intercept = temp_evidencex_Intercept; 
  // actual population-level intercept 
  real b_evidencey_Intercept = temp_evidencey_Intercept; 
  matrix[nresp, nresp] Rescor = multiply_lower_tri_self_transpose(Lrescor); 
  vector<lower=-1,upper=1>[nrescor] rescor; 
  // take only relevant parts of residual correlation matrix 
  rescor[1] = Rescor[1, 2]; 
} 
