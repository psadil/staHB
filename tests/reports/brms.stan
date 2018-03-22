// generated with brms 2.0.1
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
  vector[N] Y_y1;  // response variable 
  int<lower=1> Kmo_mu1_y1;  // number of monotonic effects
  int<lower=1> Imo_mu1_y1;  // number of monotonig variables
  int<lower=2> Jmo_mu1_y1[Imo_mu1_y1];  // length of simplexes
  // monotonic variables 
  int Xmo_mu1_y1_1[N];
  // prior concentration of monotonic simplexes
  vector[Jmo_mu1_y1[1]] con_simo_mu1_y1_1;
  int<lower=1> Kmo_mu2_y1;  // number of monotonic effects
  int<lower=1> Imo_mu2_y1;  // number of monotonig variables
  int<lower=2> Jmo_mu2_y1[Imo_mu2_y1];  // length of simplexes
  // monotonic variables 
  int Xmo_mu2_y1_1[N];
  // prior concentration of monotonic simplexes
  vector[Jmo_mu2_y1[1]] con_simo_mu2_y1_1;
  vector[2] con_theta_y1;  // prior concentration 
  vector[N] Y_y2;  // response variable 
  int<lower=1> Kmo_mu1_y2;  // number of monotonic effects
  int<lower=1> Imo_mu1_y2;  // number of monotonig variables
  int<lower=2> Jmo_mu1_y2[Imo_mu1_y2];  // length of simplexes
  // monotonic variables 
  int Xmo_mu1_y2_1[N];
  // prior concentration of monotonic simplexes
  vector[Jmo_mu1_y2[1]] con_simo_mu1_y2_1;
  int<lower=1> Kmo_mu2_y2;  // number of monotonic effects
  int<lower=1> Imo_mu2_y2;  // number of monotonig variables
  int<lower=2> Jmo_mu2_y2[Imo_mu2_y2];  // length of simplexes
  // monotonic variables 
  int Xmo_mu2_y2_1[N];
  // prior concentration of monotonic simplexes
  vector[Jmo_mu2_y2[1]] con_simo_mu2_y2_1;
  vector[2] con_theta_y2;  // prior concentration 
  int prior_only;  // should the likelihood be ignored? 
} 
transformed data { 
} 
parameters { 
  // scale of monotonic effects 
  vector[Kmo_mu1_y1] bmo_mu1_y1; 
  // simplexes of monotonic effects 
  simplex[Jmo_mu1_y1[1]] simo_mu1_y1_1; 
  real<lower=0> sigma1_y1;  // residual SD 
  // scale of monotonic effects 
  vector[Kmo_mu2_y1] bmo_mu2_y1; 
  // simplexes of monotonic effects 
  simplex[Jmo_mu2_y1[1]] simo_mu2_y1_1; 
  real<lower=0> sigma2_y1;  // residual SD 
  simplex[2] theta_y1;  // mixing proportions 
  ordered[2] ordered_Intercept_y1;  // to identify mixtures 
  // scale of monotonic effects 
  vector[Kmo_mu1_y2] bmo_mu1_y2; 
  // simplexes of monotonic effects 
  simplex[Jmo_mu1_y2[1]] simo_mu1_y2_1; 
  real<lower=0> sigma1_y2;  // residual SD 
  // scale of monotonic effects 
  vector[Kmo_mu2_y2] bmo_mu2_y2; 
  // simplexes of monotonic effects 
  simplex[Jmo_mu2_y2[1]] simo_mu2_y2_1; 
  real<lower=0> sigma2_y2;  // residual SD 
  simplex[2] theta_y2;  // mixing proportions 
  ordered[2] ordered_Intercept_y2;  // to identify mixtures 
} 
transformed parameters { 
  // identify mixtures via ordering of the intercepts 
  real temp_mu1_y1_Intercept = ordered_Intercept_y1[1]; 
  // identify mixtures via ordering of the intercepts 
  real temp_mu2_y1_Intercept = ordered_Intercept_y1[2]; 
  // mixing proportions 
  real<lower=0,upper=1> theta1_y1 = theta_y1[1]; 
  real<lower=0,upper=1> theta2_y1 = theta_y1[2]; 
  // identify mixtures via ordering of the intercepts 
  real temp_mu1_y2_Intercept = ordered_Intercept_y2[1]; 
  // identify mixtures via ordering of the intercepts 
  real temp_mu2_y2_Intercept = ordered_Intercept_y2[2]; 
  // mixing proportions 
  real<lower=0,upper=1> theta1_y2 = theta_y2[1]; 
  real<lower=0,upper=1> theta2_y2 = theta_y2[2]; 
} 
model { 
  vector[N] mu1_y1 = rep_vector(0, N) + temp_mu1_y1_Intercept; 
  vector[N] mu2_y1 = rep_vector(0, N) + temp_mu2_y1_Intercept; 
  vector[N] mu1_y2 = rep_vector(0, N) + temp_mu1_y2_Intercept; 
  vector[N] mu2_y2 = rep_vector(0, N) + temp_mu2_y2_Intercept; 
  for (n in 1:N) { 
    mu1_y1[n] = mu1_y1[n] + (bmo_mu1_y1[1]) * mo(simo_mu1_y1_1, Xmo_mu1_y1_1[n]); 
    mu2_y1[n] = mu2_y1[n] + (bmo_mu2_y1[1]) * mo(simo_mu2_y1_1, Xmo_mu2_y1_1[n]); 
    mu1_y2[n] = mu1_y2[n] + (bmo_mu1_y2[1]) * mo(simo_mu1_y2_1, Xmo_mu1_y2_1[n]); 
    mu2_y2[n] = mu2_y2[n] + (bmo_mu2_y2[1]) * mo(simo_mu2_y2_1, Xmo_mu2_y2_1[n]); 
  } 
  // priors including all constants 
  target += student_t_lpdf(temp_mu1_y1_Intercept | 3, 1, 10); 
  target += dirichlet_lpdf(simo_mu1_y1_1 | con_simo_mu1_y1_1); 
  target += student_t_lpdf(sigma1_y1 | 3, 0, 10); 
  target += student_t_lpdf(temp_mu2_y1_Intercept | 3, 1, 10); 
  target += dirichlet_lpdf(simo_mu2_y1_1 | con_simo_mu2_y1_1); 
  target += student_t_lpdf(sigma2_y1 | 3, 0, 10); 
  target += dirichlet_lpdf(theta_y1 | con_theta_y1); 
  target += student_t_lpdf(temp_mu1_y2_Intercept | 3, 0, 10); 
  target += dirichlet_lpdf(simo_mu1_y2_1 | con_simo_mu1_y2_1); 
  target += student_t_lpdf(sigma1_y2 | 3, 0, 10); 
  target += student_t_lpdf(temp_mu2_y2_Intercept | 3, 0, 10); 
  target += dirichlet_lpdf(simo_mu2_y2_1 | con_simo_mu2_y2_1); 
  target += student_t_lpdf(sigma2_y2 | 3, 0, 10); 
  target += dirichlet_lpdf(theta_y2 | con_theta_y2); 
  // likelihood including all constants 
  if (!prior_only) { 
    for (n in 1:N) { 
      real ps[2]; 
      ps[1] = log(theta1_y1) + normal_lpdf(Y_y1[n] | mu1_y1[n], sigma1_y1); 
      ps[2] = log(theta2_y1) + normal_lpdf(Y_y1[n] | mu2_y1[n], sigma2_y1); 
      target += log_sum_exp(ps); 
    } 
    for (n in 1:N) { 
      real ps[2]; 
      ps[1] = log(theta1_y2) + normal_lpdf(Y_y2[n] | mu1_y2[n], sigma1_y2); 
      ps[2] = log(theta2_y2) + normal_lpdf(Y_y2[n] | mu2_y2[n], sigma2_y2); 
      target += log_sum_exp(ps); 
    } 
  } 
} 
generated quantities { 
  // actual population-level intercept 
  real b_mu1_y1_Intercept = temp_mu1_y1_Intercept; 
  // actual population-level intercept 
  real b_mu2_y1_Intercept = temp_mu2_y1_Intercept; 
  // actual population-level intercept 
  real b_mu1_y2_Intercept = temp_mu1_y2_Intercept; 
  // actual population-level intercept 
  real b_mu2_y2_Intercept = temp_mu2_y2_Intercept; 
} 