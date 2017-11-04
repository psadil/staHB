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
  int N;
  int P;
  int P2; // the number of instruments
  matrix[N,2] Y; // outcomes, binary
  matrix[N, P] X; // first column is a column of ones
  matrix[N, P2] Z; // matrix of instruments
}
transformed data {
  matrix[N, P+P2] X1;
  matrix[N, P + 1] X2;
  X1 = append_col(X, Z);
  X2 = append_col(Y[:,1], X);
}
parameters {
  vector[P + P2] beta;
  vector[P + 1] gamma;
  real<lower = -1, upper = 1> rho;
}
model {
  // priors
  beta ~ normal(0, 1);
  gamma ~ normal(0, 1);

  // likelihood
  for(i in 1:N) {
    Y[i] ~ biprobit(X1[i] * beta, X2[i] * gamma, rho);
  }
}
