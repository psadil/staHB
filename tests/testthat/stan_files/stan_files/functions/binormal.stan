row_vector binormal_cdf_vector(matrix z, vector rho) {
  // NOTE: this funciton has a hard time computing values very close to 0.
  // In those cases, it sometimes incorrectly outputs a negative number.
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

row_vector biprobit_lpdf_real(matrix Y, matrix mu, real rho) {
  int n = rows(Y);
  matrix[n, 2] q = 2.0 * Y - 1.0;
  matrix[n, 2] z = q .* mu;
  vector[n] rho1 = col(q,1) .* col(q,2) * rho;

  return log(binormal_cdf_vector(z, rho1));
}

row_vector biprobit_lpdf_vector(matrix Y, matrix mu, vector rho) {
  int n = rows(Y);
  matrix[n, 2] q = 2.0 * Y - 1.0;
  matrix[n, 2] z = q .* mu;
  vector[n] rho1 = col(q,1) .* col(q,2) .* rho;

  return log(binormal_cdf_vector(z, rho1));
}
