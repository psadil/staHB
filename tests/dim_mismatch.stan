data {
  int a;
  int b;
  int c;
  int d;
  vector[c] y;
  matrix[c,d] X[b,a];
}
parameters{
  vector[d] beta;
}
model {
  beta ~ normal(0, 1);
  y ~ normal(X[1,1]*beta, 1);
}
