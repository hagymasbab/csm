data {
  int<lower=0> N;
	
  int<lower=0> d_x;
	int<lower=0> d_u;

  real<lower=0> z_shape;
  real<lower=0> z_scale;

  cov_matrix[d_u] C;
  matrix[d_u, d_x] A;
  cov_matrix[d_x] sigma_x;

  vector[d_x] x[N];
}

transformed data {
  real z_rate; 
  vector[d_u] mu_z;   

  z_rate <- 1 / z_scale;

  for (i in 1:d_u)
    mu_z[i] <- 0.0;
}

parameters {
  real<lower=0.001> z[N];
  vector[d_u] u[N];
}

model {
  z ~ gamma(z_shape,z_rate);
  u ~ multi_normal(mu_z,C);
  for (i in 1:N)
    x[i] ~ multi_normal(z[i]*A*u[i],sigma_x);
}