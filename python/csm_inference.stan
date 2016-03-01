data {
  int<lower=1> N;
  int<lower=1> k;
    
  int<lower=1> d_x;
  int<lower=1> d_u;

  real<lower=0> z_shape;
  real<lower=0> z_scale;

  real<lower=0> g_shape;
  real<lower=0> g_scale;

  cov_matrix[d_u] C[k];
  matrix[d_u, d_x] A;
  cov_matrix[d_x] sigma_x;

  vector[d_x] x[N];
}

transformed data {
  real z_rate; 
  real g_rate; 
  vector[d_u] mu_u;   

  z_rate <- 1 / z_scale;
  g_rate <- 1 / g_scale;

  for (i in 1:d_u)
    mu_u[i] <- 0.0;
}

parameters {
  vector[k] g[N];
  real<lower=0.001> z[N];
  vector[d_u] u[N];
}

model {
  matrix[d_u,d_u] C_u;
    
  z ~ gamma(z_shape,z_rate);  
  for (i in 1:N) {
    for (j in 1:k)
      g[i,j] ~ gamma(g_shape,g_rate);

    C_u <- g[i,1] * C[1];
    for (j in 2:k)
      C_u <- C_u + g[i,j] * C[j];
      
    u ~ multi_normal(mu_u,C_u);
    x[i] ~ multi_normal(z[i]*A*u[i],sigma_x);
  }
}