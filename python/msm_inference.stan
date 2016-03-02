data {
  int<lower=1> N;
  int<lower=1> k;
    
  int<lower=1> d_x;
  int<lower=1> d_u;

  real<lower=0> z_shape;
  real<lower=0> z_scale;

  real<lower=0> g_shape;
  real<lower=0> g_scale;

  matrix[d_u,k] C;
  cov_matrix[d_u] Var_u;
  matrix[d_x, d_u] A;
  cov_matrix[d_x] sigma_x;

  vector[d_x] x[N];
}

transformed data {
  real z_rate; 
  real g_rate; 

  z_rate <- 1 / z_scale;
  g_rate <- 1 / g_scale;
}

parameters {
  vector[k] g[N];
  real<lower=0.001> z[N];
  vector[d_u] u[N];
}

model {
  vector[d_u] mu_u;
    
  z ~ gamma(z_shape,z_rate);  
  
  for (i in 1:N) {
    for (j in 1:k){
      g[i,j] ~ gamma(g_shape,g_rate);
    }
      
    mu_u <- C * g[i,:];

    u[i] ~ multi_normal(mu_u,Var_u);
    x[i] ~ multi_normal(z[i]*A*u[i],sigma_x);
  }
}