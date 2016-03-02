data {
  int<lower=1> N;
  int<lower=1> k;
    
  int<lower=1> d_x;
  int<lower=1> d_u;

  real<lower=0> z_shape;
  real<lower=0> z_scale;

  real<lower=0> g_shape;
  real<lower=0> g_scale;

  matrix[d_u,d_u] C[k];
  cov_matrix[d_u] Var_u;
  matrix[d_x, d_u] A;
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
  vector<lower=0, upper=1>[k] g[N];
  real<lower=0.001> z[N];
  vector[d_u] u[N];
}

model {
  matrix[d_u,d_u] C_u;
  vector[2] actEl;
  actEl[1] <-  1.0;
    
  z ~ gamma(z_shape,z_rate);  
  
  for (i in 1:N) {
    for (j in 1:k)
      # g[i,j] ~ gamma(g_shape,g_rate);
      g[i,j] ~ beta(g_shape,g_scale);

    # add weighed components
    C_u <- g[i,1] * C[1];
    for (j in 2:k)
      C_u <- C_u + g[i,j] * C[j];

    # add unit diagonal and max out at 1
    for (l in 1:d_u){
      for (m in 1:d_u){
        actEl[2] <- C_u[l,m];
        C_u[l,m] <- if_else (l==m, 1.0, min(actEl));
      }
    }

    # transform with variance
    C_u <- Var_u * C_u * Var_u;
      
    u ~ multi_normal(mu_u,C_u);
    x[i] ~ multi_normal(z[i]*A*u[i],sigma_x);
  }
}