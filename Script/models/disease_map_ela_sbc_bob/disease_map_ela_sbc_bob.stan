functions {
  matrix K_functor (vector phi,
                    vector[] x,
                    real[] delta, int[] delta_int) {
    int n_obs = delta_int[1];
    matrix[n_obs, n_obs] K = cov_exp_quad(x, phi[1], phi[2]);
    for (i in 1:n_obs) K[i, i] += 1e-8;
    return K;
  }
}

data {
  int <lower = 0> n_obs;
  int <lower = 0> n_covariates;
  vector[n_covariates] x[n_obs];
  vector[n_obs] ye;
}

transformed data{
  real tol = 1e-6;
  real delta[0];
  int delta_int[1] = {n_obs};
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  int n_phi = 2;
  int n_samples[n_obs] = rep_array(1, n_obs);
  
  real <lower = 0> rho_alpha_prior = 2.42393;
  real <lower = 0> rho_beta_prior = 14.8171;
  real<lower = 0> rho_sim = inv_gamma_rng(rho_alpha_prior, rho_beta_prior);
  real<lower = 0> alpha_sim = inv_gamma_rng(10, 10);
  vector[n_obs] eta_sim = to_vector(normal_rng(rep_vector(0, n_obs), rep_vector(1, n_obs)));
  matrix[n_obs, n_obs] K_sim = cov_exp_quad(x, alpha_sim, rho_sim);
  for (n in 1:n_obs) K_sim[n, n] = K_sim[n,n] + 1e-8; // delta;
  matrix[n_obs, n_obs] Sigma_sim = cholesky_decompose(K_sim);
  vector[n_obs] theta_sim = Sigma_sim * eta_sim;
  int y_sim [n_obs] = poisson_log_rng(log(ye) + theta_sim);
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
}

transformed parameters {
  vector[n_phi] phi;
  phi[1] = alpha;
  phi[2] = rho; 
}

model {
  rho ~ inv_gamma(rho_alpha_prior, rho_beta_prior);
  alpha ~ inv_gamma(10, 10);
  target += laplace_marginal_poisson(y_sim, n_samples, ye, K_functor,
                                     phi, x, delta, delta_int, theta_0);
}

generated quantities {
  int<lower = 0, upper = 1> I_lt_sim[2] = {alpha < alpha_sim, rho < rho_sim};
}