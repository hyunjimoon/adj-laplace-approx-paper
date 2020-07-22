data {
  int <lower = 0> n_obs;
  int <lower = 0> n_covariates;
  vector[n_covariates] x[n_obs];
  vector[n_obs] ye;
}

transformed data{
  real tol = 1e-6;
  real delta = 1e-8;
  int n_phi = 2;
  int n_samples[n_obs] = rep_array(1, n_obs);
  
  real <lower = 0> rho_alpha_prior = 2.42393;
  real <lower = 0> rho_beta_prior = 14.8171;
  real<lower = 0> rho_sim = inv_gamma_rng(rho_alpha_prior, rho_beta_prior);
  real<lower = 0> alpha_sim = inv_gamma_rng(10, 10);
  vector[n_obs] eta_sim = to_vector(normal_rng(rep_vector(0, n_obs), rep_vector(1, n_obs)));
  matrix[n_obs, n_obs] K_sim = cov_exp_quad(x, alpha_sim, rho_sim);
  for (n in 1:n_obs) K_sim[n, n] = K_sim[n,n] + delta;
  matrix[n_obs, n_obs] Sigma_sim = cholesky_decompose(K_sim);
  vector[n_obs] theta_sim = Sigma_sim * eta_sim;
  int y_sim [n_obs] = poisson_log_rng(log(ye) + theta_sim);
}

parameters {
  real<lower = 0> alpha;
  real<lower = 0> rho;
  vector[n_obs] eta;
}

transformed parameters {
   vector[n_obs] theta;
   {
     matrix[n_obs, n_obs] L_Sigma;
     matrix[n_obs, n_obs] Sigma;
     Sigma = cov_exp_quad(x, alpha, rho);
     for (n in 1:n_obs) Sigma[n, n] = Sigma[n,n] + delta;
     L_Sigma = cholesky_decompose(Sigma);
     theta = L_Sigma * eta;
   }
}

model {
  rho ~ inv_gamma(rho_alpha_prior, rho_beta_prior);
  alpha ~ inv_gamma(10, 10);  // CHECK -- is this a good prior?
  eta ~ normal(0, 1);
  y_sim ~ poisson_log(log(ye) + theta);
}

generated quantities {
  int<lower = 0, upper = 1> I_sim_hp[2] = {alpha < alpha_sim, rho < rho_sim};
  int<lower = 0, upper = 1> I_sim_lp[n_obs];
  for (n in 1:n_obs) I_sim_lp[n] = (theta[n] < theta_sim[n]);
}