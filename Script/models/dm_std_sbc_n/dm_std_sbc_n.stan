data {
  int <lower = 0> n_obs;
  int <lower = 0> n_covariates;
  vector[n_covariates] x[n_obs];
  vector[n_obs] ye;
  real <lower = 0> alpha_mu_prior;
  real <lower = 0> alpha_sd_prior;
  real <lower = 0> rho_mu_prior;
  real <lower = 0> rho_sd_prior;
}

transformed data{
  real tol = 1e-6;
  real delta = 1e-8;
  int n_phi = 2;
  int n_par = 3;
  real<lower = 0> alpha_ = normal_rng(alpha_mu_prior, alpha_sd_prior);
  real<lower = 0> rho_ = normal_rng(rho_mu_prior, rho_sd_prior);
  vector[n_obs] eta_ = to_vector(normal_rng(rep_vector(0, n_obs), rep_vector(1, n_obs)));
  matrix[n_obs, n_obs] K_ = cov_exp_quad(x, alpha_, rho_);
  for (n in 1:n_obs) K_[n, n] = K_[n,n] + delta;
  matrix[n_obs, n_obs] Sigma_ = cholesky_decompose(K_);
  vector[n_obs] theta_ = Sigma_ * eta_;
  int y [n_obs] = poisson_log_rng(log(ye) + theta_);
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
  alpha ~ normal(alpha_mu_prior, alpha_sd_prior); 
  rho ~ normal(rho_mu_prior, rho_sd_prior);
  eta ~ normal(0, 1);
  y ~ poisson_log(log(ye) + theta);
}

generated quantities {
  int y_[n_obs] = y;
  vector[n_par] pars_;
  // vector[n_obs] log_lik;
  pars_[1] = alpha_;
  pars_[2] = rho_;
  pars_[3] = theta_[1];
  int ranks_[n_par] = {alpha < alpha_, rho < rho_, theta[1] < theta_[1]};
}
