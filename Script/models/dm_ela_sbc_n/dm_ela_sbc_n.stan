
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
  real <lower = 0> alpha_mu_prior;
  real <lower = 0> alpha_sd_prior;
  real <lower = 0> rho_mu_prior;
  real <lower = 0> rho_sd_prior;
}

transformed data{
  real tol = 1e-8;
  real delta[0];
  int delta_int[1] = {n_obs};
  vector[n_obs] theta_0 = rep_vector(0, n_obs);
  int n_phi = 2;
  int n_par = 3;
  int n_samples[n_obs] = rep_array(1, n_obs);
  
  real<lower = 0> alpha_ = abs(normal_rng(alpha_mu_prior, alpha_sd_prior));
  real<lower = 0> rho_ = abs(normal_rng(rho_mu_prior, rho_sd_prior));
  vector[n_obs] eta_ = to_vector(normal_rng(rep_vector(0, n_obs), rep_vector(1, n_obs)));
  matrix[n_obs, n_obs] K_ = cov_exp_quad(x, alpha_, rho_);
  matrix[n_obs, n_obs] Sigma_;
  vector[n_obs] theta_;
  int y [n_obs];
  for (n in 1:n_obs) K_[n, n] = K_[n,n] + tol;
  Sigma_ = cholesky_decompose(K_);
  theta_ = Sigma_ * eta_;
  y = poisson_log_rng(log(ye) + theta_);
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
  rho ~ normal(rho_mu_prior, rho_sd_prior);
  alpha ~ normal(alpha_mu_prior, alpha_sd_prior);

  target += laplace_marginal_poisson(y, n_samples, ye, K_functor,
                                     phi, x, delta, delta_int, theta_0);
}

generated quantities {
  int y_[n_obs] = y;
  vector[n_par] pars_;
  
  // vector[n_obs] log_lik;
  pars_[1] = alpha_;
  pars_[2] = rho_;
  pars_[3] = theta_[1];
  vector[n_obs] theta = laplace_approx_poisson_rng(y, n_samples, ye, K_functor,
                               phi, x, delta, delta_int, theta_0);
  int ranks_[n_par] = {alpha < alpha_, rho < rho_, theta[1] < theta_[1]};
  //print("alpha_", alpha_, "rho_", rho_);
  // for (n in 1:n_obs) log_lik[n] = laplace_marginal_poisson(y, n_samples, ye, K_functor
                                    // phi, x, delta, delta_int, theta_0); 
}