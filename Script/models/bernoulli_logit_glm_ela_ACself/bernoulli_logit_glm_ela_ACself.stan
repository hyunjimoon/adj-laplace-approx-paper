
functions {
  matrix K_functor (vector parm, matrix X,
                    real[] scale_icept, int[] delta_int) {
    int d = delta_int[1];
    vector[d + 1] k_diag;
    k_diag[1] = scale_icept[1];
    k_diag[2:(d + 1)] = parm[1] * parm[2:(d + 1)];

    return X * diag_pre_multiply(k_diag, X');
  }
}

data {
  int<lower=0> n_obs;				      // number of observations
  int<lower=0> n_covariates;             // number of predictors
  matrix[n_obs,n_covariates] x;
  real<lower=0> scale_icept;    // prior std for the intercept
  real<lower=0> scale_global;	// scale for the half-t prior for tau
  real<lower=1> nu_global;	  // degrees of freedom for the half-t priors for tau
  real<lower=1> nu_local;		  // degrees of freedom for the half-t priors for lambdas
                              // (nu_local = 1 corresponds to the horseshoe)
  real<lower=0> slab_scale;   // for the regularized horseshoe
  real<lower=0> slab_df;
}

transformed data {
  real beta0;
  int delta_int[1] = {n_covariates};
  real delta[1] = {scale_icept};
  int<lower=0,upper=1> y[n_obs];	// outputs
  matrix[n_obs, n_covariates + 1] X;  // design matrix with intercept
  X[, 1] = rep_vector(1.0, n_obs);
  X[, 2:(n_covariates + 1)] = x;
  int n_samples[n_obs] = rep_array(1, n_obs);
  int n_par = 3;                      
  vector[n_obs] theta0 = rep_vector(0, n_obs);
  vector[n_covariates] z_;              // for non-centered parameterization
  real <lower=0> tau_ = abs(student_t_rng(nu_global, 0, scale_global*2));      // global shrinkage parameter
  vector <lower=0> [n_covariates] lambda_;
  for (i in 1:n_covariates) lambda_[i] = abs(student_t_rng(nu_local, 0, 1));
  real<lower=0> caux_ = inv_gamma_rng(0.5*slab_df, 0.5*slab_df);
  real beta0_ = normal_rng(0, scale_icept);
  for (i in 1:n_covariates) z_[i]  = normal_rng(0,1);
  vector[n_covariates] beta_;                     // regression coefficients
  {
    vector[n_covariates] lambda_tilde_;   // 'truncated' local shrinkage parameter
    real c_ = slab_scale * sqrt(caux_); // slab scale
    //for (i in 1:n_covariates) lambda_tilde_[i] = sqrt( c_^2 * square(lambda_[i]) / (c_^2 + tau_^2*square(lambda_[i])));
    lambda_tilde_ = sqrt( c_^2 * square(lambda_) ./ (c_^2 + tau_^2*square(lambda_)));
    beta_ = z_ .* lambda_tilde_*tau_;
  }
  y = bernoulli_logit_rng(beta0_ + x * beta_);
}

parameters {
  real <lower=0> tau;         // global shrinkage parameter
  vector <lower=0>[n_covariates] lambda; // local shrinkage parameter
  real<lower=0> caux;
}

transformed parameters {
  vector[n_covariates + 1] parm;        // regression coefficients including icept
  {
    vector[n_covariates] lambda_tilde;   // 'truncated' local shrinkage parameter
    real c = slab_scale * sqrt(caux); // slab scale
    //for (i in 1:n_covariates) lambda_tilde[i] = sqrt( c^2 * square(lambda[i]) / (c^2 + tau^2*square(lambda[i])));
    lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)));
    parm[1] = tau;
    parm[2:(n_covariates + 1)] = lambda_tilde;
  }
}

model {
  // half-t priors for lambdas and tau, and inverse-gamma for c^2theta0v
  lambda ~ student_t(nu_local, 0, 1);
  tau ~ student_t(nu_global, 0, scale_global*2);
  caux ~ inv_gamma(0.5*slab_df, 0.5*slab_df);

  target += laplace_marginal_bernoulli(y, n_samples, K_functor, parm,
                                       X, delta, delta_int, theta0);
}

generated quantities {
  int y_[n_obs] = y;
  vector[n_par] pars_;
  pars_[1] = caux_;
  pars_[2] = tau_;
  pars_[3] = lambda_[2586];
  vector[n_obs] f = laplace_approx_bernoulli_rng(y, n_samples, K_functor, parm,
                                             X, delta, delta_int, theta0);
  int ranks_[n_par] = {caux < caux_, tau < tau_, lambda[2586] < lambda_[2586]};
}
