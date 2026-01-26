functions {
  real performance_partial_sum_asym(
    array[] real y_slice,
    int start, int end,
    array[] int id_full,
    array[] real age_full,
    vector p_i,
    vector logkL_i,
    vector logkR_i,
    vector T_i,
    real sigma_y
  ) {
    real lp = 0;
    int N_slice = end - start + 1;
    for (n in 1:N_slice) {
      int idx = start + n - 1;
      int i = id_full[idx];
      real a = age_full[idx];
      real d2 = square(a - T_i[i]);
      real k = (a <= T_i[i]) ? exp(logkL_i[i]) : exp(logkR_i[i]);
      real mu = p_i[i] - k * d2;
      lp += normal_lpdf(y_slice[n] | mu, sigma_y);
    }
    return lp;
  }
}

data {
  int<lower=1> N;
  int<lower=1> I;
  array[N] int<lower=1, upper=I> id;
  vector[N] age;
  array[N] real y;

  int<lower=1> grainsize;

  // AFT (logT) covariates
  int<lower=1> K_T;
  matrix[I, K_T] X_T;

  // p covariates (allow none)
  int<lower=0> K_p;
  matrix[I, K_p] X_p;

  // curvature covariates (allow none)
  int<lower=0> K_k;
  matrix[I, K_k] X_k;

  // -------- priors (hyperparameters provided by R) --------
  real p_pop_loc;
  real<lower=0> p_pop_scale;

  real logkL_pop_loc;
  real<lower=0> logkL_pop_scale;

  real logkR_pop_loc;
  real<lower=0> logkR_pop_scale;

  vector<lower=0>[3] tau_u_rate;  // exponential rates for tau_u

  real<lower=0> sigma_y_scale;

  real<lower=0> beta_T_scale;
  real<lower=0> beta_p_scale;
  real<lower=0> beta_kL_scale;
  real<lower=0> beta_kR_scale;

  real<lower=0> sigma_T_scale;
}

transformed data {
  array[N] real age_array;
  for (n in 1:N) age_array[n] = age[n];
}

parameters {
  real p_pop;
  real logkL_pop;
  real logkR_pop;

  vector<lower=0>[3] tau_u;
  array[I] vector[3] z_u;

  real<lower=1e-6> sigma_y;

  vector[I] logT_i;
  vector[K_T] beta_T;
  real<lower=0> sigma_T;

  vector[K_p] beta_p;
  vector[K_k] beta_kL;
  vector[K_k] beta_kR;
}

transformed parameters {
  vector[I] p_i;
  vector[I] logkL_i;
  vector[I] logkR_i;
  vector[I] T_i;
  vector[I] mu_logT;

  vector[I] fe_p  = X_p * beta_p;
  vector[I] fe_kL = X_k * beta_kL;
  vector[I] fe_kR = X_k * beta_kR;

  for (i in 1:I) {
    real u_p  = tau_u[1] * z_u[i][1];
    real u_kL = tau_u[2] * z_u[i][2];
    real u_kR = tau_u[3] * z_u[i][3];

    p_i[i]     = p_pop     + fe_p[i]  + u_p;
    logkL_i[i] = logkL_pop + fe_kL[i] + u_kL;
    logkR_i[i] = logkR_pop + fe_kR[i] + u_kR;
  }

  T_i = exp(logT_i);
  mu_logT = X_T * beta_T;
}

model {
  // priors: baselines
  p_pop     ~ normal(p_pop_loc, p_pop_scale);
  logkL_pop ~ normal(logkL_pop_loc, logkL_pop_scale);
  logkR_pop ~ normal(logkR_pop_loc, logkR_pop_scale);

  // priors: random effects
  tau_u[1] ~ exponential(tau_u_rate[1]);
  tau_u[2] ~ exponential(tau_u_rate[2]);
  tau_u[3] ~ exponential(tau_u_rate[3]);
  for (i in 1:I) z_u[i] ~ std_normal();

  // priors: noise
  sigma_y ~ normal(0, sigma_y_scale);

  // priors: regressions
  beta_T   ~ normal(0, beta_T_scale);
  beta_p   ~ normal(0, beta_p_scale);
  beta_kL  ~ normal(0, beta_kL_scale);
  beta_kR  ~ normal(0, beta_kR_scale);
  sigma_T  ~ normal(0, sigma_T_scale);

  // AFT likelihood
  logT_i ~ normal(mu_logT, sigma_T);

  // performance likelihood
  target += reduce_sum(
    performance_partial_sum_asym,
    y,
    grainsize,
    id,
    age_array,
    p_i,
    logkL_i,
    logkR_i,
    T_i,
    sigma_y
  );
}

generated quantities {
  vector[I] T_hat = exp(logT_i);
  vector[I] p_hat = p_i;
  vector[I] logkL_hat = logkL_i;
  vector[I] logkR_hat = logkR_i;
}
