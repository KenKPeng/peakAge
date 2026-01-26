functions {
  real performance_partial_sum(
    array[] real y_slice,
    int start, int end,
    array[] int id_full,
    array[] real age_full,
    vector p_i,
    vector logk_i,
    vector T_i,
    real sigma_y
  ) {
    real lp = 0;
    int N_slice = end - start + 1;
    for (n in 1:N_slice) {
      int idx = start + n - 1;
      int i = id_full[idx];
      real mu_n = p_i[i] - exp(logk_i[i]) * square(age_full[idx] - T_i[i]);
      lp += normal_lpdf(y_slice[n] | mu_n, sigma_y);
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

  // performance covariates (allow none)
  int<lower=0> K_p;
  matrix[I, K_p] X_p;

  int<lower=0> K_k;
  matrix[I, K_k] X_k;

  // -------- priors (hyperparameters provided by R) --------
  real p_pop_loc;
  real<lower=0> p_pop_scale;

  real logk_pop_loc;
  real<lower=0> logk_pop_scale;

  vector<lower=0>[2] tau_u_rate;  // exponential rates for tau_u

  real<lower=0> sigma_y_scale;

  real<lower=0> beta_T_scale;
  real<lower=0> beta_p_scale;
  real<lower=0> beta_k_scale;

  real<lower=0> sigma_T_scale;
}

transformed data {
  array[N] real age_array;
  for (n in 1:N) age_array[n] = age[n];
}

parameters {
  // performance baselines
  real p_pop;
  real logk_pop;

  // athlete random effects
  vector<lower=0>[2] tau_u;
  array[I] vector[2] z_u;

  real<lower=1e-6> sigma_y;

  // AFT for log T
  vector[I] logT_i;
  vector[K_T] beta_T;
  real<lower=0> sigma_T;

  // covariate effects for p and logk
  vector[K_p] beta_p;
  vector[K_k] beta_k;
}

transformed parameters {
  vector[I] p_i;
  vector[I] logk_i;
  vector[I] T_i;
  vector[I] mu_logT;

  vector[I] mu_p;
  vector[I] mu_logk;

  mu_p    = X_p * beta_p;
  mu_logk = X_k * beta_k;

  for (i in 1:I) {
    p_i[i]    = p_pop    + mu_p[i]    + tau_u[1] * z_u[i][1];
    logk_i[i] = logk_pop + mu_logk[i] + tau_u[2] * z_u[i][2];
  }

  T_i = exp(logT_i);
  mu_logT = X_T * beta_T;
}

model {
  // priors: performance baselines
  p_pop    ~ normal(p_pop_loc, p_pop_scale);
  logk_pop ~ normal(logk_pop_loc, logk_pop_scale);

  // priors: random effects
  tau_u[1] ~ exponential(tau_u_rate[1]);
  tau_u[2] ~ exponential(tau_u_rate[2]);
  for (i in 1:I) z_u[i] ~ std_normal();

  // priors: noise
  sigma_y  ~ normal(0, sigma_y_scale);

  // priors: regressions
  beta_T   ~ normal(0, beta_T_scale);
  beta_p   ~ normal(0, beta_p_scale);
  beta_k   ~ normal(0, beta_k_scale);
  sigma_T  ~ normal(0, sigma_T_scale);

  // AFT likelihood
  logT_i ~ normal(mu_logT, sigma_T);

  // performance likelihood
  target += reduce_sum(
    performance_partial_sum,
    y,
    grainsize,
    id,
    age_array,
    p_i,
    logk_i,
    T_i,
    sigma_y
  );
}

generated quantities {
  vector[I] T_hat = T_i;
  vector[I] p_hat = p_i;
  vector[I] logk_hat = logk_i;
}
