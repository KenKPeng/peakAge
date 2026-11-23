functions {
  // W has left and right widths 2 inv_logit(-shape) and 2 inv_logit(shape).
  // Z = (W - E[W]) / SD(W) has mean zero and unit variance.
  vector split_t_error_constants(real shape, real nu) {
    vector[7] constants;
    real asymmetry = tanh(0.5 * shape);
    real half_t_mean = sqrt(nu / pi())
      * exp(lgamma(0.5 * (nu - 1)) - lgamma(0.5 * nu));
    real t_second_moment = nu / (nu - 2);
    real variance = t_second_moment
      + square(asymmetry) * (3 * t_second_moment - 4 * square(half_t_mean));
    constants[1] = 2 * asymmetry * half_t_mean;
    constants[2] = sqrt(variance);
    constants[3] = log(2.0) + log_inv_logit(-shape) + 0.5 * log(nu);
    constants[4] = log(2.0) + log_inv_logit(shape) + 0.5 * log(nu);
    constants[5] = exp(-constants[3]);
    constants[6] = exp(-constants[4]);
    constants[7] = lgamma(0.5 * (nu + 1)) - lgamma(0.5 * nu)
      - 0.5 * (log(nu) + log(pi())) + log(constants[2]);
    return constants;
  }

  // Slice athletes to calculate trajectory means and the vectorized t kernel
  // on the same thread. The separate log-domain kernel handles extreme tails.
  real split_t_partial_sum_lpdf(array[] vector theta_slice, int start, int end,
      vector age, vector y, array[] int n_obs, array[] int first_obs,
      array[] int observation_index, vector observation_shift, real nu,
      vector constants, vector inverse_width) {
    vector[sum(n_obs[start:end])] scaled;
    int n_regular = 0;
    real lp = 0;
    for (j in 1:size(theta_slice)) {
      int i = start + j - 1;
      if (n_obs[i] > 0) {
        real alpha = theta_slice[j, 1];
        real rho_L = theta_slice[j, 2];
        real rho_R = theta_slice[j, 3];
        real tau = theta_slice[j, 4];
        real inv_sigma = theta_slice[j, 6];
        lp -= n_obs[i] * theta_slice[j, 5];
        for (k in 0:(n_obs[i] - 1)) {
          int n = observation_index[first_obs[i] + k];
          real rho = age[n] <= tau ? rho_L : rho_R;
          real mu = alpha - rho * square(age[n] - tau) + observation_shift[n];
          real w = (y[n] - mu) * inv_sigma * constants[2] + constants[1];
          int side = w < 0 ? 1 : 2;
          real x = w * inverse_width[side];
          // Avoid overflow in the built-in t kernel and its derivatives.
          if (!is_nan(x) && abs(x) <= 1e140 && nu <= 1e140) {
            n_regular += 1;
            scaled[n_regular] = x;
          } else {
            real scaled_square = square(w * constants[4 + side]);
            real log_kernel = log1p(scaled_square);
            if (is_inf(scaled_square) || is_nan(scaled_square)) {
              log_kernel = w == 0 ? 0
                : log1p_exp(2 * (log(abs(w)) - constants[2 + side]));
            }
            lp += constants[7] - 0.5 * (nu + 1) * log_kernel;
          }
        }
      }
    }
    if (n_regular > 0) {
      // The two-piece law has no inverse-width factor. Standardization
      // contributes log(SD(W)) per record, and log(sigma) was subtracted above.
      lp += student_t_lpdf(scaled[1:n_regular] | nu, 0, 1)
        + n_regular * log(constants[2]);
    }
    return lp;
  }
}

data {
  int<lower=1> I;
  int<lower=1> N;
  array[4] int<lower=1> K;
  matrix[I, K[1]] X_alpha;
  matrix[I, K[2]] X_L;
  matrix[I, K[3]] X_R;
  matrix[I, K[4]] X_tau;
  int<lower=0> P;
  matrix[N, P] Z;
  vector[N] y;
  vector[N] age;
  array[N] int<lower=1, upper=I> id;
  int<lower=1> G;
  array[I] int<lower=1, upper=G> group;

  vector[K[1]] beta_prior_mean_alpha;
  vector[K[2]] beta_prior_mean_L;
  vector[K[3]] beta_prior_mean_R;
  vector[K[4]] beta_prior_mean_tau;
  vector<lower=0>[K[1]] beta_prior_sd_alpha;
  vector<lower=0>[K[2]] beta_prior_sd_L;
  vector<lower=0>[K[3]] beta_prior_sd_R;
  vector<lower=0>[K[4]] beta_prior_sd_tau;
  vector<lower=0>[P] delta_prior_sd;
  vector<lower=0>[4] scale_prior_sd;
  real<lower=0> residual_scale_prior_sd;
  real<lower=0> residual_heterogeneity_prior_sd;
  real<lower=0> lkj_eta;
  real<lower=0> error_shape_prior_sd;
  real<lower=0> error_nu_prior_shape;
  real<lower=0> error_nu_prior_rate;
  int<lower=1> grainsize;
}

transformed data {
  array[I] int n_obs = rep_array(0, I);
  array[I] int first_obs;
  array[N] int observation_index;
  vector[I] age_min = rep_vector(0, I);
  vector[I] age_max = rep_vector(0, I);
  for (n in 1:N) {
    int i = id[n];
    if (n_obs[i] == 0) {
      age_min[i] = age[n];
      age_max[i] = age[n];
    } else {
      age_min[i] = fmin(age_min[i], age[n]);
      age_max[i] = fmax(age_max[i], age[n]);
    }
    n_obs[i] += 1;
  }
  {
    int position = 1;
    array[I] int next_obs;
    for (i in 1:I) {
      first_obs[i] = position;
      next_obs[i] = position;
      position += n_obs[i];
    }
    for (n in 1:N) {
      int i = id[n];
      observation_index[next_obs[i]] = n;
      next_obs[i] += 1;
    }
  }
}

parameters {
  vector[K[1]] b_alpha;
  vector[K[2]] b_L;
  vector[K[3]] b_R;
  vector[K[4]] b_tau;
  vector[P] delta;

  vector<lower=0>[G] sigma_alpha_group;
  real<lower=0> sigma_L;
  real<lower=0> sigma_R;
  real<lower=0> sigma_tau;
  cholesky_factor_corr[4] L_Omega;
  matrix[4, I] z_b;

  real<lower=1e-6> sigma_y;
  real<lower=0> omega_y;
  vector[I] z_sigma;
  real error_shape;
  real<lower=2> error_nu;
}

transformed parameters {
  // The common correlation order is peak level, log left curvature,
  // log right curvature, and log peak age.
  vector[I] p_hat;
  vector[I] logkL_hat;
  vector[I] logkR_hat;
  vector[I] T_hat;
  vector[I] sigma_y_i = sigma_y * exp(omega_y * z_sigma);
  {
    matrix[4, I] b = diag_pre_multiply([1, sigma_L, sigma_R, sigma_tau]', L_Omega) * z_b;
    for (i in 1:I) {
      b[1, i] *= sigma_alpha_group[group[i]];
    }
    p_hat = X_alpha * b_alpha + b[1, ]';
    logkL_hat = X_L * b_L + b[2, ]';
    logkR_hat = X_R * b_R + b[3, ]';
    T_hat = exp(X_tau * b_tau + b[4, ]');
  }
}

model {
  vector[I] rho_L = rep_vector(0, I);
  vector[I] rho_R = rep_vector(0, I);
  vector[N] observation_shift = rep_vector(0, N);
  vector[7] constants = split_t_error_constants(error_shape, error_nu);
  vector[2] inverse_width = sqrt(error_nu) * constants[5:6];
  vector[I] log_sigma_y_i = log(sigma_y) + omega_y * z_sigma;
  vector[I] inv_sigma_y_i = exp(-log_sigma_y_i);
  array[I] vector[6] athlete_theta;

  b_alpha ~ normal(beta_prior_mean_alpha, beta_prior_sd_alpha);
  b_L ~ normal(beta_prior_mean_L, beta_prior_sd_L);
  b_R ~ normal(beta_prior_mean_R, beta_prior_sd_R);
  b_tau ~ normal(beta_prior_mean_tau, beta_prior_sd_tau);
  if (P > 0) {
    delta ~ normal(0, delta_prior_sd);
    observation_shift = Z * delta;
  }

  sigma_alpha_group ~ normal(0, scale_prior_sd[1]);
  sigma_L ~ normal(0, scale_prior_sd[2]);
  sigma_R ~ normal(0, scale_prior_sd[3]);
  sigma_tau ~ normal(0, scale_prior_sd[4]);
  L_Omega ~ lkj_corr_cholesky(lkj_eta);
  to_vector(z_b) ~ std_normal();

  sigma_y ~ normal(0, residual_scale_prior_sd);
  omega_y ~ normal(0, residual_heterogeneity_prior_sd);
  z_sigma ~ std_normal();
  error_shape ~ normal(0, error_shape_prior_sd);
  // The Gamma prior is truncated to nu > 2 by the parameter bound.
  error_nu ~ gamma(error_nu_prior_shape, error_nu_prior_rate);

  for (i in 1:I) {
    // Avoid evaluating the curvature exponential on an unused side.
    if (n_obs[i] > 0) {
      if (age_min[i] <= T_hat[i]) rho_L[i] = exp(logkL_hat[i]);
      if (age_max[i] > T_hat[i]) rho_R[i] = exp(logkR_hat[i]);
    }
    athlete_theta[i] = [p_hat[i], rho_L[i], rho_R[i], T_hat[i],
                       log_sigma_y_i[i], inv_sigma_y_i[i]]';
  }
  target += reduce_sum(split_t_partial_sum_lupdf, athlete_theta, grainsize,
    age, y, n_obs, first_obs, observation_index, observation_shift,
    error_nu, constants, inverse_width);
}

generated quantities {
  corr_matrix[4] Omega = multiply_lower_tri_self_transpose(L_Omega);
}
