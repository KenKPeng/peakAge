# peakAge

> Hierarchical Bayesian joint modeling of athlete performance trajectories and latent peak ages

`peakAge` is an R package for modeling athlete performance trajectories
and estimating latent peak ages using Bayesian hierarchical models.
The package is designed for athlete-level longitudinal performance data
and supports both **symmetric** and **asymmetric** parametric peak-performance curves.

The modeling framework jointly captures

- individual performance trajectories over age,
- latent athlete-specific peak ages,
- athlete-level covariate effects,
- full Bayesian uncertainty quantification.

This package is intended for researchers and practitioners in
sports science, performance analysis, and applied statistics.

---

## Status

Version 0.2.0 adds the extended asymmetric joint model used in the revised
manuscript while retaining the existing baseline models.

The revised manuscript's asymmetric joint model uses
`model = "extended", curve = "asymmetric"`. It includes correlated athlete
effects, component-specific covariate formulas, competition-level adjustments,
grouped peak-level SDs, and athlete-specific residual SDs with standardized
two-piece Student-t errors. Symmetric trajectories remain available in the
baseline model; the extended model is asymmetric only.

See [Extended asymmetric model](#extended-asymmetric-model) for an example.

---

## Installation

### 1. Install CmdStan (required)

This package uses **CmdStan** via the `cmdstanr` interface.
CmdStan only needs to be installed **once**.

```r
install.packages("cmdstanr")
cmdstanr::install_cmdstan()
```

### 2. Install peakAge

```r
install.packages("remotes")
remotes::install_github("KenKPeng/peakAge")
```

For the fixed manuscript version, use
`remotes::install_github("KenKPeng/peakAge@v0.2.0")`.

Load the package:

```r
library(peakAge)
```

## Quick example (simulated data)

The example below uses simulated athlete data to demonstrate the workflow.

```r
set.seed(1)
# minimal fake data
I <- 8                      # athletes
n_per <- 6                  # comps per athlete
N <- I * n_per

Athlete_id <- rep(paste0("A", seq_len(I)), each = n_per)

# athlete-level covariates
x1_i <- rnorm(I)
x2_i <- rbinom(I, 1, 0.5)

# expand to long
x1 <- rep(x1_i, each = n_per)
x2 <- rep(x2_i, each = n_per)

# ages (long)
age <- unlist(lapply(seq_len(I), function(i) sort(runif(n_per, 18, 32))))

# generate latent parameters
beta0 <- log(26)
beta1 <- 0.15
beta2 <- -0.10
sigma_T <- 0.08

logT_i <- rnorm(I, mean = beta0 + beta1 * x1_i + beta2 * x2_i, sd = sigma_T)
T_i <- exp(logT_i)

p_pop <- 0.5
logk_pop <- log(0.25)
tau_u <- c(5, 0.2)

u_p <- rnorm(I, 0, tau_u[1])
u_k <- rnorm(I, 0, tau_u[2])

p_i <- p_pop + u_p
logk_i <- logk_pop + u_k

sigma_y <- 1.5

mu <- p_i[match(Athlete_id, paste0("A", seq_len(I)))] -
  exp(logk_i[match(Athlete_id, paste0("A", seq_len(I)))]) *
  (age - T_i[match(Athlete_id, paste0("A", seq_len(I)))])^2

y <- rnorm(N, mu, sigma_y)

df <- data.frame(
  Athlete_id = Athlete_id,
  age = age,
  y = y,
  x1 = x1,
  x2 = x2
)

head(df) # check data

# just for sanity check stan_data build, no need for fit model.
sd <- build_stan_data_sym(
  df,
  peak_formula = ~ x1 + x2,
  p_formula = ~ 1,     # no covariates for p
  k_formula = ~ 1,     # no covariates for logk
  grainsize = 50
)
str(sd) # all good

# run a tiny fit (keep it small)
fit <- fit_peakage(
  df,
  curve = "symmetric",
  peak_formula = ~ x1 + x2,
  prior = list(
    p = list(pop = c(0.5, 1)),      # for non-standardized y, mush specify reasonable values
    k = list(pop = c(0, 1.5))         # logk's prior(k = exp(logk))
  ),
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200,
  seed = 123,
  refresh = 1
)

fit$fit$summary(c("p_pop","logk_pop","sigma_y","sigma_T","beta_T")) # for advanced user; can directly check through MCMC sampling

s <- summarise_peakage(fit, format = "data.frame") # build in function
s
s2 <- summarise_peakage(fit, format = "html")
cat(s2$html)

p1 <- plot_athlete_curve(
  fit, df,
  athlete_id = c("A2"),
  point = "mean",
  show_band = T,
  band_level = 0.95,
  band_type = "pred",
  title = "A2: fitted curve (symmetric, mean + 95% band)"
)
print(p1)

# draw estimated tau, p, k
ath <- athlete_estimates(fit,
                         ids = c("A2", "A7"),
                         point = "mean",
                         include_ci = TRUE,
                         include_k = T,
                         include_p = F)




# non-symmetric
fit2 <- fit_peakage(
  df,
  curve = "asymmetric",
  peak_formula = ~ x1 + x2,
  p_formula    = ~ x1,
  k_formula    = ~ x1+x2,
  prior = list(
    p = list(pop = c(0.6, 0.8)),       
    k = list(pop = list(
      L = c(-1.0, 1.0),
      R = c( 0.2, 1.0)
    ))
  ),
  chains = 2,
  parallel_chains = 2,
  iter_warmup = 200,
  iter_sampling = 200,
  seed = 123
)

s <- summarise_peakage(fit2)
print(s)


# plot
p2 <- plot_athlete_curve(
  fit2,
  df,
  athlete_id = c("A1","A4"),
  point = "median",
  show_band = TRUE,
  band_level = 0.95,
  band_type = "pred",
  title = "A2: posterior predictive band (symmetric)",
  gg_opts = list(
    # observed points
    point_size  = 2.2,
    point_alpha = 0.85,
    point_color = "black",
    xlim = c(18, 40),
    # fitted curve
    line_size   = 1.1,
    curve_color = "steelblue",
    
    # predictive band
    band_alpha  = 0.25,
    band_fill   = "steelblue",
    
    # peak-age vertical line
    vline_color    = "red",
    vline_linetype = "dashed",
    
    # overall theme
    theme = ggplot2::theme_bw()
  )
)

print(p2)

# draw estimates
ath2 <- athlete_estimates(fit2,
                         ids = c("A2", "A7"),
                         point = "mean",
                         include_ci = TRUE,
                         include_k = T,
                         include_p = F)
```

## Extended asymmetric model

Use `model = "extended"` with `curve = "asymmetric"` to fit the following
components together.

- Separate athlete-level regressions for peak level, log pre-peak curvature,
  log post-peak curvature, and log peak age. If a side's formula is omitted,
  it uses `k_formula`.
- Correlated athlete deviations in all four quantities. An optional
  athlete-level column `peak_sd_group` defines groups with different
  between-athlete SDs for peak level.
- Competition-level score adjustments via `observation_formula`. Use `~ 0`
  for a trajectory with no competition adjustments.
- Athlete-specific residual SDs and standardized two-piece Student-t errors.
  The two-piece distribution can represent asymmetry, and its unit variance
  makes each athlete's residual SD interpretable on the score scale.

The existing calls above use the baseline model. The example below generates
competition scores and illustrates the extended fitting and prediction workflow.

```r
set.seed(24)
n_athletes <- 16L
n_per <- 12L
athletes <- data.frame(
  Athlete_id = sprintf("S%02d", seq_len(n_athletes)),
  event = factor(rep(c("track", "field"), each = n_athletes / 2L)),
  cohort = rep(c(-0.5, 0.5), times = n_athletes / 2L)
)
row_id <- rep(seq_len(n_athletes), each = n_per)
sim <- athletes[row_id, ]
sim$age <- rep(seq(20, 34, length.out = n_per), times = n_athletes)
sim$indoor <- rbinom(nrow(sim), 1, 0.3)

track <- as.numeric(athletes$event == "track")
height_effect <- rnorm(n_athletes)
peak_level <- 0.8 + 0.05 * track + (0.04 + 0.04 * track) * height_effect
peak_age <- exp(log(27) + 0.03 * track + 0.01 * athletes$cohort +
                  0.05 * (0.5 * height_effect + sqrt(0.75) * rnorm(n_athletes)))
left <- exp(log(0.009) + 0.2 * track + rnorm(n_athletes, sd = 0.15))
right <- exp(log(0.006) + rnorm(n_athletes, sd = 0.15))
residual_sd <- 0.06 * exp(rnorm(n_athletes, sd = 0.3))

# Center and scale a two-piece Student-t draw to mean zero and variance one.
shape <- 0.7
nu <- 8
left_mass <- plogis(-shape)
half_t_mean <- sqrt(nu / pi) * exp(lgamma((nu - 1) / 2) - lgamma(nu / 2))
raw_mean <- 2 * tanh(shape / 2) * half_t_mean
raw_variance <- nu / (nu - 2) * (1 + 3 * tanh(shape / 2)^2) - raw_mean^2
width <- ifelse(runif(nrow(sim)) < left_mass, -2 * left_mass, 2 * (1 - left_mass))
error <- (width * abs(rt(nrow(sim), df = nu)) - raw_mean) / sqrt(raw_variance)
distance <- sim$age - peak_age[row_id]
curvature <- ifelse(distance <= 0, left[row_id], right[row_id])
sim$y <- peak_level[row_id] - curvature * distance^2 +
  0.04 * sim$indoor + residual_sd[row_id] * error

fit_ext <- fit_peakage(
  sim,
  curve = "asymmetric",
  model = "extended",
  peak_formula = ~ event + cohort,
  p_formula = ~ event,
  left_formula = ~ event,
  right_formula = ~ 1,
  observation_formula = ~ 0 + indoor,
  peak_sd_group = "event",
  prior = list(
    beta_mean = list(alpha = c(0.8, 0), L = c(log(0.01), 0), R = log(0.01)),
    beta_sd = list(alpha = 0.2, L = 0.5, R = 0.5, tau = 0.1),
    delta_sd = 0.1,
    scale_sd = c(alpha = 0.15, L = 0.5, R = 0.5, tau = 0.1),
    residual_scale_sd = 0.15
  ),
  chains = 2,
  iter_warmup = 300,
  iter_sampling = 300,
  seed = 24
)

summarise_peakage(fit_ext)
plot_athlete_curve(fit_ext, sim, athlete_id = "S01", band_type = "pred")
```

The formulas for athlete-level parameters include an intercept. Athlete-level
covariates and peak-level SD groups are constant within athlete. Competition
covariates may vary between records. The plot shows the fitted trajectory and
predictive band at zero competition adjustment. The observed scores are adjusted
to that same reference using the fitted competition effects.

The extended model uses fixed priors. Set `prior` for the scale of your scores
and covariates. The default peak-level and residual-scale priors are intended
for scores roughly on a zero-to-one scale. `beta_mean` and `beta_sd` are lists
with entries named `alpha`, `L`, `R`, and `tau`. Each supplied entry can be a
scalar or a vector aligned with that component's design matrix. Named vectors
must match the design column names. See `?fit_peakage` for the full list of
prior arguments.

To reproduce the World Athletics analysis, use its exact design matrices,
contrast coding, birth-year scaling, peak-level SD groups, and two fixed
prior specifications. The synthetic example above illustrates the same model
components with its own design and priors.

The installed example `aoas_v8.R` prepares the full-cohort inputs from
your copy of the original v8 workflow. It supports both fixed prior profiles.

```r
source(system.file("examples", "aoas_v8.R", package = "peakAge"))
inputs <- prepare_aoas_v8("path/to/v8", profile = "primary")
# fit <- fit_aoas_v8("path/to/v8", profile = "primary")
```
