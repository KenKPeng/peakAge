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

🚧 **Under active development**  
The API may change as new features are added.
Feedback and contributions are welcome.

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
