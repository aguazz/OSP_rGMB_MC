##############################
# data_generator.R
#
# This script runs multiple test scenarios ("cases") that illustrate how to:
#   - Set up parameters for various Gauss–Markov processes,
#   - Solve an optimal stopping problem (OSP) using osp_solver(...),
#   - Save the results as .RData files in a 'data' directory.
#
# Each "case" below defines its own parameter list, calls osp_solver(...),
# then saves the resulting solution object. Adjust as needed.
##############################

# 1) Ensure 'data' directory exists
if(!dir.exists("data")) dir.create("data")

# 2) Source the main GMP/OSP functions
#    This file must provide: gmp_integrals, simulate_GMP_bridge,
#                           sample_z_discrete, sample_z_cont (MH), osp_solver
source("functions.R")


##############################
# CASE 1: TWO-POINTS BROWNIAN BRIDGE
##############################
# In this scenario, we assume a standard Brownian Motion (beta=0, alpha=0, zeta=1)
# but with a discrete prior on X(T) taking values { -1, +1 } each with probability 0.5.
# We vary the initial condition x0 from -1 to +1 in increments of 0.05, solve the OSP,
# and save each solution.

x0_values <- seq(-1, 1, by = 0.05)
parameter_list <- lapply(x0_values, function(x0) {
  list(
    t_grid = log(seq(1, exp(1), length.out = 1000)), 
    x_grid = seq(-3, 3, length.out = 1000), 
    slope_fun = function(t) rep(0, length(t)), 
    trend_fun = function(t) rep(0, length(t)), 
    vol_fun   = function(t) rep(1, length(t)), 
    n_gmp_grid = 5000, 
    x0 = x0, 
    discrete_prior = TRUE, 
    support = c(-1, 1), 
    probs = c(0.5, 0.5), 
    M_samples = 10000,
    proposal_sd = 1, 
    burn_in = 1000
  )
})

for (i in seq_along(parameter_list)) {
  params <- parameter_list[[i]]
  
  # Solve the OSP using the discrete prior approach
  osp_solution <- osp_solver(
    t_grid         = params$t_grid, 
    x_grid         = params$x_grid, 
    slope_fun      = params$slope_fun, 
    trend_fun      = params$trend_fun, 
    vol_fun        = params$vol_fun, 
    n_gmp_grid     = params$n_gmp_grid, 
    x0             = params$x0, 
    discrete_prior = params$discrete_prior, 
    support        = params$support, 
    probs          = params$probs, 
    M_samples      = params$M_samples,
    proposal_sd    = params$proposal_sd, 
    burn_in        = params$burn_in
  )
  
  # Save solution as RData file
  file_name <- sprintf("%03d_osp_solution_2pBB_x0_%.2f.RData", i, params$x0)
  save(osp_solution, file = file.path("data", file_name))
  cat("Saved solution to", file_name, "\n")
}


##############################################
# CASE 2: SINUSOIDAL MEAN REVERSION
##############################################
#   slope(t) = -1 (pulling force)
#   trend(t) = 2*sin(10*pi*t)  (oscillatory)
#   vol(t)   = 1
# We use a *truncated normal* prior for X(T).
# The prior's mean = 0, and std dev is some fraction of the final-time variance.

params <- list(
  regime = "sin_mean_rev",
  t_grid = log(seq(1, exp(1), length.out = 1000)),
  x_grid = seq(-1, 3, length.out = 1000),
  slope_fun = function(t) rep(-1, length(t)),
  trend_fun = function(t) 2*sin(10*pi*t),
  vol_fun   = function(t) rep(1, length(t)),
  n_gmp_grid = 20000,
  x0 = 0,
  discrete_prior = FALSE,    # continuous prior sampling
  initial_val = 0.2,
  M_samples = 10000,
  proposal_sd = 1,
  burn_in = 1000
)

# Compute final-time variance from the GMP integrals
last_var <- var_gmp(
  t0 = 0, t1 = max(params$t_grid),
  gmp = gmp_integrals(params$slope_fun, params$trend_fun, params$vol_fun,
                          T = max(params$t_grid), n_grid = 5000)
)

# Truncated-normal prior: [0, +Inf)
params$prior_type <- "truncated_normal"
params$prior_params <- list(mu = 0, sd = 0.25*sqrt(last_var), lower = 0, upper = Inf)

# Solve the OSP
osp_solution <- osp_solver(
  t_grid         = params$t_grid,
  x_grid         = params$x_grid,
  slope_fun      = params$slope_fun,
  trend_fun      = params$trend_fun,
  vol_fun        = params$vol_fun,
  n_gmp_grid     = params$n_gmp_grid,
  x0             = params$x0,
  discrete_prior = params$discrete_prior,
  prior_type     = params$prior_type,
  prior_params   = params$prior_params,
  M_samples      = params$M_samples,
  burn_in        = params$burn_in,
  proposal_sd    = params$proposal_sd,
  initial_val    = params$initial_val
)

file_name <- sprintf("osp_solution_%s_x0_%.2f.RData", params$regime, params$x0)
save(osp_solution, file = file.path("data", file_name))
cat("Saved solution to", file_name, "\n")


##############################################
# CASE 3: RAPID SWITCHING OF PULLING STRENGTH
##############################################
#   slope(t) = slope_1 + (slope_2 - slope_1)/2 * (1 + tanh(k*(t - t0)))
#   trend(t) = 0
#   vol(t)   = 1
# Again, we use a truncated-normal prior on [0, +Inf).

params <- list(
  regime = "slope_switch",
  t_grid = log(seq(1, exp(1), length.out = 1000)),
  x_grid = seq(-1, 1, length.out = 1000),
  slope_fun = function(t) {
    slope_1 <- -10
    slope_2 <- -0.5
    k       <- 100
    t0      <- 0.5
    slope_1 + (slope_2 - slope_1)/2 * (1 + tanh(k*(t - t0)))
  },
  trend_fun  = function(t) rep(0, length(t)),
  vol_fun    = function(t) rep(1, length(t)),
  n_gmp_grid = 20000,
  x0 = 0,
  discrete_prior = FALSE,
  initial_val = 0.2,
  M_samples = 10000,
  proposal_sd = 1,
  burn_in = 1000
)

# Compute final-time variance for the prior
last_var <- var_gmp(
  t0 = 0, t1 = max(params$t_grid),
  gmp = gmp_integrals(params$slope_fun, params$trend_fun, params$vol_fun,
                          T = max(params$t_grid), n_grid = 5000)
)
params$prior_type <- "truncated_normal"
params$prior_params <- list(mu = 0, sd = 0.25*sqrt(last_var), lower = 0, upper = Inf)

# Solve and save
osp_solution <- osp_solver(
  t_grid         = params$t_grid,
  x_grid         = params$x_grid,
  slope_fun      = params$slope_fun,
  trend_fun      = params$trend_fun,
  vol_fun        = params$vol_fun,
  n_gmp_grid     = params$n_gmp_grid,
  x0             = params$x0,
  discrete_prior = params$discrete_prior,
  prior_type     = params$prior_type,
  prior_params   = params$prior_params,
  M_samples      = params$M_samples,
  burn_in        = params$burn_in,
  proposal_sd    = params$proposal_sd,
  initial_val    = params$initial_val
)

file_name <- sprintf("osp_solution_%s_x0_%.2f.RData", params$regime, params$x0)
save(osp_solution, file = file.path("data", file_name))
cat("Saved solution to", file_name, "\n")


##############################################
# CASE 4: PARABOLIC VOLATILITY SMILE
##############################################
#   slope(t) = -1
#   trend(t) = 0
#   vol(t)   = v0 + (a*(t - t0))^k
#   => demonstrates time-varying volatility
##############################################

params <- list(
  regime = "vol_smile",
  t_grid = log(seq(1, exp(1), length.out = 1000)),
  x_grid = seq(-8, 8, length.out = 1000),
  slope_fun = function(t) rep(-1, length(t)),
  trend_fun = function(t) rep(0, length(t)),
  vol_fun = function(t) {
    t0 <- 0.5
    v0 <- 0.25
    k  <- 4
    a  <- 4
    v0 + (a*(t - t0))^k
  },
  n_gmp_grid = 20000,
  x0 = 0,
  discrete_prior = FALSE,
  initial_val = 0.2,
  M_samples = 10000,
  proposal_sd = 1,
  burn_in = 1000
)

# Final-time variance
last_var <- var_gmp(
  t0 = 0, t1 = max(params$t_grid),
  gmp = gmp_integrals(params$slope_fun, params$trend_fun, params$vol_fun,
                          T = max(params$t_grid), n_grid = 5000)
)
params$prior_type <- "truncated_normal"
params$prior_params <- list(mu = 0, sd = 0.25*sqrt(last_var), lower = 0, upper = Inf)

osp_solution <- osp_solver(
  t_grid         = params$t_grid,
  x_grid         = params$x_grid,
  slope_fun      = params$slope_fun,
  trend_fun      = params$trend_fun,
  vol_fun        = params$vol_fun,
  n_gmp_grid     = params$n_gmp_grid,
  x0             = params$x0,
  discrete_prior = params$discrete_prior,
  prior_type     = params$prior_type,
  prior_params   = params$prior_params,
  M_samples      = params$M_samples,
  burn_in        = params$burn_in,
  proposal_sd    = params$proposal_sd,
  initial_val    = params$initial_val
)

file_name <- sprintf("osp_solution_%s_x0_%.2f.RData", params$regime, params$x0)
save(osp_solution, file = file.path("data", file_name))
cat("Saved solution to", file_name, "\n")


##############################################
# CASE 5: GREATER FINAL VARIANCE
##############################################
#   slope(t) = 0
#   trend(t) = 0
#   vol(t)   = 1 (or user can modify)
#   => we use a normal prior with an inflated variance factor.
##############################################

params <- list(
  regime = "greater_final_var",
  t_grid = log(seq(1, exp(1), length.out = 1000)),
  x_grid = seq(-20, 3, length.out = 1000),
  slope_fun = function(t) rep(0, length(t)),
  trend_fun = function(t) rep(0, length(t)),
  vol_fun   = function(t) rep(1, length(t)),
  n_gmp_grid = 20000,
  x0 = 0,
  discrete_prior = FALSE,
  initial_val = 0.2,
  M_samples = 10000,
  proposal_sd = 1,
  burn_in = 1000
)

# Compute final-time variance
last_var <- var_gmp(
  t0 = 0, t1 = max(params$t_grid),
  gmp = gmp_integrals(params$slope_fun, params$trend_fun, params$vol_fun,
                          T = max(params$t_grid), n_grid = 5000)
)
# Use normal prior with variance 1.02 * last_var
params$prior_type <- "normal"
params$prior_params <- list(mu = 0, sd = sqrt(1*last_var), lower = -Inf, upper = Inf)

osp_solution <- osp_solver(
  t_grid         = params$t_grid,
  x_grid         = params$x_grid,
  slope_fun      = params$slope_fun,
  trend_fun      = params$trend_fun,
  vol_fun        = params$vol_fun,
  n_gmp_grid     = params$n_gmp_grid,
  x0             = params$x0,
  discrete_prior = params$discrete_prior,
  prior_type     = params$prior_type,
  prior_params   = params$prior_params,
  M_samples      = params$M_samples,
  burn_in        = params$burn_in,
  proposal_sd    = params$proposal_sd,
  initial_val    = params$initial_val
)

file_name <- sprintf("osp_solution_%s_x0_%.2f.RData", params$regime, params$x0)
save(osp_solution, file = file.path("data", file_name))
cat("Saved solution to", file_name, "\n")
