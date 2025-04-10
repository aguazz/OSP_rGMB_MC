############################################################################
# functions.R
#
# This script contains functions for:
#  1) Building splines for partial integrals associated with a
#     Gauss-Markov process (GMP),
#  2) Computing mean, variance, and covariance for the GMP,
#  3) Simulating from bridges (conditional on a terminal value),
#  4) Sampling terminal values based on user-specified priors,
#  5) Solving an optimal stopping problem (OSP) via a basic
#     backward dynamic programming approach.
############################################################################

library(compiler)
library(truncnorm)

############################################################################
# 1. BUILDING GAUSS–MARKOV PARTIAL INTEGRALS
############################################################################
# These functions construct spline approximations to certain integrals:
#   - partial_slope(t) = ∫(0 to t) beta(u) du
#   - partial_trend(t) = ∫(0 to t) alpha(u) * exp(- partial_slope(u)) du
#   - partial_vol(t)   = ∫(0 to t) zeta^2(u) * exp(-2 * partial_slope(u)) du
# We use a fine time-grid and a trapezoidal rule, then store results in splinefun.
############################################################################

#-------------------------------------------------------------------
# gmp_integrals:
#   slope_fun  -> beta(t)
#   trend_fun  -> alpha(t)
#   vol_fun    -> zeta(t)
#   T          -> final time
#   n_grid     -> number of points for the trapezoid rule
#
# Returns a list with three spline functions for partial_slope, partial_trend,
# and partial_vol, plus the chosen final time T.
#-------------------------------------------------------------------
gmp_integrals <- function(slope_fun, trend_fun, vol_fun, T, n_grid = 5000) {
  t_grid <- seq(0, T, length.out = n_grid)
  dt     <- t_grid[2] - t_grid[1]
  
  # Evaluate user-provided functions on the time grid
  slope_vals <- slope_fun(t_grid)  # beta(t)
  trend_vals <- trend_fun(t_grid)  # alpha(t)
  vol_vals   <- vol_fun(t_grid)    # zeta(t)
  
  # partial_slope(t): integrates slope_fun from 0 to t
  partial_slope_vals <- numeric(n_grid)
  for (k in 2:n_grid) {
    partial_slope_vals[k] <- partial_slope_vals[k - 1] +
      0.5 * (slope_vals[k - 1] + slope_vals[k]) * dt
  }
  
  # partial_trend(t): ∫ alpha(u) * exp(- partial_slope(u)) du
  integrand_trend_vals <- trend_vals * exp(-partial_slope_vals)
  partial_trend_vals <- numeric(n_grid)
  for (k in 2:n_grid) {
    partial_trend_vals[k] <- partial_trend_vals[k - 1] +
      0.5 * (integrand_trend_vals[k - 1] + integrand_trend_vals[k]) * dt
  }
  
  # partial_vol(t): ∫ zeta^2(u) * exp(-2*partial_slope(u)) du
  integrand_vol_vals <- vol_vals^2 * exp(-2 * partial_slope_vals)
  partial_vol_vals <- numeric(n_grid)
  for (k in 2:n_grid) {
    partial_vol_vals[k] <- partial_vol_vals[k - 1] +
      0.5 * (integrand_vol_vals[k - 1] + integrand_vol_vals[k]) * dt
  }
  
  # Create spline approximations
  partial_slope_spline <- splinefun(t_grid, partial_slope_vals, method="natural")
  partial_trend_spline <- splinefun(t_grid, partial_trend_vals, method="natural")
  partial_vol_spline   <- splinefun(t_grid, partial_vol_vals,   method="natural")
  
  # Compile for speed and return as a list
  list(
    T                    = T,
    partial_slope_spline = cmpfun(partial_slope_spline),
    partial_trend_spline = cmpfun(partial_trend_spline),
    partial_vol_spline   = cmpfun(partial_vol_spline)
  )
}

############################################################################
# 2. GAUSS–MARKOV MEAN, VARIANCE, AND COVARIANCE
############################################################################
# Using the partial integrals, we define:
#   phi(t0, t1) = exp( partial_slope(t1) - partial_slope(t0) )
#
#   mean_gmp(t0, x, t1) = phi(t0,t1)* x
#                         + exp(partial_slope(t1)) *
#                           [ partial_trend(t1) - partial_trend(t0) ]
#
#   var_gmp(t0, t1) = exp(2*partial_slope(t1)) *
#                     [ partial_vol(t1) - partial_vol(t0) ]
#
#   cov_gmp(t1, t2) uses the factorization r1(min)*r2(max).
############################################################################

phi <- function(t0, t1, gmp) {
  ps <- gmp$partial_slope_spline
  exp(ps(t1) - ps(t0))
}

mean_gmp <- function(t0, x, t1, gmp) {
  ps   <- gmp$partial_slope_spline
  ptr  <- gmp$partial_trend_spline
  
  phi_val <- exp(ps(t1) - ps(t0))
  phi_val * x + exp(ps(t1)) * (ptr(t1) - ptr(t0))
}

var_gmp <- function(t0, t1, gmp) {
  ps  <- gmp$partial_slope_spline
  pvo <- gmp$partial_vol_spline
  exp(2 * ps(t1)) * (pvo(t1) - pvo(t0))
}

cov_gmp <- function(t1, t2, gmp) {
  # We define r1(t), r2(t) so that Cov(X_{t1}, X_{t2}) = r1(min)*r2(max)
  r1 <- function(t) {
    e_slope <- exp(gmp$partial_slope_spline(t))
    p_vol   <- gmp$partial_vol_spline(t)
    e_slope * p_vol
  }
  r2 <- function(t) {
    exp(gmp$partial_slope_spline(t))
  }
  
  t_min <- min(t1, t2)
  t_max <- max(t1, t2)
  r1(t_min) * r2(t_max)
}

############################################################################
# 3. SIMULATIONS
############################################################################
# simulate_gmp_bridge:
#    Simulate X(t1) given X(t0)=x, and a forced terminal value z = X(T).
#    This yields a normal conditional distribution:
#      Mean = m1 + [cov12 / vTsq] * (z - mT)
#      Var  = v1sq - (cov12^2 / vTsq)
#    where:
#      m1 = mean_gmp(t0, x, t1),
#      v1sq = var_gmp(t0, t1),
#      mT = mean_gmp(t0, x, T),
#      vTsq = var_gmp(t0, T),
#      cov12 = phi(t1, T)*v1sq.
############################################################################

simulate_gmp_bridge <- function(t0, t1, x, z, gmp) {
  # 1) unconditional from t0->t1
  m1    <- mean_gmp(t0, x, t1, gmp)
  v1sq  <- var_gmp(t0, t1, gmp)
  
  # 2) unconditional from t0->T
  mT    <- mean_gmp(t0, x, gmp$T, gmp)
  vTsq  <- var_gmp(t0, gmp$T, gmp)
  
  # 3) cov(X_{t1}, X_{T})
  cov12 <- phi(t1, gmp$T, gmp) * v1sq
  
  # 4) bridging formula
  condMean <- m1 + (cov12 / vTsq) * (z - mT)
  condVar  <- v1sq - (cov12^2 / vTsq)
  
  rnorm(length(z), mean = condMean, sd = sqrt(pmax(condVar, 1e-14)))
}
simulate_gmp_bridge <- cmpfun(simulate_gmp_bridge)

############################################################################
# 4. POSTERIOR SAMPLING OF z = X(T)
############################################################################
# sample_z_discrete:
#   For a discrete prior on z, we incorporate the likelihood factor
#   exp( - (z - mT)^2 / (2 * vTsq) ) * prior.
############################################################################

sample_z_discrete <- function(t0, x, gmp, support_z, probs_z, n = 1) {
  mT   <- mean_gmp(t0, x, gmp$T, gmp)
  vTsq <- var_gmp(t0, gmp$T, gmp)
  
  # Unnormalized posterior log-weights
  log_lik   <- - (support_z - mT)^2 / (2 * vTsq)
  log_prior <- log(probs_z)
  log_post  <- log_lik + log_prior
  
  # Normalize
  cmax <- max(log_post)
  w <- exp(log_post - cmax)
  w <- w / sum(w)
  
  sample(support_z, size = n, replace = TRUE, prob = w)
}
sample_z_discrete <- cmpfun(sample_z_discrete)

############################################################################
# sample_z_cont:
#   For a continuous prior, we either:
#    - Solve directly if prior is normal or truncated normal,
#    - Fall back to Metropolis-Hastings if user chooses "mh".
############################################################################

sample_z_cont <- function(
    t0, x, x0, gmp,
    prior_type = c("normal", "truncated_normal", "mh"),
    prior_params = list(mu = 0, sd = 1, lower = -Inf, upper = Inf),
    n_samples  = 1000,
    n_iter     = 5000,
    burn_in    = 1000,
    proposal_sd = 0.5,
    initial_val = NULL
) {
  prior_type <- match.arg(prior_type)
  
  m_t   <- mean_gmp(t0, x, gmp$T, gmp)
  v_tsq <- var_gmp(t0, gmp$T, gmp)
  
  # Bridge correction from 0->T:
  m_0   <- mean_gmp(0, x0, gmp$T, gmp)
  v_0sq <- var_gmp(0, gmp$T, gmp)
  
  # ---- Normal prior
  if (prior_type == "normal") {
    mu0 <- prior_params$mu
    sd0 <- prior_params$sd
    
    A <- (1 / v_tsq) - (1 / v_0sq) + (1 / (sd0^2))
    B <- (m_t / v_tsq) - (m_0 / v_0sq) + (mu0 / (sd0^2))
    if (A <= 0) stop("Posterior precision A <= 0; check parameters.")
    post_mu    <- B / A
    post_sigma <- sqrt(1 / A)
    return(rnorm(n_samples, mean = post_mu, sd = post_sigma))
  }
  
  # ---- Truncated normal prior
  if (prior_type == "truncated_normal") {
    mu0   <- prior_params$mu
    sd0   <- prior_params$sd
    lower <- prior_params$lower
    upper <- prior_params$upper
    
    A <- (1 / v_tsq) - (1 / v_0sq) + (1 / (sd0^2))
    B <- (m_t / v_tsq) - (m_0 / v_0sq) + (mu0 / (sd0^2))
    if (A <= 0) stop("Posterior precision A <= 0; check parameters.")
    post_mu    <- B / A
    post_sigma <- sqrt(1 / A)
    
    return(truncnorm::rtruncnorm(
      n_samples, a = lower, b = upper, mean = post_mu, sd = post_sigma
    ))
  }
  
  # ---- Metropolis-Hastings fallback
  if (is.null(initial_val)) {
    initial_val <- m_t
  }
  
  chain <- numeric(n_iter)
  chain[1] <- initial_val
  
  log_target <- function(z) {
    # Combined log target = [likelihood * prior * bridging correction] in log form
    term1 <- - (z - m_t)^2 / (2 * v_tsq)
    term2 <- log(prior_params$dens(z))  # user-provided prior density
    term3 <- - (z - m_0)^2 / (2 * v_0sq)
    term1 + term2 + term3
  }
  log_target <- cmpfun(log_target)
  
  Lcur <- log_target(initial_val)
  cur  <- initial_val
  
  for (i in 2:n_iter) {
    cand <- rnorm(1, mean = cur, sd = proposal_sd)
    cand_log_target <- log_target(cand)
    if (log(runif(1)) < (cand_log_target - Lcur)) {
      cur  <- cand
      Lcur <- cand_log_target
    }
    chain[i] <- cur
  }
  
  chain[(burn_in + 1):n_iter]
}
sample_z_cont <- cmpfun(sample_z_cont)

############################################################################
# 5. DYNAMIC PROGRAMMING APPROACH TO THE OSP
############################################################################
# osp_solver:
#   Uses a backward recursion over time (t_grid) and state (x_grid).
#   At each point (t, x), we compare:
#     - immediate payoff (here, identity = x)
#     - expected payoff of continuing (estimated via Monte Carlo).
#
#   The user can choose discrete vs continuous priors for X(T).
############################################################################

osp_solver <- function(
    t_grid = seq(0, 1, length.out = 20),
    x_grid = seq(-1, 3, length.out = 50),
    slope_fun = function(t) rep(0, length(t)),    # beta(t)
    trend_fun = function(t) rep(0, length(t)),    # alpha(t)
    vol_fun   = function(t) rep(1, length(t)),    # zeta(t)
    n_gmp_grid = 2000,
    x0 = 0,
    # Prior info:
    discrete_prior = TRUE,
    support = 0, probs = 1,   # For discrete prior
    prior_type = "normal",
    prior_params = list(mu = 0, sd = 0.5, lower = -Inf, upper = Inf),
    initial_val = NULL,
    M_samples = 2000,
    burn_in = 1000,
    proposal_sd = 1.0
) {
  # 1) Build partial integrals for the GMP
  gmp <- gmp_integrals(slope_fun, trend_fun, vol_fun,
                           T = max(t_grid), n_grid = n_gmp_grid)
  
  N_t <- length(t_grid)
  N_x <- length(x_grid)
  
  # 2) Initialize DP arrays
  V_matrix <- matrix(NA, nrow = N_t, ncol = N_x)
  decision_matrix <- matrix(NA, nrow = N_t, ncol = N_x)
  
  # 3) Terminal condition: immediate payoff = x at final time
  for (j in seq_len(N_x)) {
    V_matrix[N_t, j] <- x_grid[j]
    decision_matrix[N_t, j] <- 1  # Stop
  }
  
  # 4) Posterior sampler for the terminal value X(T)
  sample_z_func <- if (discrete_prior) {
    function(t0, x, M) {
      sample_z_discrete(t0, x, gmp, support, probs, n = M)
    }
  } else {
    function(t0, x, M) {
      sample_z_cont(t0, x, x0, gmp,
                    prior_type = prior_type,
                    prior_params = prior_params,
                    n_samples = M,
                    n_iter = (M + burn_in), burn_in = burn_in,
                    proposal_sd = proposal_sd,
                    initial_val = initial_val)
    }
  }
  
  # 5) Backward recursion over t_grid
  for (i in seq(N_t - 1, 1, by = -1)) {
    t0_val <- t_grid[i]
    t1_val <- t_grid[i + 1]
    
    # Build an interpolation of V_matrix at the next time step
    spline_V <- splinefun(x_grid, V_matrix[i + 1, ], method = "natural")
    cont_value <- function(x_val) spline_V(x_val)
    
    for (j in seq_len(N_x)) {
      x_cur <- x_grid[j]
      immediate_payoff <- x_cur  # Identity payoff
      
      # Sample possible X(T) from the posterior given X(t0_val)=x_cur
      z_samples <- sample_z_func(t0_val, x_cur, M_samples)
      
      # For each z_sample, simulate X(t1_val) consistent with X(T)=z
      x_next_vec <- simulate_gmp_bridge(t0_val, t1_val, x_cur, z_samples, gmp)
      
      # Estimate expected continuation at (t1_val, x_next)
      cont_est <- mean(cont_value(x_next_vec))
      
      # Compare immediate vs continuation
      if (immediate_payoff >= cont_est) {
        V_matrix[i, j] <- immediate_payoff
        decision_matrix[i, j] <- 1  # stop
      } else {
        V_matrix[i, j] <- cont_est
        decision_matrix[i, j] <- 0  # continue
      }
    }
    cat(sprintf("Time step i=%d, t=%.3f done.\n", i, t0_val))
  }
  
  list(
    t_grid          = t_grid,
    x_grid          = x_grid,
    V_matrix        = V_matrix,
    decision_matrix = decision_matrix
  )
}
osp_solver <- compiler::cmpfun(osp_solver)


############################################################################
# EXAMPLE USAGE
############################################################################
# slope_fun  <- function(t) rep(-5, length(t))
# trend_fun  <- function(t) rep(5,  length(t))
# vol_fun    <- function(t) rep(1,  length(t))
# t_grid     <- log(seq(1, exp(1), length.out = 1000))
# x_grid     <- seq(-20, 20, length.out = 1000)
#
# sol <- osp_solver(
#   t_grid = t_grid,
#   x_grid = x_grid,
#   slope_fun = slope_fun,
#   trend_fun = trend_fun,
#   vol_fun   = vol_fun,
#   n_gmp_grid = 50000,
#   x0 = 0,
#   discrete_prior = FALSE,
#   prior_type     = "normal",
#   prior_params   = list(mu=0, sd=0.75, lower=-Inf, upper=Inf),
#   M_samples      = 10000,
#   burn_in        = 10,
#   proposal_sd    = 1.0
# )
#
# image(t_grid, x_grid, sol$decision_matrix,
#       col = c("darkolivegreen3","firebrick3"),
#       xlab = "Time", ylab = "State (x)",
#       main = "Stopping(1) vs Continuing(0)")
# lines(c(0,1), c(0,0), lty = 2)
