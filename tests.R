############################################################################
# test.R
#
# This script uses 'testthat' to validate the correctness of core functions
# from functions.R, focusing on:
#   1) Comparison of Gauss???Markov mean/variance/covariance to a closed-form
#      Ornstein???Uhlenbeck (OU) model,
#   2) Visual checks of Metropolis???Hastings samplers for normal and truncated
#      normal priors.
############################################################################

library(testthat)
source("functions.R")  # Ensure this points to your main GMP functions file

# -----------------------------------------------------
# Closed-form for an OU SDE:
#   dX_t = (a + b X_t) dt + c dB_t
# with times t0 < t1 and X_{t0} = x0.
# This is used to check the numerical results of the spline-based approach.
# -----------------------------------------------------
ou_mean_var_cov_closed_form <- function(t0, x0, t1, a, b, c) {
  dt <- t1 - t0
  if (dt < 0) stop("t1 < t0 not allowed")
  
  # Mean
  if (abs(b) < 1e-14) {
    mval <- x0 + a * dt
  } else {
    e_bdt <- exp(b * dt)
    mval  <- e_bdt * x0 + (a / b) * (e_bdt - 1)
  }
  
  # Variance
  if (abs(b) < 1e-14) {
    vval <- c^2 * dt
  } else {
    e_2bdt <- exp(2 * b * dt)
    vval   <- (c^2 / (2 * b)) * (e_2bdt - 1)
  }
  
  # Covariance( X_{t0}, X_{t1} ) under an OU
  if (abs(b) < 1e-14) {
    cov_val <- c^2 * t0
  } else {
    cov_val <- (c^2 / (2 * b)) * (exp(2 * b * t0) - 1) * exp(b * (t1 - t0))
  }
  
  list(mean = mval, var = vval, cov = cov_val)
}


# -----------------------------------------------------
# TEST 1: Compare GMP-based results to closed-form OU
# -----------------------------------------------------
test_that("Ornstein???Uhlenbeck (OU) closed-form matches spline-based GMP", {
  
  set.seed(123)  # reproducibility
  
  n_runs <- 50
  T_max  <- 1   # maximum horizon
  
  for (i in seq_len(n_runs)) {
    # Draw random parameters for a, b, c, times, and initial condition
    a <- runif(1, min = -5, max = 1)
    b <- runif(1, min = -5, max = 1)
    c <- runif(1, min = 0.1, max = 2)
    t0 <- runif(1, min = 0, max = T_max - 0.5)
    x0 <- runif(1, min = -5, max = 5)
    
    # Build GMP splines with constant slope/trend/vol
    slope_fun <- function(t) rep(b, length(t))
    trend_fun <- function(t) rep(a, length(t))
    vol_fun   <- function(t) rep(c, length(t))
    
    gmp_ou <- gmp_integrals(slope_fun, trend_fun, vol_fun,
                                T = T_max, n_grid = 5000)
    
    # Evaluate from t0 to T_max in 100 steps
    t1_seq <- seq(t0, T_max, length.out = 100)
    
    # Numeric (GMP) approximations
    m_approx  <- sapply(t1_seq, function(t1) mean_gmp(t0, x0, t1, gmp_ou))
    v_approx  <- sapply(t1_seq, function(t1) var_gmp(t0, t1, gmp_ou))
    cv_approx <- sapply(t1_seq, function(t1) cov_gmp(t0, t1, gmp_ou))
    
    # Closed-form OU
    exact_vals <- sapply(t1_seq, function(t1) {
      res <- ou_mean_var_cov_closed_form(t0, x0, t1, a, b, c)
      c(mean = res$mean, var = res$var, cov = res$cov)
    })
    m_exact  <- exact_vals["mean", ]
    v_exact  <- exact_vals["var", ]
    cv_exact <- exact_vals["cov", ]
    
    # Tolerance check
    tol <- 1e-4
    expect_true(max(abs(m_approx  - m_exact))  < tol, "Mean mismatch too large.")
    expect_true(max(abs(v_approx  - v_exact))  < tol, "Variance mismatch too large.")
    expect_true(max(abs(cv_approx - cv_exact)) < tol, "Covariance mismatch too large.")
  }
})


# -----------------------------------------------------
# TEST 2: Visual tests for random OU parameters
# -----------------------------------------------------
# This section generates quick plots to visually compare mean, var, cov from
# the GMP splines vs. closed-form. Not an automated test; use testthat mostly
# for convenience of organization.

# Random parameters
set.seed(456)
a <- runif(1, min = -5, max = 2)
b <- runif(1, min = -5, max = 2)
c <- runif(1, min = 0.1, max = 2)
T_max  <- 1
t0 <- runif(1, min = 0, max = T_max - 0.5)
x0 <- runif(1, min = -5, max = 5)

slope_fun <- function(t) rep(b, length(t))
trend_fun <- function(t) rep(a, length(t))
vol_fun   <- function(t) rep(c, length(t))

gmp_ou <- gmp_integrals(slope_fun, trend_fun, vol_fun, T = T_max, n_grid = 5000)
t1_seq <- seq(t0, T_max, length.out = 100)

# GMP-based
m_approx  <- sapply(t1_seq, function(t1) mean_gmp(t0, x0, t1, gmp_ou))
v_approx  <- sapply(t1_seq, function(t1) var_gmp(t0, t1, gmp_ou))
cv_approx <- sapply(t1_seq, function(t1) cov_gmp(t0, t1, gmp_ou))

# Closed-form OU
exact_vals <- sapply(t1_seq, function(t1) {
  res <- ou_mean_var_cov_closed_form(t0, x0, t1, a, b, c)
  c(mean = res$mean, var = res$var, cov = res$cov)
})
m_exact  <- exact_vals["mean", ]
v_exact  <- exact_vals["var", ]
cv_exact <- exact_vals["cov", ]

# Plots
op <- par(no.readonly = TRUE)
par(mfrow = c(1,3))

# Mean
plot(t1_seq, m_exact, type = "l", col = "blue", lwd = 2,
     xlab = "t1", ylab = "Mean", main = "Mean Comparison")
lines(t1_seq, m_approx, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Closed Form", "GMP"), col = c("blue","red"),
       lty = c(1,2), lwd = 2)

# Variance
plot(t1_seq, v_exact, type = "l", col = "blue", lwd = 2,
     xlab = "t1", ylab = "Variance", main = "Variance Comparison")
lines(t1_seq, v_approx, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Closed Form", "GMP"), col = c("blue","red"),
       lty = c(1,2), lwd = 2)

# Covariance
plot(t1_seq, cv_exact, type = "l", col = "blue", lwd = 2,
     xlab = "t1", ylab = "Covariance", main = "Covariance Comparison")
lines(t1_seq, cv_approx, col = "red", lwd = 2, lty = 2)
legend("topleft", legend = c("Closed Form", "GMP"), col = c("blue","red"),
       lty = c(1,2), lwd = 2)

par(op)


########################################
# TEST 3: Metropolis???Hastings for Normal Prior
########################################
# Here we check that the MCMC chain matches the known normal posterior.

# 1) Define the prior
sigma0 <- 0.5 # standard deviation
prior_type <- "normal"

# 2) OU model parameters
a <- -5
b <- -5
c <- 1
T_val <- 1
x0 <- 0
t0 <- 0.1
x  <- 1

# 3) Build GMP partial integrals (constant coefficients)
slope_fun <- function(t) rep(b, length(t))
trend_fun <- function(t) rep(a, length(t))
vol_fun   <- function(t) rep(c, length(t))

gmp_ou <- gmp_integrals(slope_fun, trend_fun, vol_fun, T = T_val, n_grid = 10000)

# 4) MCMC
chain <- sample_z_cont(t0, x, x0, gmp_ou,
                       prior_type = prior_type,
                       prior_params = list(mu = 0, sd = sigma0, lower = -Inf, upper = Inf),
                       n_samples  = 5000, proposal_sd = 0.2)

# 5) Exact posterior parameters for a bridging-based normal posterior
m_t   <- mean_gmp(t0, x, gmp_ou$T, gmp_ou)
v_tsq <- var_gmp(t0, gmp_ou$T, gmp_ou)
m_0   <- mean_gmp(0, x0, gmp_ou$T, gmp_ou)
v_0sq <- var_gmp(0, gmp_ou$T, gmp_ou)

A <- (1 / v_tsq) - (1 / v_0sq) + (1 / sigma0^2)
B <- (m_t / v_tsq) - (m_0 / v_0sq)
post_mu <- B / A
post_sigma <- sqrt(1 / A)

# 6) Visual check
par(mfrow = c(1,1))
hist(chain, breaks = 50, probability = TRUE,
     main = "Metropolis???Hastings vs. Exact Normal Posterior",
     xlab = "z", col = rgb(0,0,1,0.3), border = "white")
curve(dnorm(x, mean = post_mu, sd = post_sigma), add = TRUE,
      col = "red", lwd = 2)
legend("topright", legend = c("Exact Posterior", "MCMC samples"),
       col = c("red", rgb(0,0,1,0.3)), lwd = c(2, NA), pch = c(NA, 15))


#-------------------------------------------------------------------
# TEST 4: Metropolis???Hastings for a Truncated Normal Prior
#-------------------------------------------------------------------

# 1) Define the prior
sigma0 <- 0.5 # standard deviation
prior_type <- "truncated_normal"

# OU parameters
a <- 0
b <- -1
c <- 1
T_val <- 1
x0 <- 0
t0 <- 0.3
x  <- 0

# 3) Build GMP partial integrals (constant coefficients)
slope_fun <- function(t) rep(b, length(t))
trend_fun <- function(t) rep(a, length(t))
vol_fun   <- function(t) rep(c, length(t))

gmp_ou <- gmp_integrals(slope_fun, trend_fun, vol_fun, T = T_val, n_grid = 5000)

# 4) MCMC
chain <- sample_z_cont(t0, x, x0, gmp_ou,
                       prior_type = prior_type,
                       prior_params = list(mu = 0, sd = sigma0, lower = 0, upper = Inf),
                       n_samples  = 5000, proposal_sd = 0.2)

# The posterior is a truncated normal. We compute parameters, then scale by
# the truncation factor at 0.
m_t   <- mean_gmp(t0, x, gmp_ou$T, gmp_ou)
v_tsq <- var_gmp(t0, gmp_ou$T, gmp_ou)
m_0   <- mean_gmp(0, x0, gmp_ou$T, gmp_ou)
v_0sq <- var_gmp(0, gmp_ou$T, gmp_ou)

A <- (1 / v_tsq) - (1 / v_0sq) + (1 / sigma0^2)
B <- (m_t / v_tsq) - (m_0 / v_0sq)
post_mu <- B / A
post_sigma <- sqrt(1 / A)

# Exact truncated normal posterior
exact_post_density <- Vectorize(function(z) {
  if (z < 0) {
    0
  } else {
    # Normal density / normalizing factor for z >= 0
    denom <- 1 - pnorm(0, mean = post_mu, sd = post_sigma)
    dnorm(z, mean = post_mu, sd = post_sigma) / denom
  }
})

par(mfrow = c(1,1))
hist(chain, breaks = 50, probability = TRUE,
     main = "MCMC Samples vs. Exact Posterior (Truncated)",
     xlab = "z", col = rgb(0,0,1,0.3), border = "white")
curve(exact_post_density(x), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = c("Exact Posterior", "MCMC Density"),
       col = c("red", rgb(0,0,1,0.3)), lty = c(1,NA), pch = c(NA,15))

