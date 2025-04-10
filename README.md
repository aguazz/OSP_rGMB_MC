# Optimal Stopping for Random Gauss–Markov Bridges

This repository contains R scripts that implement and test numerical methods for solving an Optimal Stopping Problem (OSP) for Gauss–Markov processes (GMPs). In particular, the approach involves transforming the underlying process, simulating bridges, and then solving an OSP by backward dynamic programming. The repository provides the following main scripts:

- **functions.R**: Contains core functions for creating Gauss–Markov process (GMP) splines, computing mean/variance/covariance, simulating bridges, and solving the optimal stopping problem.
- **data_generator.R**: Uses the functions defined in *functions.R* to generate simulation data under several parameter regimes (cases) and saves the OSP solution objects in the `data/` directory.
- **image_generator.R**: Loads the data files produced by *data_generator.R* and generates decision region plots. A parameter allows toggling between a green/red color scheme and a white/gray scheme.
- **test.R**: Contains unit tests (via the `testthat` package) and visual checks to verify that the implemented functions behave as expected. This includes comparing the numerical mean/variance/covariance produced by the splines against closed-form values for an Ornstein–Uhlenbeck (OU) process and testing the performance of the Metropolis–Hastings (MH) samplers.

## Repository Structure

```
.
├── data_generator.R       # Script to generate simulation data (OSP solutions)
├── functions.R            # Core functions for the GMP, simulation, and OSP solver
├── image_generator.R      # Script for generating plots (decision regions) from simulation results
├── test.R                 # Unit tests and visual checks for the core functions
├── data/                  # (Automatically created) Directory where simulation data (.RData) is saved
├── img/                   # (Automatically created) Directory where image files (PNG) are saved
└── README.md              # This file
```

## Prerequisites

Make sure you have R installed, along with the following packages:

- `compiler`
- `truncnorm`
- `scales`
- `tools`
- `testthat`

You can install missing packages using:

```r
install.packages(c("compiler", "truncnorm", "scales", "tools", "testthat"))
```

## Scripts Overview

### 1. functions.R

This is the main script containing functions which:

- **build_gmp_splines(...)**: Computes partial integrals (for slope, trend, and volatility) on a fine time grid using the trapezoidal rule and returns spline functions.
- **phi(...)**, **mean_gmp(...)**, **var_gmp(...)**, **cov_gmp(...)**: Compute the transition kernel components (exponential factor, mean, variance, covariance) for the GMP.
- **simulate_gmp_bridge(...)**: Simulates a bridge from time `t0` to `t1` given an initial condition and a forced terminal value using the numerical approximations.
- **sample_z_discrete(...)** and **sample_z_cont(...)**: Sample terminal values from the posterior distribution of the process at time `T` given the current state; supports both discrete and continuous (including MH-based) sampling.
- **osp_solver(...)**: Implements a backward dynamic programming approach to solve the OSP in the original `(t, x)` coordinates. It uses the above functions along with user-defined grids and prior information.

Each function has in-line comments to explain the main steps and assumptions. Functions are optimized using `compiler::cmpfun` where appropriate.

### 2. data_generator.R

This script generates simulation data under multiple scenarios:
- **Case 1**: Two-points Brownian Bridge (discrete prior with terminal values \{-1, 1\}).
- **Case 2**: Sinusoidal Mean Reversion.
- **Case 3**: Rapid Switching of Pulling Strength.
- **Case 4**: Parabolic Volatility Smile.
- **Case 5**: Greater Final Variance.

For each case, a set of parameters is defined and the OSP is solved by calling `osp_solver()`. The resulting solutions (which include value matrices and decision regions) are saved into the `data/` folder with descriptive filenames.

### 3. image_generator.R

This script extracts the saved OSP solution data and plots the decision regions. It includes the following features:
- A Boolean parameter `use_green_red` that toggles the color scheme between green/red (when `TRUE`) and white/gray (when `FALSE`).
- Images are saved as PNG files in the `img/` directory if `save_img` is set to `TRUE`; otherwise, they are displayed on-screen.
- Only decision region images are generated (plots of the value function have been removed).

### 4. test.R

This script uses the `testthat` package to perform automated tests and visual checks:
- **Test 1**: Compares the numerical mean, variance, and covariance computed using the GMP splines with closed-form expressions for an OU process under random parameter choices.
- **Test 2**: Generates plots comparing the spline-based results with the closed-form solutions for visual verification.
- **Tests 3 & 4**: Evaluate the performance of the Metropolis–Hastings sampler for both a standard normal prior and a truncated normal prior by comparing the MCMC sample density with the exact posterior density.

## How to Run

From your R console or RStudio, run the whole data generator script (`source("data_generator.R")`) or just isolated cases within. Run the image generator script (`source("image_generator.R")`) to obtain the associated images. You can also create/add your own cases following the template in these two scripts.
