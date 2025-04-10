##############################
# image_generator.R
#
# This script generates plots (decision regions) from the OSP solutions saved as .RData files.
# It includes a boolean parameter 'use_green_red' to switch color schemes:
#   - TRUE  => green/red
#   - FALSE => white/gray
##############################

# 1) Set to TRUE to save images (PNG) or FALSE to display on-screen only.
save_img <- FALSE 

# 2) Color scheme parameter:
#    Set this to TRUE for green/red or FALSE for white/gray.
use_green_red <- FALSE 

# 3) Ensure output directory for plots
if(!dir.exists("img")) dir.create("img")

# Libraries
library(scales)  # alpha(), etc.
library(tools)

# Helper function to pick colors based on use_green_red
decide_colors <- function(use_gr) {
  if (use_gr) {
    # Original scheme
    c("darkolivegreen3", "firebrick3")
  } else {
    # White/gray fallback
    c("white", "grey70")
  }
}
chosen_colors <- decide_colors(use_green_red)

##############################
# CASE 1: TWO-POINTS BROWNIAN BRIDGE
##############################
files <- list.files(
  "data", 
  pattern = "^[0-9]{3}_osp_solution_2pBB_x0_(-?\\d+\\.\\d+)\\.RData$", 
  full.names = TRUE
)
pattern <- "^[0-9]{3}_osp_solution_2pBB_x0_(-?\\d+\\.\\d+)\\.RData$"

for (file in files) {
  load(file)  # loads 'osp_solution'
  base <- file_path_sans_ext(basename(file))
  
  # Extract x0 from filename
  m <- regexec(pattern, basename(file))
  reg_matches <- regmatches(basename(file), m)
  if (length(reg_matches[[1]]) >= 2) {
    x0_val <- reg_matches[[1]][2]
  } else {
    x0_val <- "NA"
  }
  
  # Extract variables
  decision_matrix <- osp_solution$decision_matrix
  t_grid          <- osp_solution$t_grid
  x_grid          <- osp_solution$x_grid
  
  # Plot Decision Regions
  if(save_img){
    png(file.path("img", paste0(base, "_decision_regions.png")), width = 800, height = 700)
  }
  par(mar = c(3, 4, 2, 4), mgp = c(3, 2, 0), cex.axis = 2)
  image(t_grid, x_grid, decision_matrix,
        col = chosen_colors,
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = seq(min(t_grid), max(t_grid), length.out = 5))
  axis(2)
  if(save_img) dev.off()
  
  cat("Generated images for", base, "\n")
}


##############################
# CASE 2: SINUSOIDAL MEAN REVERSION
##############################
files <- list.files(
  "data",
  pattern = "^osp_solution_sin_mean_rev_x0_[+-]?[0-9]+\\.[0-9]{2}\\.RData$",
  full.names = TRUE
)
for (file in files) {
  load(file)
  base <- file_path_sans_ext(basename(file))
  
  m <- regexec("x0_([+-]?\\d+\\.\\d+)", basename(file))
  reg_matches <- regmatches(basename(file), m)
  if (length(reg_matches[[1]]) >= 2) {
    x0_val <- reg_matches[[1]][2]
  } else {
    x0_val <- "NA"
  }
  
  decision_matrix <- osp_solution$decision_matrix
  t_grid          <- osp_solution$t_grid
  x_grid          <- osp_solution$x_grid
  
  # Decision Regions
  if(save_img){
    png(file.path("img", paste0(base, "_decision_regions.png")), width = 800, height = 700)
  }
  par(mar = c(3, 4, 2, 4), mgp = c(3, 2, 0), cex.axis = 2)
  image(t_grid, x_grid, decision_matrix,
        col = chosen_colors,
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = seq(min(t_grid), max(t_grid), length.out = 5))
  axis(2)
  lines(c(0,1), c(0,0), lty = 2)
  if(save_img) dev.off()
  
  cat("Generated images for", base, "\n")
}


##############################
# CASE 3: EXPONENTIAL PULLING TOWARD ZERO
##############################
files <- list.files(
  "data",
  pattern = "^osp_solution_slope_switch_x0_[+-]?[0-9]+\\.[0-9]{2}\\.RData$",
  full.names = TRUE
)
for (file in files) {
  load(file)
  base <- file_path_sans_ext(basename(file))
  
  m <- regexec("x0_([+-]?\\d+\\.\\d+)", basename(file))
  reg_matches <- regmatches(basename(file), m)
  if (length(reg_matches[[1]]) >= 2) {
    x0_val <- reg_matches[[1]][2]
  } else {
    x0_val <- "NA"
  }
  
  decision_matrix <- osp_solution$decision_matrix
  t_grid          <- osp_solution$t_grid
  x_grid          <- osp_solution$x_grid
  
  # Decision Regions
  if(save_img){
    png(file.path("img", paste0(base, "_decision_regions.png")), width = 800, height = 700)
  }
  par(mar = c(3, 4, 2, 4), mgp = c(3, 2, 0), cex.axis = 2)
  image(t_grid, x_grid, decision_matrix, ylim = c(-0.2, 1),
        col = chosen_colors,
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = seq(min(t_grid), max(t_grid), length.out = 5))
  axis(2)
  lines(c(0,1), c(0,0), lty = 2)
  if(save_img) dev.off()
  
  cat("Generated images for", base, "\n")
}


##############################
# CASE 4: PARABOLIC VOLATILITY SMILE
##############################
files <- list.files(
  "data",
  pattern = "^osp_solution_vol_smile_x0_[+-]?[0-9]+\\.[0-9]{2}\\.RData$",
  full.names = TRUE
)
for (file in files) {
  load(file)
  base <- file_path_sans_ext(basename(file))
  
  m <- regexec("x0_([+-]?\\d+\\.\\d+)", basename(file))
  reg_matches <- regmatches(basename(file), m)
  if (length(reg_matches[[1]]) >= 2) {
    x0_val <- reg_matches[[1]][2]
  } else {
    x0_val <- "NA"
  }
  
  decision_matrix <- osp_solution$decision_matrix
  t_grid          <- osp_solution$t_grid
  x_grid          <- osp_solution$x_grid
  
  # Decision Regions
  if(save_img){
    png(file.path("img", paste0(base, "_decision_regions.png")), width = 800, height = 700)
  }
  par(mar = c(3, 4, 2, 4), mgp = c(3, 2, 0), cex.axis = 2)
  image(t_grid, x_grid, decision_matrix, ylim = c(-0.2, 8),
        col = chosen_colors,
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = seq(min(t_grid), max(t_grid), length.out = 5))
  axis(2)
  lines(c(0,1), c(0,0), lty = 2)
  if(save_img) dev.off()
  
  cat("Generated images for", base, "\n")
}


##############################
# CASE 5: GREATER FINAL VARIANCE
##############################
files <- list.files(
  "data",
  pattern = "^osp_solution_greater_final_var_x0_[+-]?[0-9]+\\.[0-9]{2}\\.RData$",
  full.names = TRUE
)
for (file in files) {
  load(file)
  base <- file_path_sans_ext(basename(file))
  
  m <- regexec("x0_([+-]?\\d+\\.\\d+)", basename(file))
  reg_matches <- regmatches(basename(file), m)
  if (length(reg_matches[[1]]) >= 2) {
    x0_val <- reg_matches[[1]][2]
  } else {
    x0_val <- "NA"
  }
  
  decision_matrix <- osp_solution$decision_matrix
  t_grid          <- osp_solution$t_grid
  x_grid          <- osp_solution$x_grid
  
  # Decision Regions
  if(save_img){
    png(file.path("img", paste0(base, "_decision_regions.png")), width = 800, height = 700)
  }
  par(mar = c(3, 4, 2, 4), mgp = c(3, 2, 0), cex.axis = 2)
  image(t_grid, x_grid, decision_matrix, ylim = c(-8, 8),
        col = chosen_colors,
        xlab = "", ylab = "", axes = FALSE)
  axis(1, at = seq(min(t_grid), max(t_grid), length.out = 5))
  axis(2)
  lines(c(0,1), c(0,0), lty = 2)
  if(save_img) dev.off()
  
  cat("Generated images for", base, "\n")
}

