#########################################################
# Simulation Study: Misspecified F in CPCM Graph Estimation
#
# This script runs a grid of simulation settings:
#   - Different families of distributions
#   - Different functional forms for θ(X1)
#
# For each combination, we:
#   1. Generate n samples
#   2. Fit CPCM_graph_estimate
#   3. Check if it recovers the correct causal edge (1 → 2)
#
#
# Author: Juraj Bodik
#########################################################
library(MASS)   # for mvrnorm
library(dplyr)  # for data manipulation
library(purrr)  # for functional mapping
library(tidyr)  # for reshaping results

source('Main_function.R')

set.seed(1)
n = 1000
reps = 100
#########################################################
# Helper functions to generate random smooth θ(x)
#########################################################

# Draws a smooth function from a squared-exponential GP
# Used for random θ(x) in the simulation
random_function <- function(draw_function = FALSE, max = 10){
  cov_matrix <- function(x, kernel_fn, ...) {
    outer(x, x, function(a, b) kernel_fn(a, b, ...))
  }
  draw_samples <- function(x, N, kernel_fn, ...) {
    Y <- matrix(NA, nrow = length(x), ncol = N)
    for (n in 1:N) {
      K <- cov_matrix(x, kernel_fn, ...)
      Y[, n] <- mvrnorm(1, mu = rep(0, length(x)), Sigma = K)
    }
    Y
  }
  se_kernel <- function(x, y, sigma = 1, length = 1) {
    sigma^2 * exp(- (x - y)^2 / (2 * length^2))
  }
  
  sekvencia <- seq(0, max, length.out = 201)
  Y <- draw_samples(sekvencia, 1, kernel_fn = se_kernel, length = 0.2)
  Y <- Y - min(Y) + 0.5  # shift so θ > 0
  return(Y)
}

# Matches each x to the closest value in the GP draw and returns θ(x)
random_theta <- function(x, f, max){
  sekvencia <- seq(0, max, length.out = 201)
  theta <- numeric(length(x))
  for (i in seq_along(x)) {
    theta[i] <- f[which.min(abs(x[i] - sekvencia))]
  }
  theta
}

#########################################################
# Functional forms for θ(x)
#########################################################

theta_fns <- list(
  "theta(x)=x"         = function(x) x,
  "theta(x)=x^2+1"     = function(x) x^2 + 1,
  "theta(x)=exp(x)/2"  = function(x) exp(x) / 2,
  "theta(x)=random GP" = function(x) {
    f <- random_function(max = ceiling(max(x)))
    random_theta(x, f, max = ceiling(max(x)))
  }
)



#########################################################
# Families of distributions to test
#########################################################

families <- c(
  "Gamma with fixed scale",
  "Gamma",
  "Pareto",
  "Gumbel",
  "Gumbel with fixed scale",
  "Gaussian with fixed sigma",
  "Gaussian", 
  "Exponential"
)



#########################################################
# Safe per-cell runner
#########################################################
run_one_cell <- function(family, theta_label, n = 100, reps = 100) {

    successes <- 0L

    # Compute θ function once
    theta_fn <- theta_fns[[theta_label]]
    if (is.null(theta_fn)) stop("Unknown theta function label: ", theta_label)
    
    for (k in 1:reps) {
      tryCatch({
        # Generate X1 ~ N(2, 1) truncated to positive values
        X1 <- c()
        while (length(X1) < n) {
          value <- rnorm(1, 2, 1)
          if (value > 0) X1 <- c(X1, value)
        }
        
        # θ(X1)
        theta <- theta_fn(X1)
        
        # X2 | X1 ~ Exp(rate = θ)
        X2 <- rexp(n, rate = theta)
        X <- data.frame(X1 = X1, X2 = X2)
        if (family %in% c("Pareto", "Pareto2")) X <- X + 1
        
        graph <- CPCM_graph_estimate(X, family_of_distributions = family)
        graph
        successes <- successes + 2 - as.numeric(substring(graph[[1]][6], 1, 1))
      },
      error = function(e) {
        # On error, count as 0 successes for this replicate
        successes <- successes + 0
      })
    }
    
    acc <- successes / reps
    
   
  
    tibble(
      family      = family,
      theta       = theta_label,
      accuracy    = round(acc, 3)
    )

}

#########################################################
# Full grid runner: survives bad cells + checkpoints
#########################################################


checkpoint_file = "sim3_checkpoint_long.csv"

  grid <- expand.grid(family = families, theta = names(theta_fns), stringsAsFactors = FALSE)
  total_cells <- nrow(grid)
  results <- vector("list", total_cells)
  start_time <- Sys.time()
  
  for (i in seq_len(total_cells)) {
    cell_start <- Sys.time()
    fam <- grid$family[i]
    th  <- grid$theta[i]
    
    cat(sprintf("\n[%d/%d] Running family='%s', theta='%s'...\n", 
                i, total_cells, fam, th))
    
    # Run safely
    res <- run_one_cell(fam, th, n = n, reps = reps)
    results[[i]] <- res
    
    # Progress + ETA
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    avg_time_per_cell <- elapsed / i
    remaining_time <- avg_time_per_cell * (total_cells - i)
    cell_time <- as.numeric(difftime(Sys.time(), cell_start, units = "secs"))
    cat(sprintf("   Cell done in %.1f sec. Estimated time left: ~%.1f min\n",
                cell_time, remaining_time / 60))
    
    # Checkpoint after each cell so nothing is lost
    suppressWarnings({
      tmp_df <- bind_rows(results[!vapply(results, is.null, logical(1))])
      if (nrow(tmp_df) > 0) {
        write.csv(tmp_df, checkpoint_file, row.names = FALSE)
      }
    })
  }
  
results_df <- bind_rows(results)
  
all_results = results_df %>%
      select(family, theta, accuracy) %>%
      tidyr::pivot_wider(names_from = theta, values_from = accuracy) %>%
      arrange(factor(family, levels = families))
  



print(all_results)

# write.csv(all_results$long, "sim_misspecF_long.csv", row.names = FALSE)
# write.csv(all_results$wide, "sim_misspecF_table.csv", row.names = FALSE)

