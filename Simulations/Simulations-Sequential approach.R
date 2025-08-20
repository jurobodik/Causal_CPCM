###########################################################################
# Sequential approach simulations
#
# This script demonstrates the performance of the Sequential approach 
# compared to the Oracle approach.
#
# Author: Juraj Bodik
# Reference:  https://arxiv.org/abs/2303.15376
###########################################################################
library(splines)
library(dplyr)    # for experiment loop
library(ggplot2)
library(tidyr)
library(patchwork)  # for side-by-side plots
source('CPCM_function.R')

generate_data <- function(n, family_of_distributions = 1, dist_choice = NULL) {
  
  repeat {
    X <- rnorm(n)
    basis <- ns(X, df = 3)
    theta <- as.vector(basis %*% runif(ncol(basis), -2, 2))  # smooth function
    
    # Clip theta to avoid extreme or tiny values
    theta <- 0.5 + 4.5 * (tanh(theta) + 1) / 2  # maps smoothly to [0.5, 5]    
    # Choose a distribution (random if not provided)
    if (is.null(dist_choice)) {
      dist_choices_1 <- c("Gaussian with fixed sigma", "Poisson", "Pareto", "Exponential")
      dist_choices_2 <- c("Gaussian", "Negative_binomial", "Pareto2", "Gamma")
      dist_choice <- if (family_of_distributions == 1) {
        sample(dist_choices_1, 1)
      } else {
        sample(dist_choices_2, 1)
      }
    }
    
    # Try to sample Y, with error handling
    try_result <- try({
      Y <- switch(dist_choice,
                  "Gaussian with fixed sigma" = rnorm(n, mean = theta, sd = 1),
                  "Poisson" = rpois(n, lambda = pmax((theta), 1e-3)),
                  "Pareto" = 1 / runif(n)^(1 / theta),
                  "Exponential" = rexp(n, rate = pmax((theta), 1e-3)),
                  
                  "Gaussian" = rnorm(n, mean = theta, sd = abs(X)+1),
                  "Negative_binomial" = rnbinom(n, mu = pmax((theta), 1e-3), size = abs(X)+1),
                  "Pareto2" = exp(rgamma(n, shape = theta, rate = abs(X)+1)),
                  "Gamma" = rgamma(n, shape = theta, rate = abs(X)+1),
                  
                  stop("Invalid distribution choice.")
      )
    }, silent = TRUE)
    
    # Check for errors or non-finite values
    if (inherits(try_result, "try-error")) next
    if (any(!is.finite(X)) || any(!is.finite(theta)) || any(!is.finite(Y))) next
    
    # Return clean data
    return(list(data = data.frame(X = X, Y = Y),
                X = X, Y = Y, theta = theta,
                distribution_used = dist_choice))
  }
}

########################### Example of usage ########################


example <- generate_data(n = 500, family_of_distributions = 2)
estimate = CPCM_graph_estimate(example$data, family_of_distributions = 'Sequential choice')

sub(".*;", "", estimate[[1]][6])=="1 --> 2" # Is the graph estimated correctly '1 --> 2'?
sub(".*;", "", estimate[[1]][7])==example$distribution_used   #Is the chosen distribution correct?

####################### Repeat experiment and plot it ##########################
set.seed(1)  # for reproducibility
ns <- c(100, 250, 500, 1000, 1500, 2000)
families <- c(1, 2)
reps <- 100

results <- expand.grid(n = ns, fam = families) %>%
  rowwise() %>%
  mutate(
    correct_pct = list({
      count_correct <- 0
      count_graph <- 0
      count_oracle_graph <- 0
      
      for (i in 1:reps) {
        if(i%%10==0) cat(sprintf("Running: n = %d | family = %d | repetition = %d\n", n, fam, i))
        
        example <- generate_data(n = n, family_of_distributions = fam)
        
        estimate <- try(CPCM_graph_estimate(example$data, family_of_distributions = 'Sequential choice'), silent = TRUE)
        if (inherits(estimate, "try-error")) next
        
        estimated_dist <- sub(".*;", "", estimate[[1]][7])
        true_dist <- example$distribution_used
        
        if (estimated_dist == true_dist) {
          count_correct <- count_correct + 1
        }
        if (estimate[[1]][6] == "1 --> 2") {
          count_graph <- count_graph + 1
        }
        
        if (estimated_dist == true_dist && estimate[[1]][6] == "1 --> 2") {
          count_oracle_graph <- count_oracle_graph + 1
        }
        
        if(estimated_dist != true_dist){
          estimate_oracle <- try(CPCM_graph_estimate(example$data, family_of_distributions = fam), silent = TRUE)
          if (inherits(estimate_oracle, "try-error")) next
          
          if(estimate_oracle[[1]][6] == "1 --> 2"){  count_oracle_graph <- count_oracle_graph + 1  }
        }
        
        if (i %% 50 == 0) {
          cat(sprintf("  Intermediate: dist = %.1f%%, graph = %.1f%% (%d/%d)\n",
                      100 * count_correct / i, 100 * count_graph / i, i, reps))
        }
      }
      
      pct_dist <- 100 * count_correct / reps
      pct_graph <- 100 * count_graph / reps
      pct_oracle_graph <- 100 * count_oracle_graph / reps
      if (i %% 10 == 0){cat(sprintf("✅ Finished: n = %d | family = %d | dist = %.1f%% | graph = %.1f%%\n\n",
                                    n, fam, pct_dist, pct_graph))}
      
      list(pct_dist = pct_dist, pct_graph = pct_graph, pct_oracle_graph=pct_oracle_graph)
    })
  ) %>%
  unnest_wider(correct_pct)

# Clean for plotting
results_long <- results %>%
  mutate(family_label  = factor(fam, levels = c(1, 2),
                                labels = c("Family 1", "Family 2"))) %>%
  pivot_longer(cols = c(pct_dist, pct_graph, pct_oracle_graph),
               names_to = "metric", values_to = "accuracy") %>%
  mutate(metric = recode(metric,
                         pct_dist = "Correct distribution",
                         pct_graph = "Sequential",
                         pct_oracle_graph = "Oracle"))

# Plot 1: Correct distribution (no legend)
p1 <- results_long %>%
  filter(metric == "Correct distribution") %>%
  ggplot(aes(x = n, y = accuracy, color = family_label)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Sequential approach:\nCorrect family choice",
    x = "Sample size (n)",
    y = "% correct"
  ) +
  ylim(0, 100) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

# Plot 2: Graph = '1 --> 2' (legend updated)
p2 <- results_long %>%
  filter(metric %in% c("Sequential", "Oracle")) %>%
  ggplot(aes(x = n, y = accuracy, color = family_label, linetype = metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Accuracy:\nSequential vs Oracle",
    x = "Sample size (n)",
    y = "% correct",
    linetype = "Method",
    color = "Generating \nfamily"
  ) +
  ylim(50, 100) +
  theme_minimal(base_size = 14)


p2 + p1  # patchwork combines them horizontally

# Save the combined plot
ggsave("Sequ.appr.pdf", width = 8.3, height = 3.3, dpi = 500)
