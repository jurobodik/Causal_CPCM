library(ggplot2)
library(pcalg)
library(igraph)
library(graph)      # for graphNEL
library(RBGL)       # needed by igraph.from.graphNEL()
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)


source('CPCM_function.R')
source('utils_for_graphs.R')
#############################  Example #################################
# ----------- DAG + SCM generation -----------
n <- 1000
d= 4; prob = 2/(d-1) # Probability of edge inclusion in the random DAG
true_dag <- generate_random_dag(d=d, prob = prob)
X <- generate_random_scm(n = n, dag = true_dag, scenario = 'LINGAM_linear') #scenarios= c('Additive_Gaussian_sin', 'CPCM_exponential_linear', 'LINGAM_linear', 'CPCM_exp_gauss')
# ----------- graph estimation -----------
cpcm_estimate = CPCM_wrapper(X, 's')
pc_estimate=PC_wrapper(X, 'hsic') #change to 'gauss' for a different independence test
ges_estimate=GES_wrapper(X)
lingam_estimate = LINGAM_wrapper(X)
anm_estimate = ANM_RESIT_wrapper(X)
random_estimate <- generate_random_dag(d = d, prob = 0.5)

# -----------  SID  --------------------
true_dag_adj =  amat(true_dag)
sid_cpcm <- sid_distance(true_dag_adj, cpcm_estimate)
sid_pc  <- sid_distance(true_dag_adj, pc_estimate)
sid_ges <- sid_distance(true_dag_adj, ges_estimate)
sid_lingam <- sid_distance(true_dag_adj, lingam_estimate)
sid_anm <- sid_distance(true_dag_adj, anm_estimate)
sid_random <- sid_distance(true_dag_adj, random_estimate)

# ----------- Report results -----------
cat(
  "\nSID Results:\n",
  "CPCM     - SID:", sid_cpcm, "\n",
  "PC       - SID:", sid_pc, "\n",
  "GES      - SID:", sid_ges, "\n",
  "LINGAM   - SID:", sid_lingam, "\n",
  "ANM      - SID:", sid_anm, "\n",
  "Random   - SID:", sid_random, "\n"
)


############################## Repeating the above for multiple scenarios and methods ##############################
scenarios <- c('Additive_Gaussian_sin', 'CPCM_exponential_linear', 
               'LINGAM_linear', 'CPCM_exp_gauss')
methods <- c('PC', 'GES', 'LINGAM', 'ANM', 'ANM_greedy', 'CPCM', 'Random')
n_iter <- 50
n <- 1000
d <- 5
prob <- 2/(d-1)

# Store results
all_results <- data.frame()

total_tasks <- length(scenarios) * n_iter
task_counter <- 0
start_time <- Sys.time()

safe_run <- function(expr) {
  result <- tryCatch(expr, error = function(e) NULL)
  return(result)
}

for (sc in scenarios) {
  cat("Starting scenario:", sc, "\n")
  
  for (i in 1:n_iter) {
    task_counter <- task_counter + 1
    
    repeat {  # retry loop if any method fails
      success <- TRUE
      true_dag <- generate_random_dag(d = d, prob = prob)
      X <- generate_random_scm(n = n, dag = true_dag, scenario = sc)
      true_dag_adj <- amat(true_dag)
      
      # Run all methods safely
      estimates <- list(
        PC        = safe_run(PC_wrapper(X, 'gauss')), #Change to 'hsic' for a different independence test
        GES       = safe_run(GES_wrapper(X)),
        LINGAM    = safe_run(LINGAM_wrapper(X)),
        ANM       = safe_run(ANM_RESIT_wrapper(X)),
        ANM_greedy= safe_run(ANM_RESIT_greedy_wrapper(X)),
        CPCM      = safe_run(CPCM_wrapper(X, 's')),
        Random    = safe_run(generate_random_dag(d = d, prob = 0.5))
      )
      
      # If any method failed, redo iteration
      if (any(sapply(estimates, is.null))) {
        success <- FALSE
        next
      }
      
      # Otherwise compute SIDs
      for (m in names(estimates)) {
        sid <- sid_distance(true_dag_adj, estimates[[m]])
        all_results <- rbind(all_results, data.frame(
          Scenario = sc,
          Iteration = i,
          Method = m,
          SID = sid
        ))
      }
      
      break # leave repeat loop if success
    }
    
    # Progress report
    elapsed <- Sys.time() - start_time
    avg_time <- as.numeric(elapsed, units = "secs") / task_counter
    remaining_time <- (total_tasks - task_counter) * avg_time
    cat(sprintf("Completed %d of %d (%.1f%%) — Elapsed: %.1f sec — ETA: %.1f sec\n",
                task_counter, total_tasks, 100 * task_counter / total_tasks,
                as.numeric(elapsed, units = "secs"), remaining_time))
  }
}


avg_results <- all_results %>%
  group_by(Scenario, Method) %>%
  summarise(mean_SID = mean(SID), .groups = "drop") %>%
  mutate(Method = factor(Method, levels = c("CPCM", "LINGAM", "ANM", "ANM_greedy", "PC", "GES", "Random"))) %>%
  arrange(Method) %>%
  pivot_wider(names_from = Scenario, values_from = mean_SID)

print(avg_results)

