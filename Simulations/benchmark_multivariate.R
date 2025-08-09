library(ggplot2)
library(pcalg)
library(igraph)
library(graph)      # for graphNEL
library(RBGL)       # needed by igraph.from.graphNEL()
library(ggplot2)
library(gridExtra)

source('utils_for_graphs.R')
#############################  Example #################################
# ----------- DAG + SCM generation -----------
n <- 1000
d=5;p=3
true_dag <- generate_random_dag(d=d, p=p)
X <- generate_random_scm(n = n, dag = true_dag, scenario = 'LINGAM_linear') #scenarios= c('Additive_Gaussian_sin', 'CPCM_exponential_linear', 'LINGAM_linear', 'CPCM_exp_gauss')
# ----------- graph estimation -----------
pc_estimate=PC_wrapper(X)
ges_estimate=GES_wrapper(X)
lingam_estimate = LINGAM_wrapper(X)
anm_estimate = ANM_RESIT_wrapper(X)
cpcm_estimate = CPCM_wrapper(X, 1)
random_estimate <- generate_random_dag_matrix_with_equal_prob(d = d, edge_prob = 0.5)

# -----------  SID  --------------------
true_dag_adj =  amat(true_dag)
sid_pc  <- sid_distance(true_dag_adj, pc_estimate)
sid_ges <- sid_distance(true_dag_adj, ges_estimate)
sid_lingam <- sid_distance(true_dag_adj, lingam_estimate)
sid_anm <- sid_distance(true_dag_adj, anm_estimate)
sid_cpcm <- sid_distance(true_dag_adj, cpcm_estimate)
sid_random <- sid_distance(true_dag_adj, random_estimate)

# ----------- Report results -----------
cat(
  "\nSID Results:\n",
  "PC       - SID:", sid_pc, "\n",
  "GES      - SID:", sid_ges, "\n",
  "LINGAM   - SID:", sid_lingam, "\n",
  "ANM      - SID:", sid_anm, "\n",
  "CPCM     - SID:", sid_cpcm, "\n",
  "Random   - SID:", sid_random, "\n"
)




############################## Repeating the above for multiple scenarios and methods ##############################



# Define scenarios and methods
scenarios <- c('Additive_Gaussian_sin', 'CPCM_exponential_linear', 'LINGAM_linear', 'CPCM_exp_gauss')
methods <- c('PC', 'GES', 'LINGAM', 'ANM', 'CPCM', 'Random')
n_iter <- 50
n <- 1000
d <- 8

# Store all results
all_results <- data.frame()

total_tasks <- length(scenarios) * n_iter
task_counter <- 0
start_time <- Sys.time()

for (sc in scenarios) {
  cat("Starting scenario:", sc, "\n")
  
  for (i in 1:n_iter) {
    task_counter <- task_counter + 1
    p <- sample(1:5, 1)
    true_dag <- generate_random_dag(d = d, p = p)
    X <- generate_random_scm(n = n, dag = true_dag, scenario = sc)
    true_dag_adj <- amat(true_dag)
    
    # Run all methods
    estimates <- list(
      PC = PC_wrapper(X),
      GES = GES_wrapper(X),
      LINGAM = LINGAM_wrapper(X),
    #  ANM = ANM_RESIT_wrapper(X),
      CPCM = CPCM_wrapper(X, 1),
      Random = generate_random_dag_matrix_with_equal_prob(d = d, edge_prob = 0.5)
    )
    
    # Compute SIDs and store
    for (m in names(estimates)) {
      sid <- sid_distance(true_dag_adj, estimates[[m]])
      all_results <- rbind(all_results, data.frame(
        Scenario = sc,
        Iteration = i,
        Method = m,
        SID = sid
      ))
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

# Summary table
sid_summary <- aggregate(SID ~ Scenario + Method, data = all_results, FUN = mean)
print(sid_summary)

# Plot boxplots (one per scenario)
for (sc in scenarios) {
  p <- ggplot(subset(all_results, Scenario == sc), aes(x = Method, y = SID)) +
    geom_boxplot() +
    ggtitle(paste("SID Boxplot - Scenario:", sc)) +
    theme_minimal()
  print(p)
}
