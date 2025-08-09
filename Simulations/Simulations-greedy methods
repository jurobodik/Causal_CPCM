##################### Compare Greedy Methods for DAG Learning #####################
library(bnlearn)
library(SID)
library(pcalg)
library(igraph)
library(graph)      # for graphNEL
library(RBGL)       # needed by igraph.from.graphNEL()
library(ggplot2)
library(gridExtra)
library(patchwork)


######################### Utils for Graph Operations #########################
compare_dags <- function(learned_dag, true_dag) {
  
  old_names <- nodes(learned_dag)
  new_names <- paste0("X", seq_along(old_names))
  name_map <- setNames(new_names, old_names)
  learned_dag <- rename.nodes(learned_dag, name_map)
  
  
  old_names <- nodes(true_dag)
  new_names <- paste0("X", seq_along(old_names))
  name_map <- setNames(new_names, true_dag)
  true_dag <- rename.nodes(true_dag, name_map)
  
  
  learned_mat <- amat(learned_dag)
  true_mat <- amat(true_dag)
  
  # Union of all nodes
  all_nodes <- union(rownames(learned_mat), rownames(true_mat))
  
  # Ensure both matrices are in the same node order
  learned_mat <- learned_mat[all_nodes, all_nodes, drop = FALSE]
  true_mat <- true_mat[all_nodes, all_nodes, drop = FALSE]
  
  # Extract directed edges
  learned_edges <- which(learned_mat == 1, arr.ind = TRUE)
  true_edges <- which(true_mat == 1, arr.ind = TRUE)
  
  # Convert to edge lists
  edge_list <- function(edges, nodes) {
    apply(edges, 1, function(x) paste0(nodes[x[1]], "->", nodes[x[2]]))
  }
  learned_list <- edge_list(learned_edges, all_nodes)
  true_list <- edge_list(true_edges, all_nodes)
  
  # Compare
  correct_edges <- intersect(learned_list, true_list)
  missing_edges <- setdiff(true_list, learned_list)
  extra_edges <- setdiff(learned_list, true_list)
  
  # Output
  cat("=== DAG Comparison ===\n")
  cat("Correct edges: ", length(correct_edges), "\n")
  cat("Missing edges: ", length(missing_edges), "\n")
  if (length(missing_edges)) cat("  ", paste(missing_edges, collapse = ", "), "\n")
  cat("Extra edges: ", length(extra_edges), "\n")
  if (length(extra_edges)) cat("  ", paste(extra_edges, collapse = ", "), "\n")
  
  # Return lists for further use
  invisible(list(correct = correct_edges, missing = missing_edges, extra = extra_edges))
}

generate_random_dag <- function(d, p) {
  # Ensure p is within valid range
  max_edges <- d * (d - 1) / 2  # Maximum edges for a DAG
  if (p > max_edges) {
    stop("Too many edges for a DAG with this number of nodes.")
  }
  
  # Create node names as character strings
  nodes <- as.character(seq_len(d))
  
  # Generate a random topological ordering of nodes
  node_order <- sample(nodes)  # Random order ensures no cycles
  
  # Generate all possible edges following this order
  possible_edges <- expand.grid(from = node_order, to = node_order, stringsAsFactors = FALSE)
  possible_edges <- possible_edges[possible_edges$from != possible_edges$to, ]
  
  # Keep only edges that respect the topological order
  possible_edges <- possible_edges[match(possible_edges$from, node_order) < 
                                     match(possible_edges$to, node_order), ]
  
  # Randomly select p edges
  selected_edges <- possible_edges[sample(nrow(possible_edges), p), ]
  
  # Create an empty DAG
  dag <- empty.graph(nodes)
  
  # Add edges to the DAG
  for (i in 1:nrow(selected_edges)) {
    dag <- set.arc(dag, from = as.character(selected_edges$from[i]), 
                   to = as.character(selected_edges$to[i]))
  }
  
  old_names <- nodes(dag)
  new_names <- paste0("X", seq_along(old_names))
  name_map <- setNames(new_names, old_names)
  dag <- rename.nodes(dag, name_map)
  
  return(dag)
}

generate_random_dag_matrix_with_equal_prob <- function(d, edge_prob = 0.5) {
  mat <- matrix(0, nrow = d, ncol = d)
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      mat[i, j] <- rbinom(1, 1, edge_prob)
    }
  }
  
  # Shuffle the vertices (rows and columns)
  perm <- sample(d)
  mat <- mat[perm, perm]
  
  # Restore variable names to X1, ..., Xd (in order, not permuted)
  var_names <- paste0("X", 1:d)
  rownames(mat) <- colnames(mat) <- var_names
  
  return(mat)
}

generate_random_scm <- function(n = 500, dag,
                                scenario = c('Additive_Gaussian_sin',
                                             'CPCM_exponential_linear',
                                             'LINGAM_linear',
                                             'CPCM_exp_gauss')) {
  if (!inherits(dag, "bn")) stop("dag must be a bnlearn DAG object.")
  
  scenario <- match.arg(scenario)
  
  # Define basic building blocks
  scenario_map <- list(
    Additive_Gaussian_sin = list(
      f = function(x) rowSums(sin(x) + 0.5 * x),
      gen = function(fx) fx + rnorm(n, sd = 0.5)
    ),
    CPCM_exponential_linear = list(
      f = function(x) 1/rowSums(abs(x)),
      gen = function(fx) rexp(n, rate = pmax(0.1, fx))
    ),
    LINGAM_linear = list(
      f = function(x) rowSums(x),
      gen = function(fx) fx + rexp(n, rate = 1)
    )
  )
  
  # Prepare DAG structure
  nodes <- nodes(dag)
  parents_list <- lapply(nodes, function(v) bnlearn::parents(dag, v))
  names(parents_list) <- nodes
  ordering <- bnlearn::node.ordering(dag)
  data <- matrix(NA, nrow = n, ncol = length(nodes))
  colnames(data) <- nodes
  
  for (node in ordering) {
    pa <- parents_list[[node]]
    if (length(pa) == 0) {
      data[, node] <- rnorm(n)
    } else {
      pa_data <- data[, pa, drop = FALSE]
      
      if (scenario != "CPCM_exp_gauss") {
        fx <- scenario_map[[scenario]]$f(pa_data)
        data[, node] <- scenario_map[[scenario]]$gen(fx)
        
      } else {
        # Hybrid case: choose randomly between exponential and Gaussian
        use_exp <- rbinom(1, 1, 0.5) == 1
        if (use_exp) {
          fx <- scenario_map$CPCM_exponential_linear$f(pa_data)
          data[, node] <- scenario_map$CPCM_exponential_linear$gen(fx)
        } else {
          fx <- scenario_map$Additive_Gaussian_sin$f(pa_data)
          data[, node] <- scenario_map$Additive_Gaussian_sin$gen(fx)
        }
      }
    }
  }
  
  colnames(data) <- paste0("X", seq_len(ncol(data)))
  as.data.frame(data)
}

sid_distance <- function(true, estimate){
  if(inherits(true, "bn")){true = amat(true)}
  if(inherits(estimate, "bn")){estimate = amat(estimate)}
  
  m=abs(true-estimate)
  # Count where both m[i,j] and m[j,i] == 1
  bidirected_count_estimate <- sum(estimate[upper.tri(estimate)] & t(estimate)[upper.tri(estimate)])
  bidirected_count_true <- sum(true[upper.tri(true)] & t(true)[upper.tri(true)])
  
  return(sum(m) - bidirected_count_estimate -  bidirected_count_true)
}


######################### ######################### ######################### ######################### 
######################### ######################### ######################### ######################### 
######################### ######################### ######################### ######################### 
######################### ######################### ######################### ######################### 
# Main script to compare greedy methods for DAG learning
n <- 1000
ds <- c(3, 4, 6, 8, 10)
n_reps <- 50
methods <- c("RESIT", "edge_greedy", "RESIT_greedy", "exact")
scenario = "CPCM_exponential_linear"
results <- list()

for (d in ds) {
  for (method in methods) {
    
    # Skip 'exact' method for d > 4
    if (method == "exact" && !(d %in% c(3, 4))) {
      next
    }
    
    cat("d =", d, " | method =", method, "\n")
    
    sid_vals <- c()
    times <- c()
    
    for (rep in 1:n_reps) {
      cat("Rep", rep, "\n")
      p = sample(1:(d*(d-1)/2), 1, prob=exp(-0.1*(1:(d*(d-1)/2))))
      dag <- generate_random_dag(d = d, p = p )
      X <- generate_random_scm(n, dag, scenario = scenario)
      
      result <- try({
        start_time <- Sys.time()
        est <- CPCM_graph_estimate(X, method = method)
        end_time <- Sys.time()
        
        sid <- sid_distance(est, dag)
        time_taken <- as.numeric(difftime(end_time, start_time, units = "secs"))
        
        sid_vals <- c(sid_vals, sid)
        times <- c(times, time_taken)
      }, silent = TRUE)
      
      if (inherits(result, "try-error")) {
        break
      }
    }
    
    results[[paste0("d=", d, "_", method)]] <- list(
      sid = mean(sid_vals),
      time = mean(times)
    )
  }
}

# Convert to data frame for plotting
df <- data.frame(
  d = numeric(),
  sid = numeric(),
  time = numeric(),
  method = character()
)

for (name in names(results)) {
  parts <- strsplit(name, "_")[[1]]
  d_val <- as.numeric(gsub("d=", "", parts[1]))
  method <- parts[2]
  df <- rbind(df, data.frame(
    d = d_val,
    sid = results[[name]]$sid,
    time = results[[name]]$time,
    method = method
  ))
}

df[19,] = df[1,]
df[19,]$method = "Trivial"
df[19,]$sid = 1
df[19,]$time = 0
df[19,]$d = 3
df[20,] = df[19,]
df[20,]$d = 10


p1 <- ggplot(df, aes(x = d, y = sid, color = method)) +
  geom_line() +
  geom_point() +
  labs(title = "Normalized SID vs. number of nodes", y = "SID/d", x = "Number of nodes") +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )


p2 <- ggplot(df, aes(x = d, y = time, color = method)) +
  geom_line() +
  geom_point() +
  scale_y_log10() +
  coord_cartesian(ylim = c(NA, 2000)) +
  labs(title = "Runtime vs. number of nodes", y = "Time (s, log scale)", x = "Number of nodes") +
  theme(
    plot.title = element_text(size = 15),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 13),
    legend.position = "none"
  )

p1 + p2

#save 
ggsave("greedy_comp.pdf", width = 10, height = 4, dpi = 500)



