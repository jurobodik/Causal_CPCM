# just some functions for graph manipulation, 
# and PC, GES, LINGAM... definition wrappers
library(bnlearn)
library(SID)
library(pcalg)
library(igraph)
library(graph)      # for graphNEL
library(RBGL)       # needed by igraph.from.graphNEL()
library(ggplot2)
library(gridExtra)

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

generate_random_dag <- function(d, number_of_edges = NULL, prob = NULL) # Use either `number_of_edges` (fixed number of edges) or `prob` (independent inclusion probability) to generate the DAG.
{
  # --- checks ---
  if (!is.null(prob) && (!is.numeric(prob) || length(prob) != 1 || prob < 0 || prob > 1)) {
    stop("`prob` must be a single number in [0, 1].")
  }
  if (!is.null(number_of_edges) && (!is.numeric(number_of_edges) || number_of_edges < 0 || floor(number_of_edges) != number_of_edges)) {
    stop("`number_of_edges` must be a non-negative integer.")
  }
  if (is.null(number_of_edges) && is.null(prob)) {
    stop("Provide either `number_of_edges` or `prob`.")
  }
  if (!is.null(number_of_edges) && !is.null(prob)) {
    warning("Both `number_of_edges` and `prob` supplied; using `prob` and ignoring `number_of_edges`.")
  }
  
  # node names
  nodes_chr <- as.character(seq_len(d))
  
  # random topological order (uniform over permutations)
  order <- sample(nodes_chr)
  
  # all edges consistent with this order
  possible_edges <- do.call(rbind, lapply(seq_len(d - 1), function(i) {
    from <- order[i]
    to   <- order[(i + 1):d]
    data.frame(from = from, to = to, stringsAsFactors = FALSE)
  }))
  max_edges <- nrow(possible_edges) # equals d*(d-1)/2
  
  # pick edges
  if (!is.null(prob)) {
    keep <- runif(max_edges) < prob
    selected_edges <- possible_edges[keep, , drop = FALSE]
  } else {
    if (number_of_edges > max_edges) stop("Too many edges for a DAG with this number of nodes.")
    pick <- sample.int(max_edges, number_of_edges)
    selected_edges <- possible_edges[pick, , drop = FALSE]
  }
  
  # build DAG (bnlearn)
  dag <- empty.graph(nodes_chr)
  if (nrow(selected_edges) > 0) {
    for (i in seq_len(nrow(selected_edges))) {
      dag <- set.arc(dag, from = selected_edges$from[i], to = selected_edges$to[i])
    }
  }
  
  # rename nodes to X1,...,Xd
  old_names <- nodes(dag)
  new_names <- paste0("X", seq_along(old_names))
  dag <- rename.nodes(dag, setNames(new_names, old_names))
  
  return(dag)
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

PC_wrapper = function(X){
  pc_fit <- pc(list(C = cor(X), n = nrow(X)), indepTest = gaussCItest,
               alpha = 0.05, labels = colnames(X))
  pc_bn <- as.bn(as(pc_fit@graph, "graphNEL"))
  pc_mat   <- as(as.graphNEL(pc_bn), "matrix")
  return(pc_mat)
}

GES_wrapper=function(X){
  ges_fit <- ges(new("GaussL0penObsScore", data = X))
  ges_bn <- as.bn(as(ges_fit$essgraph, "graphNEL"))
  ges_mat  <- as(as.graphNEL(ges_bn), "matrix")
  return(ges_mat)
}

LINGAM_wrapper = function(X){
  lingam_mat <- lingam(X)$B != 0  # adjacency matrix
  return(t(lingam_mat))
}

CPCM_wrapper = function(X, family_of_distributions = 1){
  cpcm_bn = CPCM_graph_estimate(X, family_of_distributions = family_of_distributions)
  return(cpcm = as(as.graphNEL(cpcm_bn), "matrix"))
}

ANM_RESIT_wrapper = function(X){
  anm = CPCM_graph_estimate(X, family_of_distributions = 'Gaussian with fixed sigma', greedy_method = 'RESIT')
  return(cpcm = as(as.graphNEL(anm), "matrix")) 
  
}
