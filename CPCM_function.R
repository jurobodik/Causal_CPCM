#This is the code for the main function that estimates the causal graph via CPCM, called CPCM_graph_estimate(X, family_of_distributions). 

#family_of_distributions correspond to the models we use. If you want to use CPCM(F) model, the choices for 'family_of_distributions' are
#are the following: "Gaussian", "Gaussian with fixed sigma", "Pareto", "Pareto2", "Exponential", "Gamma", "Gamma with fixed scale", "Gumbel", "Gumbel with fixed scale", "Poisson", "Negative_binomial"


#We implemented a joint CPCM(F1...Fk) model with the following:
#family_of_distributions = 1 if we use family1 = c('Gaussian with fixed sigma', 'Poisson', 'Exponential', 'Pareto')
#family_of_distributions = 2 if we use family2 = c('Gaussian', 'Negative_binomial','Gamma', 'Pareto2')
#family_of_distributions = 'Sequential choice' if we first try family_of_distributions = 1, test plausibility and if unplausible then family_of_distributions = 2
#however, it is more recommended to use family_of_distributions = 1 if n<=1000 as a rule of thumb


#If you want different family than already coded or different method than implemented, it is easy! For that,
#simply add into function "estimate_epsilon(Y, X, family="YOUR DISTRIBUTION FUNCTION")" a code how you want to estimate your epsilons
#We always used "gam" function from "MGCV" package for the (smooth) parameter estimation. 

library("mgcv") 
library(dHSIC)
library(bnlearn)
library(MASS) #Used indirectly via nb() in gam() for Negative Binomial family
library("gamlss") #	Provides additional distribution families (gaulss, gammals, gumbls) used in gam()
library(stringr) #Used for sub() string operations in the "Sequential choice" logic
library(dplyr) #I dont think I used it, but its always useful for data manipulation

##################################   Example   ###########################################
#n=1000
#X1 = rnorm(n)
#X2 = rnorm(n,X1^2,1)
#X  = data.frame(X1, X2)
#CPCM_graph_estimate(X, family_of_distributions = 'Sequential choice') 

#X3 = rexp(n, rate = 1/X2^2)
#X  = data.frame(X1, X2, X3)
#CPCM_graph_estimate(X, family_of_distributions = 1) #Should return DAG (bnlearn type) with edges [X1][X2|X1][X3|X2]
############################ CPCM estimation of the causal graph ##########################


CPCM_graph_estimate <- function(X, family_of_distributions = 1, greedy_method = 'RESIT_greedy', lambda = 1, 
                                quiet = TRUE){  #exact, RESIT_greedy, edge_greedy, RESIT
  
  n = length(X[,1])
  d = ncol(X)
  
  if(family_of_distributions=='S' || family_of_distributions=='s')family_of_distributions = 'Sequential choice'
  
  #Estimation \hat{S}
  estimate_support_of_Y_and_pair_with_family<-function(Y, family_of_distributions = 1){
    
    determine_support <- function(Y) {
      # Check if the data is discrete (integer values only)
      if (all(Y == floor(Y)) & length(unique(Y)) < length(Y) / 10) {
        return('discrete')  # Poisson distribution
      }
      
      # Check if the data is within the interval [0, 1]
      if (all(Y >= 0 & Y <= 1)) {
        return('interval')  # Beta distribution
      }
      
      # Check for support on [1, ∞) and estimate tail index
      if (all(Y >= 1)) {
        skewness= (quantile(Y, 0.90) + quantile(Y, 0.10) - 2 * quantile(Y, 0.50)) / (quantile(Y, 0.90) - quantile(Y, 0.10))
        if(skewness > 0.2){  
          # Is power-decay is better fit than exponential decay?
          tail_frac = 0.1
          threshold <- round((1 - tail_frac) * length(Y)):length(Y)
          Y_tail <- sort(Y)[threshold]
          S_tail <- 1 - (threshold / length(Y))
          # Avoid log(0) issues
          S_tail <- ifelse(S_tail <= 0, min(S_tail[S_tail > 0]), S_tail)
          # Fit log-log (power law)
          fit_pl <- lm(log(S_tail) ~ log(Y_tail))
          R2_pl <- summary(fit_pl)$r.squared
          # Fit log-linear (exponential)
          fit_exp <- lm(log(S_tail) ~ Y_tail)
          R2_exp <- summary(fit_exp)$r.squared
          # Decision rule
          if (R2_pl > R2_exp) {return("power tail")}
        }
      }
      
      # Check for Gamma-type characteristics (positive skewed)
      if (all(Y > 0)) {
        skewness= (quantile(Y, 0.90) + quantile(Y, 0.10) - 2 * quantile(Y, 0.50)) / (quantile(Y, 0.90) - quantile(Y, 0.10))
        if (skewness > 0.2) {
          return('half-line')  # Gamma or Exponential
        }
      }
      
      # Default: assume Gaussian
      return('full support')
    }
    
    
    if(family_of_distributions == 1){
      if( determine_support(Y)=='full support') return('Gaussian with fixed sigma')
      if( determine_support(Y)=='power tail') return('Pareto')
      if( determine_support(Y)=='discrete') return('Poisson')
      if( determine_support(Y)=='interval') return('Gaussian with fixed sigma') #Beta distribution is no longer supported :(
      if( determine_support(Y)=='half-line') return('Exponential')
    }
    
    if(family_of_distributions == 2){
      if( determine_support(Y)=='full support') return('Gaussian')
      if( determine_support(Y)=='power tail') return('Pareto2')
      if( determine_support(Y)=='discrete') return('Negative_binomial')
      if( determine_support(Y)=='interval') return('Gaussian') #Beta distribution is no longer supported  :(
      if( determine_support(Y)=='half-line') return('Gamma')
    }
    
    
  }
  
  estimate_epsilon <- function(Y, X, family = "Gaussian") {
    
    method = "smooth" # Default method for GAM formula construction; change if you want linear model
    # --- Helper: Rename X columns to X1, X2, ...
    rename_columns <- function(X) {
      if (is.null(nrow(X))) {
        X <- as.data.frame(X)
        names(X)[1] <- "X1"
      } else {
        names(X) <- paste0("X", seq_along(X))
      }
      return(X)
    }
    
    # --- Helper: Construct GAM formula based on method and variable type
    construct_formula <- function(X, include_Y = TRUE, method = "smooth", pb = FALSE) {
      terms <- sapply(seq_along(X), function(i) {
        xi <- X[[i]]
        if (length(unique(xi)) <= 9) {
          paste0("as.factor(X", i, ")")
        } else if (method == "smooth") {
          if(pb == FALSE){paste0("s(X", i, ")")}else{paste0("pb(X", i, ")")}
        } else {
          paste0("X", i)
        }
      })
      rhs <- paste(terms, collapse = " + ")
      as.formula(if (include_Y) paste("Y ~", rhs) else paste("~", rhs))
    }
    
    # --- Data preparation
    X <- rename_columns(X)
    Y <- setNames(as.data.frame(Y), "Y")
    data <- cbind(Y, X)
    
    # --- Reuse formulas
    formula_mu <- construct_formula(X, include_Y = TRUE, method = method)
    formula_sigma <- construct_formula(X, include_Y = FALSE, method = method)
    
    # --- Dispatcher via switch
    residuals <- switch(family,
                        
                        "Gaussian" = {
                          fit <- gam(list(formula_mu, formula_sigma), data = data, family = gaulss(), link = list("identity", "logb"))
                          mu <- fit$fitted.values[, 1]
                          inv_sigma <- fit$fitted.values[, 2]
                          pnorm((Y$Y - mu) * inv_sigma)
                        },
                        
                        "Gaussian with fixed sigma" = {
                          fit <- gam(list(formula_mu, ~1), data = data, family = gaulss(), link = list("identity", "logb"))
                          mu <- fit$fitted.values[, 1]
                          inv_sigma <- fit$fitted.values[, 2]
                          pnorm((Y$Y - mu) * inv_sigma)
                        },
                        
                        "Gamma" = {
                          fit <- gam(list(formula_mu, formula_sigma), data = data, family = gammals())
                          shape <- 1 / exp(fit$fitted.values[, 2])
                          scale <- fit$fitted.values[, 1] * exp(fit$fitted.values[, 2])
                          pgamma(Y$Y, shape = shape, scale = scale)
                        },
                        
                        "Gamma with fixed scale" = {
                          fit <- gam(list(formula_mu, ~1), data = data, family = gammals())
                          shape <- 1 / exp(fit$fitted.values[, 2])
                          scale <- fit$fitted.values[, 1] * exp(fit$fitted.values[, 2])
                          pgamma(Y$Y, shape = shape, scale = scale)
                        },
                        
                        "Exponential" = {
                          fit <- gam(formula_mu, data = data, family = Gamma(link = "log"))
                          alpha_hat = 1 / mgcv::predict.gam(fit, type = 'response')  # This gives mean: mu(x) = 1 / alpha(x)
                          pexp(data$Y, rate = alpha_hat)
                        },
                        
                        "Pareto" = {
                          data$Y <- log(Y$Y) #log(Pareto) = Exponential distribution
                          fit <- gam(formula_mu, data = data, family = Gamma(link = "log"))
                          alpha_hat = 1 / mgcv::predict.gam(fit, type = 'response')  # This gives mean: mu(x) = 1 / alpha(x)
                          pexp(data$Y, rate = alpha_hat)
                        },
                        
                        "Pareto2" = {
                          data$Y <- log(Y$Y) #log(Pareto2) = Gamma distribution
                          fit <- gam(list(formula_mu, formula_sigma), data = data, family = gammals())
                          shape <- 1 / exp(fit$fitted.values[, 2])
                          scale <- fit$fitted.values[, 1] * exp(fit$fitted.values[, 2])
                          pgamma(Y$Y, shape = shape, scale = scale)
                        },
                        
                        "Poisson" = {
                          my_ppois <- function(Y, lambda) {
                            ppois(Y - 1, lambda = lambda) + runif(length(Y)) * dpois(Y, lambda = lambda)
                          }
                          fit <- gam(formula_mu, data = data, family = poisson())
                          my_ppois(Y$Y, lambda = fitted(fit))
                        },
                        
                        "Negative_binomial" = { #A bit tricky, sometimes gives error and doesnt work very well :(
                          to_mgcv <- function(f) as.formula(gsub("pb\\(([^)]+)\\)", "s(\\1, bs='ps')", deparse(f)), env = environment(f))
                          
                          fit <- gam(to_mgcv(construct_formula(X, TRUE,  method, pb = TRUE)), family = nb(), data = data, method = "REML")
                          theta <- tryCatch({
                            th <- fit$family$getTheta(TRUE); 
                            if (length(th)!=1 || !is.finite(th)) stop("bad"); th
                          }, error = function(e) {
                            th <- try(fit$family$getTheta(TRUE), silent=TRUE)
                            if (inherits(th,"try-error")) exp(get(".Theta", envir = environment(fit$family$variance))) else th
                          })
                          
                          param <- list(mu = as.numeric(fitted(fit)),
                                        size = rep(theta, NROW(data)))
                          
                          my_pnbinom <- function(Y, mu, size) {
                            Y <- as.numeric(Y)
                            pnbinom(Y - 1, mu = mu, size = size) +
                              runif(length(Y)) * dnbinom(Y, mu = mu, size = size)
                          }
                          YY=Y$Y
                          mu <- param$mu
                          size <- param$size
                          epsilon <- my_pnbinom(YY, mu = mu, size = size)
                          return(epsilon)
                        },
                        
                        "Gumbel" = {
                          fit <- gam(list(formula_mu, formula_sigma), data = data, family = gumbls())
                          mu <- fit$fitted.values[, 1]
                          scale <- exp(fit$fitted.values[, 2])
                          pgumbel <- function(x, location = 0, scale = 1) {exp(-exp(-(x - location) / scale))}
                          pgumbel(Y$Y, mu, scale)
                        },
                        
                        "Gumbel with fixed scale" = {
                          fit <- gam(list(formula_mu, ~1), data = data, family = gumbls())
                          mu <- fit$fitted.values[, 1]
                          scale <- exp(fit$fitted.values[, 2])
                          pgumbel <- function(x, location = 0, scale = 1) {exp(-exp(-(x - location) / scale))}
                          pgumbel(Y$Y, mu, scale)
                        },
                        
                        stop("Family not implemented.")
    )
    
    return(residuals)
  }
  
  
  create_bivariate_output_from_p_values <- function(indep, z1, z2, family1, family2) {
    #Score-based graph estimate
    empty = -log(indep) 
    one = -log(z1) + lambda
    two = -log(z2) + lambda
    if (empty < one & empty < two) {res1 = "Empty graph"} else if (one < two) {res1 = "1 --> 2"} else {res1 = "2 --> 1"}
    
    #Testing estimate
    res2=res1
    if(indep>0.05) {res2 = "Empty graph"}else{
      if (z1>0.05 & z2>0.05) res2 = "Unidentifiable  (both directions are plausible)" 
      if (z1<0.05 & z2<0.05) res2 = "Assumptions not fulfilled (both directions are not plausible)"
      if (z1<0.05 & z2>0.05) res2 = "2 --> 1"
      if (z1>0.05 & z2<0.05) res2 = "1 --> 2"
    }
    
    #Forced estimate
    if(z1>=z2) {res3 = "1 --> 2"} else {res3 = "2 --> 1"}
    
    res=as.data.frame(c(round(indep, digits=6), round(z1, digits=6), round(z2, digits = 6), res1, res2, res3, paste0(family1, ";", family2)))
    colnames(res)<-c("Results")
    rownames(res)<-c("p-value Empty graph", "p-value 1-->2", "p-value 2-->1", "Score-based graph estimate", 'Testing estimate', 'Forced estimate', 'Families used' )
    return(res)}
  
  
  bivariate_CPCM_graph_estimate <- function(X, family_of_distributions = 1){
    
    #Are X1 and X2 idependent? If yes, we return an empty graph as an estimate
    X1 = X[,1]; X2 = X[,2]
    indep=dhsic.test(data.frame(X1, X2), method = 'gamma')$p.value
    
    
    #X1-->X2 direction
    Y = X2; X_=data.frame(X1 = X1)
    if(family_of_distributions == 1 || family_of_distributions==2)
    {family2 = estimate_support_of_Y_and_pair_with_family(Y, family_of_distributions = family_of_distributions)}else{family2 = family_of_distributions}
    r2=as.numeric(estimate_epsilon(Y, X_, family = family2))
    
    #X2-->X1 direction
    Y= X1; X_=data.frame(X1 = X2)
    if(family_of_distributions == 1 || family_of_distributions==2)
    {family1 = estimate_support_of_Y_and_pair_with_family(Y, family_of_distributions = family_of_distributions)}else{family1 = family_of_distributions}
    r1= as.numeric(estimate_epsilon(Y, X_, family = family1))
    
    #Which independence test should we choose? hoeffding.D.test or dhsic.test? hoeffding.D.test does not work well for discrete variables, but is faster
    z1=dhsic.test(data.frame(r2, X1), method = 'gamma')$p.value
    z2=dhsic.test(data.frame(r1, X2), method = 'gamma')$p.value
    
    
    
    return(create_bivariate_output_from_p_values(indep, z1, z2, family1, family2))
  }
  
  
  exact_CPCM_graph_estimate <- function(X, family_of_distributions = 1) {
    d <- ncol(X)
    stopifnot(d <= 4)
    
    if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(d))
    nodes <- colnames(X)
    
    generate_all_dags <- function(nodes) {
      d <- length(nodes)
      all_dags <- list()
      all_orders <- gtools::permutations(d, d)
      seen <- new.env(hash = TRUE)
      
      for (order in 1:nrow(all_orders)) {
        perm <- all_orders[order, ]
        
        for (edge_bits in 0:(2^choose(d, 2) - 1)) {
          bits <- as.integer(intToBits(edge_bits))[1:choose(d, 2)]
          mat <- matrix(0, d, d)
          k <- 1
          for (i in 2:d) {
            for (j in 1:(i - 1)) {
              mat[perm[i], perm[j]] <- bits[k]
              k <- k + 1
            }
          }
          
          key <- paste(mat, collapse = "")
          if (!exists(key, envir = seen)) {
            dag <- bnlearn::empty.graph(nodes)
            bnlearn::amat(dag) <- mat
            all_dags[[length(all_dags) + 1]] <- dag
            assign(key, TRUE, envir = seen)
          }
        }
      }
      return(all_dags)
    }
    
    dags <- generate_all_dags(nodes)
    p_values <- numeric(length(dags))
    p_values_adjusted <- numeric(length(dags))
    
    score_dag <- function(dag) {
      residuals <- matrix(0, nrow = nrow(X), ncol = d)
      colnames(residuals) <- nodes
      for (j in seq_len(d)) {
        pa <- parents(dag, nodes[j])
        if (length(pa) == 0) {
          residuals[, j] <- X[, j]
        } else {
          fam <- if (family_of_distributions %in% c(1, 2))
            estimate_support_of_Y_and_pair_with_family(X[, j], family_of_distributions)
          else family_of_distributions
          residuals[, j] <- estimate_epsilon(Y = X[, j], X = X[, pa, drop = FALSE], family = fam)
        }
      }
      p_val <- dhsic.test(as.data.frame(residuals), method = "gamma")$p.value
      score <- -log(p_val) + lambda * nrow(arcs(dag))
      return(list(score = score, p_value = p_val))
    }
    
    for (i in seq_along(dags)) {
      res <- score_dag(dags[[i]])
      p_values[i] <- res$p_value
      p_values_adjusted[i] <- res$score
      if(quiet ==FALSE) {
        cat("DAG:", i, "| Score:", round(res$score, 4), "| p-value:", round(res$p_value, 4), "\n")
      }
    }
    
    best_index <- which.min(p_values_adjusted)
    M_best <- dags[[best_index]]
    
    plausible <- which(p_values > 0.01)
    plausible_sorted <- plausible[order(p_values_adjusted[plausible])]
    M_plausible <- dags[plausible_sorted]
    
    return(list(
      plausible = M_plausible,
      p_values = p_values,
      p_values_adjusted = p_values_adjusted,
      result = M_best
    ))
  }
  
  
  
  greedy_CPCM_graph_estimate <- function(X, family_of_distributions = 1, use_RESIT_ordering = TRUE) {
    nodes <- colnames(X)
    d <- ncol(X)
    best_plausibility = 0
    # === Score function ===
    score_function <- function(X, dag) {
      epsilon <- X
      for (i in seq_along(nodes)) {
        node <- nodes[i]
        pa <- bnlearn::parents(dag, node)
        if (length(pa) > 0) {
          fam <- if (family_of_distributions %in% c(1,2)) 
            estimate_support_of_Y_and_pair_with_family(X[, node], family_of_distributions)
          else family_of_distributions
          epsilon[, node] <- estimate_epsilon(Y = X[, node], X = X[, pa, drop=FALSE], family = fam)
        }
      }
      dhsic_test = dhsic.test(epsilon, method = "gamma")
      return(data.frame(score =  -log(dhsic_test$p.value) + lambda * nrow(arcs(dag)),
                        p_value =  dhsic_test$p.value ))
    }
    
    arc_exists <- function(dag, from, to) {
      arcs_list <- arcs(dag)
      any(arcs_list[, 1] == from & arcs_list[, 2] == to)
    }
    
    multi_dhsic <- function(epsilon, X_parents) {
      min_p_val <- 1
      for (pj in colnames(X_parents)) {
        p_test <- dhsic.test(data.frame(epsilon, X_parents[[pj]]), method = "gamma")$p.value
        min_p_val <- min(min_p_val, p_test)
      }
      return(min_p_val)
    }
    
    # === Phase 1 ===
    if (use_RESIT_ordering) {
      # RESIT-style ordering
      ordering <- c()
      remaining_nodes <- nodes
      while (length(remaining_nodes) > 0) {
        p_vals <- c()
        for (target in remaining_nodes) {
          predictors <- setdiff(remaining_nodes, target)
          if (length(predictors) == 0) {
            p_vals[target] <- 1
            next
          }
          Y <- X[, target]
          X_parents <- X[, predictors, drop = FALSE]
          fam <- if (family_of_distributions %in% c(1,2)){ 
            estimate_support_of_Y_and_pair_with_family(Y, family_of_distributions)
          }else family_of_distributions
          epsilon <- estimate_epsilon(Y = Y, X = X_parents, family = fam)
          p_vals[target] <- multi_dhsic(epsilon, X_parents)
        }
        chosen <- names(which.max(p_vals))
        if (!quiet) cat("Removing:", chosen, "| p-values:", paste0(seq_along(p_vals), ": ", signif(p_vals, 4), collapse = ", "), "\n")
        ordering <- c(chosen, ordering)
        remaining_nodes <- setdiff(remaining_nodes, chosen)
      }
      
      # === Phase 2: Greedy backward pruning ===
      # Start from full DAG consistent with ordering
      current_dag <- bnlearn::empty.graph(nodes)
      for (i in seq_along(ordering)) {
        child <- ordering[i]
        if (i == 1) next
        parents <- ordering[1:(i - 1)]
        for (p in parents) {
          if (p != child) current_dag <- bnlearn::set.arc(current_dag, from = p, to = child)
        }
      }
      score = score_function(X, current_dag)
      best_score <- score$score
      best_plausibility = max( best_plausibility, score$p_value )
      if (!quiet) {
        arc_list <- apply(arcs(current_dag), 1, function(row) paste0(row[1], "->", row[2]))
        cat("Initial full DAG score =", round(best_score, 4), "| arcs:", paste(arc_list, collapse = ", "), "\n")
      }
      
      # Greedy backward pruning
      improved <- TRUE
      while (improved) {
        improved <- FALSE
        arc_list <- arcs(current_dag)
        best_temp_score <- best_score
        best_temp_dag <- current_dag
        best_arc_to_remove <- NULL
        
        for (k in seq_len(nrow(arc_list))) {
          from <- arc_list[k, 1]
          to <- arc_list[k, 2]
          temp_dag <- drop.arc(current_dag, from, to)
          score <- score_function(X, temp_dag)
          temp_score <- score$score
          best_plausibility <- max(best_plausibility, score$p_value)
          
          if (!quiet) {
            cat("Try removing", from, "->", to, "| score =", round(temp_score, 4), "\n")
          }
          
          if (temp_score < best_temp_score) {
            best_temp_score <- temp_score
            best_temp_dag <- temp_dag
            best_arc_to_remove <- c(from, to)
          }
        }
        
        if (!is.null(best_arc_to_remove)) {
          current_dag <- best_temp_dag
          best_score <- best_temp_score
          improved <- TRUE
          if (!quiet) {
            cat("Removed", best_arc_to_remove[1], "->", best_arc_to_remove[2], "| new best score =", round(best_score, 4), "\n")
          }
        }
      }
      
      
    } else {
      # Greedy edge addition
      current_dag <- empty.graph(nodes)
      score = score_function(X, current_dag)
      best_score <- score$score
      best_plausibility = max( best_plausibility, score$p_value )
      if (!quiet) cat("Initial score (empty DAG):", best_score, "\n")
      
      improved <- TRUE
      while (improved) {
        improved <- FALSE
        best_candidate_score <- Inf
        best_candidate_dag <- NULL
        
        for (i in 1:d) {
          for (j in 1:d) {
            if (i != j && !arc_exists(current_dag, nodes[i], nodes[j])) {
              temp_dag <- tryCatch(set.arc(current_dag, nodes[i], nodes[j], check.cycles = TRUE),
                                   error = function(e) NULL)
              if (is.null(temp_dag)) next
              score= score_function(X, temp_dag)
              temp_score <- score$score
              best_plausibility = max(best_plausibility, score$p_value)
              if (!quiet) {
                arc_list <- apply(arcs(temp_dag), 1, function(row) paste0(row[1], "->", row[2]))
                cat("Trying DAG: Score =", round(temp_score, 4), "| arcs:", paste(arc_list, collapse = ", "), "\n")
              }
              if (temp_score < best_candidate_score) {
                best_candidate_score <- temp_score
                best_candidate_dag <- temp_dag
              }
            }
          }
        }
        
        if (!is.null(best_candidate_dag) && best_candidate_score < best_score) {
          current_dag <- best_candidate_dag
          best_score <- best_candidate_score
          improved <- TRUE
          if (!quiet) cat("Updated DAG with score:", best_score, "\n")
        }
      }
      
      # === Phase 2: Prune phase (same as fast version)  === #
      for (node in nodes) {
        pa <- parents(current_dag, node)
        for (p in pa) {
          reduced_parents <- setdiff(pa, p)
          fam <- if (family_of_distributions %in% c(1,2)) 
            estimate_support_of_Y_and_pair_with_family(X[, node], family_of_distributions)
          else family_of_distributions
          
          eps <- if (length(reduced_parents) == 0) X[, node] else
            estimate_epsilon(Y = X[, node], X = X[, reduced_parents, drop=FALSE], family = fam)
          
          p_val <- dhsic.test(data.frame(eps, X[, p]), method = "gamma")$p.value
          if (p_val > 0.05) {
            current_dag <- drop.arc(current_dag, from = p, to = node)
            if (!quiet) cat("Pruned edge:", p, "→", node, "| p =", round(p_val, 4), "\n")
          }
        }
      }
    }
    
    return(list(dag  =  current_dag, best_plausibility= best_plausibility))
  }
  
  
  
  RESIT_CPCM_graph_estimate <- function(X, family_of_distributions = 1) {
    nodes <- colnames(X)
    d <- ncol(X)
    ordering <- c()
    remaining_nodes <- nodes
    
    multi_dhsic = function(epsilon, X_parents){
      min_p_val <- 1
      for (pj in colnames(X_parents)) {
        p_test <- dhsic.test(data.frame(epsilon, X_parents[[pj]]), method = "gamma")$p.value
        min_p_val <- min(min_p_val, p_test)
      }
      return(min_p_val)
    }
    ### === Phase 1: ordering === ###
    
    # Initialize adjacency matrix
    adj_mat <- matrix(0, nrow = d, ncol = d)
    colnames(adj_mat) <- rownames(adj_mat) <- nodes
    
    while (length(remaining_nodes) > 0) {
      p_vals <- c()
      
      for (target in remaining_nodes) {
        predictors <- setdiff(remaining_nodes, target)
        
        if (length(predictors) == 0) {
          p_vals[target] <- 1  # No parents -> independent
          next
        }
        
        Y <- X[, target]
        X_parents <- X[, predictors, drop = FALSE]
        
        # Estimate residuals
        if (family_of_distributions == 1 || family_of_distributions == 2) {
          fam <- estimate_support_of_Y_and_pair_with_family(Y, family_of_distributions)
        } else {
          fam <- family_of_distributions
        }
        epsilon <- estimate_epsilon(Y = Y, X = X_parents, family = fam)
        p_vals[target] <- multi_dhsic(epsilon, X_parents)
      }
      
      # Choose most independent node (largest p-value)
      chosen <- names(which.max(p_vals))
      if (quiet == FALSE) cat("Removing:", chosen, "| p-value:", p_vals[chosen], "\n")
      
      ordering <- c(chosen, ordering)
      remaining_nodes <- setdiff(remaining_nodes, chosen)
    }
    
    
    
    
    ### === Phase 2: build DAG with pruning === ###
    dag <- empty.graph(nodes)
    
    for (i in seq_along(ordering)) {
      node <- ordering[i]
      Y <- X[, node]
      if (family_of_distributions == 1 || family_of_distributions == 2) {
        fam <- estimate_support_of_Y_and_pair_with_family(Y, family_of_distributions)
      } else {
        fam <- family_of_distributions
      }
      
      parents <- if (i == 1) character(0) else ordering[1:(i - 1)]
      
      for (p in parents) {
        other_parents <- setdiff(parents, p)
        
        if (length(other_parents) == 0){ 
          test_df <- data.frame(Y, X[, p])
          p_val <- dhsic.test(test_df, method = "gamma")$p.value
        }else{
          X_parents = X[, parents, drop = FALSE]
          X_other = X[, other_parents, drop = FALSE]
          epsilon <- estimate_epsilon(Y = Y, X = X_other, family = fam)
          p_val <- multi_dhsic(epsilon, X_other)
        }
        
        if (p_val < 0.05) {
          dag <- set.arc(dag, from = p, to = node)
        } else if (!quiet) {
          cat("Pruned edge:", p, "→", node, "| p =", round(p_val, 4), "\n")
        }}
    }
    
    return(dag)
  }
  
  
  
  test_residual_independence <- function(X, best_dag, family_of_distributions = 1) {
    nodes <- colnames(X)
    d <- ncol(X)
    residual_matrix <- X  # initialize
    
    for (i in seq_along(nodes)) {
      node <- nodes[i]
      parents_i <- parents(best_dag, node)
      
      if (length(parents_i) == 0) {
        residual_matrix[, node] <- X[, node]  # No parents → leave as is
      } else {
        Y <- X[, node]
        X_parents <- X[, parents_i, drop = FALSE]
        
        if (family_of_distributions == 1 || family_of_distributions == 2) {
          fam <- estimate_support_of_Y_and_pair_with_family(Y, family_of_distributions)
        } else {
          fam <- family_of_distributions
        }
        residual_matrix[, node] <- estimate_epsilon(Y = Y, X = X_parents, family = fam)
      }
    }
    p_value = dhsic.test(as.data.frame(residual_matrix), method = "gamma")$p.value
    # Final joint independence test
    return(p_value)
  }
  
  
  
  if(ncol(X)==2) {
    if(family_of_distributions!='Sequential choice'){
      return(bivariate_CPCM_graph_estimate(X, family_of_distributions=family_of_distributions))}
    
    if(family_of_distributions=='Sequential choice'){
      alpha_level = 0.01
      
      result1 = bivariate_CPCM_graph_estimate(X, family_of_distributions = 1)
      if(max(as.numeric(result1[[1]][1:3]))>alpha_level){return(result1)}else{
        result2 = bivariate_CPCM_graph_estimate(X, family_of_distributions = 2)
        if(max(as.numeric(result2[[1]][1:3]))>alpha_level){return(result2)}else{
          result3 = result2
          result3[[1]][1] = max(as.numeric(result1[[1]][1]), as.numeric(result2[[1]][1]))
          result3[[1]][2] = max(as.numeric(result1[[1]][2]), as.numeric(result2[[1]][2]))
          result3[[1]][3] = max(as.numeric(result1[[1]][3]), as.numeric(result2[[1]][3]))
          
          family1 = sub(";.*", "", result2[[1]][7])
          family2 = sub(".*?;", "", result2[[1]][7])
          if(as.numeric(result1[[1]][2])>as.numeric(result2[[1]][2])){family1 = sub(";.*", "", result1[[1]][7])}
          if(as.numeric(result1[[1]][3])>as.numeric(result2[[1]][3])){family2 = sub(".*?;", "", result1[[1]][7])}
          
          result=create_bivariate_output_from_p_values(as.numeric(result3[[1]][1]), 
                                                       as.numeric(result3[[1]][2]), 
                                                       as.numeric(result3[[1]][3]),
                                                       family1, family2)
          return(result)
        }
      }
    } 
  }
  
  if(greedy_method == 'exact') {
    if (ncol(X) == 3 || ncol(X) == 4) {
      
      if (family_of_distributions != 'Sequential choice') {
        return(exact_CPCM_graph_estimate(X, family_of_distributions = family_of_distributions)$result)
        
      } else {
        alpha_level <- 0.01
        
        result1 <- exact_CPCM_graph_estimate(X, family_of_distributions = 1)
        if (max(result1$p_values) > alpha_level) {return(result1$result)}
        
        result2 <- trivariate_exact_CPCM_graph_estimate(X, family_of_distributions = 2)
        if (max(result2$p_values) > alpha_level) {return(result2$result)}
        
        # Combine the two results conservatively
        result3 <- result2  # base structure
        result3$p_values <- pmax(result1$p_values, result2$p_values)
        result3$p_values_adjusted <- pmin(result1$p_values_adjusted, result2$p_values_adjusted)
        result3$result <- if (max(result1$p_values) > max(result2$p_values)) result1$result else result2$result
        
        return(result3$result)
      }
      
    } else {
      cat("Warning: Exact estimation is not implemented for d > 4. If d > 4, using edge_greedy estimation instead.\n")
      greedy_method <- 'edge_greedy'
    }
  }
  
  if(greedy_method=='RESIT'& ncol(X)>2) {
    if(family_of_distributions!='Sequential choice'){ 
      return(RESIT_CPCM_graph_estimate(X, family_of_distributions=family_of_distributions))
    }
    if(family_of_distributions=='Sequential choice'){
      alpha_level = 0.01
      result1 = RESIT_CPCM_graph_estimate(X, family_of_distributions = 1)
      plausibility1 = test_residual_independence(X, result1, family_of_distributions = 1)
      if(quiet==FALSE){cat('Plausibility for family1:', plausibility1, '\n')}
      if(plausibility1>= alpha_level) {return(result1)}else{
        result2 = RESIT_CPCM_graph_estimate(X, family_of_distributions = 2)
        return(result2)
      }
    }
  }
  
  if ((greedy_method == 'RESIT_greedy' || greedy_method == 'edge_greedy') && ncol(X) > 2) {
    if(greedy_method == 'RESIT_greedy'){use_RESIT_ordering=TRUE}else{use_RESIT_ordering=FALSE}
    
    if(family_of_distributions!='Sequential choice'){ 
      return(greedy_CPCM_graph_estimate(X, family_of_distributions=family_of_distributions, 
                                        use_RESIT_ordering = use_RESIT_ordering)[[1]])
    }
    if(family_of_distributions=='Sequential choice'){
      alpha_level = 0.01
      
      result1 = greedy_CPCM_graph_estimate(X, family_of_distributions = 1, 
                                           use_RESIT_ordering = use_RESIT_ordering)
      plausibility1 = result1$best_plausibility
      if(quiet==FALSE){cat('Plausibility for family1:', plausibility1, '\n')}
      if(plausibility1>= alpha_level) {return(result1$dag)}else{
        result2 = greedy_CPCM_graph_estimate(X, family_of_distributions = 2, 
                                             use_RESIT_ordering = use_RESIT_ordering)
        return(result2$dag)
      }
    }
  }        
}


