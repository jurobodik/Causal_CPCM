#Data generators

library(MASS)
library(EnvStats)


# First five generators are taken from 
#N. Tagasovska, V. Chavez-Demoulin, and T. Vatter. Distinguishing cause from effect using
#quantiles: Bivariate quantile causal discovery. PMLR 2020
# Parts of these generatos are inspired by the supplementary material for 
# J. Peters, J. Mooij, D. Janzing, B. Sch\"olkopf: 
# "Causal Discovery with Continuous Additive Noise Models", JMLR 2014
# Functions are modified to adapt for our experiments, so we include their copyright notice.
# Copyright (c) 2013  Jonas Peters [peters@stat.math.ethz.ch]
# All rights reserved.
# 

sample_AN <- function(n){
  sample_generative_model(n, 2, 1)
}

sample_ANs = function(n){
  sample_generative_model(n, 1, 1)
}

sample_MN_u = function(n){
  sample_generative_model(n, 1, 2)
}

sample_LS <- function(n){
  sample_generative_model(n, 2, 3)
}

sample_LSs = function(n){
  sample_generative_model(n, 1, 3)
}

sample_Poisson = function(n){
  CPCM_generative_model(n, 'Poisson')
}

sample_pareto = function(n){
  CPCM_generative_model(n, 'Pareto')
}


sample_generative_model<- function(n, mech_type = 1, noise_type = 1){
  # mech_type - 1:injective (sigmoid), 2: non-injective
  # noise_type - 1:additive, 2:multiplicative, 3: ls

   ran <- rnorm(n)
   noise_exp <- 1
   noise_var <- runif(n, 1, 2)
   noisetmp <- (sqrt(noise_var) * abs(ran))^(noise_exp) * sign(ran)
   x_pa <- noisetmp
   
   noise_var_ch <- runif(n, 1, 2)
  
  if(mech_type == 1){
    a.sig <- runif(n=1, min=-2, max=2)
    bern <- rbinom(1,1,0.5)
    b.sig <- bern*runif(n=1, min=0.5, max=2) + (1-bern)*runif(n=1, min=-2, max=-0.5)
    c.sig <- rexp(n=1,rate=4)+1
    x_child <- c.sig*(b.sig*(x_pa+a.sig))/(1+abs(b.sig*(x_pa +a.sig)))
  } else if(mech_type == 2){
    kern_pa <- gauss_kernel(x_pa, 1, 1)
    x_child <- mvrnorm(1, rep(0, n), kern_pa)
  }
  ran <- rnorm(n)
  noisetmp <- (0.2 * sqrt(noise_var_ch) * abs(ran)) ^ (noise_exp) * sign(ran)
  
  
  if (noise_type == 1) {
    x_child <- x_child + noisetmp
  } else if (noise_type == 2) {
    x_child <- x_child * runif(n)
  } else if (noise_type == 3) {
    ran <- rnorm(n)
    x_child <- x_child + (x_child - min(x_child))*noisetmp
  } else if (noise_type == 4) {
    ran <- rnorm(n)
    sd = (x_child - min(x_child))
    x_child <- (0.2 * sqrt(sd) * abs(ran)) ^ (noise_exp) * sign(ran)
  } else {
    print("model type not implemented")
  }
  
  cbind(x_pa, x_child)
}


gauss_kernel <- function(x, sigmay, sigmax) {
  
  if (is.matrix(x) == FALSE)
    x <- as.matrix(x)
  
  n <- nrow(x)
  xnorm <- as.matrix(dist(x, method = "euclidean", diag = TRUE, upper = TRUE))
  
  sigmay * exp(-xnorm^2/(2*sigmax^2))
}


CPCM_generative_model = function(n, family){
  
  
  random_function <- function(draw_function = FALSE, max = 10){
    cov_matrix <- function(x, kernel_fn, ...) {
      outer(x, x, function(a, b) kernel_fn(a, b, ...))
    }
    # given x coordinates, take N draws from kernel function at those points
    draw_samples <- function(x, N, kernel_fn, ...) {
      Y <- matrix(NA, nrow = length(x), ncol = N)
      for (n in 1:N) {
        K <- cov_matrix(x, kernel_fn, ...)
        Y[, n] <- mvrnorm(1, mu = rep(0, times = length(x)), Sigma = K)
      }
      Y
    }
    se_kernel <- function(x, y, sigma = 1, length = 1) {
      sigma^2 * exp(- (x - y)^2 / (2 * length^2))
    }
    sekvencia <- seq(0, max, length.out = 201)  # x-coordinates
    Y <- draw_samples(sekvencia, 1, kernel_fn = se_kernel, length = 0.2)
    Y =Y-min(Y)+0.5
    if (draw_function==TRUE) {
      plot(range(sekvencia), range(Y), xlab = "sekvencia", ylab = "y", type = "n",
           main = "SE kernel, length = 0.2")
      for (n in 1:N) {
        lines(sekvencia, Y[, n])
      }
    }
    
    return(Y)
  }
  
  random_theta <- function(x,f, max){
    sekvencia <- seq(0, max, length.out = 201) 
    theta= c()
    for (i in 1:length(x)) {
      theta = c(theta, f[which.min(abs(x[i] - sekvencia))])
    }
    return(theta)
  }
  
  
  
  X1 = rnorm(n, 0, sqrt(2))
  X2=c();
  
  
  
  f = random_function(max = ceiling(max(X1)))
  theta = random_theta(X1, f,max = ceiling(max(X1)))
  
  if(family == 'Poisson'){for (i in 1:n) {X2 = c(X2, rpois(1,theta[i]))}}
  if(family == 'Pareto'){for (i in 1:n) {X2 = c(X2, rpareto(1, theta[i], 1))}}
      
  return(X = data.frame(X1, X2) )
  
  
  
}





