#Code for the Simulations 3 with misspecified F
#We ran each element of the resulting table one by one and rewrote into the table by hand
#Elements of the Table can be recreated as follows:
#Choose family and choose function of theta



family= "Gamma with fixed scale"
#family= "Gamma"
#family= "Pareto"
#family= "Gaussian with fixed sigma"
#family= "Gaussian"


function_for_theta<-function(x) return(x)
#function_for_theta<-function(x) return(x^2+1)
#function_for_theta<-function(x) return(exp(x)/2)
#Random theta is implemented via function "random_theta", you only have to erase # in the code at lines 71, 72 and erase the next line



#Run everything below, you will get one entry of the table
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


#set.seed(42)
n=500
result = c()
for (k in 1:100) {
  
  X1 = c()
  while (length(X1)<n) {
    value = rnorm(1,2,1)
    if (value>0) {X1 = c(X1, value)}
  }
  # f = random_function(max = ceiling(max(X1)))
  # theta = random_theta(X1, f,max = ceiling(max(X1)))
  
  
  theta =function_for_theta(X1) 
  
  
  X2 = c()
  for (i in 1:n) {
    X2 = c(X2, rexp(1, theta[i]))
  }
  
  X=data.frame(X1, X2)
  if(family == "Pareto"){X = X+1} #This is here because Pareto is implemented only with support (1, inf) and hence we just move it by 1
  graph = CPCM_graph_estimate(X, family_of_distributions = family)
  result = c(result, 2-as.numeric(substring(graph, 1, 1)) ) #1 if correct 1-->2, 0 otherwise
  cat(k, 2-as.numeric(substring(graph, 1, 1)), "\n")
}

sum(result)

































