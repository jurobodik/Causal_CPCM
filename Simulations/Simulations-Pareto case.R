#####################     Simulations about Pareto case and Consequence 1    ################################
#Running the code below, you can recreate the results and graphs from Section 6.2.
#You need to first upload the function 'CPCM_graph_estimate', that is located in the main file
#Lines 1-200 contain the simulations for the first plot
#Lines 200-300 contain the simulations for the second plot
#Lines 300-400 contain the simulations for the third plot

library(ggplot2)
library(ggpattern)
library(gridExtra)
library(EnvStats) #only for rpareto()



set.seed(1)
#The following function generates the random sample following the distribution from CPCM(F) where 
#p_{X_1}(x) \propto \frac{1}{ [\log(x)+1] x^{2} },    and     \theta(x) = x^\alpha\log(x)+1.
#We use MCMC approach for generating random variable X given a density function f. See https://stats.stackexchange.com/questions/86909/how-to-generate-random-variables-from-a-defined-density-via-r
# Function to generate data
generate_our_data <- function(n = 1000, alpha) {
  m <- 1000 * n
  f <- function(x) 1.676875028178700 / (x^2 * (log(x) + 1)) #density of X1
  M <- 1.7
  u <- runif(m, 1, 50)
  vysledok <- integer(0)
  
  for (i in 1:m) {
    y <- f(u[i]) / M
    if (sample(c(0, 1), prob = c(y, 1 - y), size = 1) == 0) {
      vysledok <- c(vysledok, i)
    }
  }
  
  X <- u[vysledok][1:n]
  
  Y <- numeric(n)
  for (i in 1:n) {
    theta <- X[i]^alpha * log(X[i]) + 2
    y <- rpareto(1, 1, theta)
    Y[i] <- y
  }
  
  return(data.frame(X = X, Y = Y))
}

# Example of usage:
X = generate_our_data(n=500, alpha=0)
plot(X)
CPCM_graph_estimate(X, family_of_distributions = 'Pareto') 

# Repeat the process multiple times for different alphas
number_of_repetitions <- 100
n <- 500
final_result <- matrix(NA, nrow = number_of_repetitions, ncol = length(c(-2, 1, 0, 1, 2)))

for (j in 1:length(c(-2, 1, 0, 1, 2))) {
  alpha <- c(-2, 1, 0, 1, 2)[j]
  result <- numeric(number_of_repetitions)
  
  for (i in 1:number_of_repetitions) {
    X <- generate_our_data(n = n, alpha)
    graph <- CPCM_graph_estimate(X, family_of_distributions = "Pareto")
    result[i] <- graph[5,]
    cat("Time remaining: alpha =", alpha, "and repetitions remaining =", number_of_repetitions - i, "\n")
  }
  
  final_result[, j] <- result
}

# Draw the results on a bar plot
percentages <- matrix(0, nrow = 5, ncol = 5)

for (j in 1:5) {
  percentages[j, 1] <- sum(final_result[, j] == "1 --> 2")
  percentages[j, 2] <- sum(final_result[, j] == "2 --> 1")
  percentages[j, 3] <- sum(final_result[, j] == "Empty graph")
  percentages[j, 4] <- sum(final_result[, j] == "Unidentifiable  (both directions are plausible)")
  percentages[j, 5] <- sum(final_result[, j] == "Assumptions not fulfilled (both directions are not plausible)")
}

percentage1 <- percentages[, 1]
percentage2 <- percentages[, 2]
percentage3 <- percentages[, 3]
percentage4 <- percentages[, 4]
percentage5 <- percentages[, 5]


# Create the data with five categories
data <- data.frame(
  Category = factor(rep(c("-2", "-1", "0", "1", "2"), each = 5), levels = c("-2", "-1", "0", "1", "2")),
  Segment = rep(c("1 --> 2", "2 --> 1", "Empty graph", "Both directions are plausible", "Neither direction is plausible"), times = 5),
  Percentage = c(percentage1,     
                 percentage3,     
                 percentage2,    
                 percentage4, 
                 percentage5) 
)

# Create the bar plot
ggplot(data, aes(x = Category, y = Percentage, fill = Segment)) +
  geom_bar(stat = "identity", position = position_stack(), width = 0.3) +
  theme_minimal(base_size = 15) +
  labs(title = "Identifiability of Pareto CPCM model", x = "alpha", y = "Percentage") +
  scale_fill_manual(values = c("1 --> 2" = "#1f78b4", "2 --> 1" = "#ff7f00", "Both directions are plausible" = "#e31a1c", "Neither direction is plausible" = "#6a3d9a", "Empty graph" = "#33a02c")) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank(),
        legend.position = "bottom")



























































































######Graph number two
number_of_repetitions = 100

alpha = 0
result1 = c()
result2 = c()

for (i in 1:number_of_repetitions) {
  X = generate_our_data(n=500, alpha)
  graph = CPCM_graph_estimate(X, family_of_distributions = "Pareto")
  result1 = c(result1, graph[1,])
  result2 = c(result2, graph[2,])
  cat("Time remaining: ", number_of_repetitions-i, "\n")
}


number_of_repetitions = 100

alpha = 2
result3 = c()
result4 = c()

for (i in 1:number_of_repetitions) {
  X = generate_our_data(n=500, alpha)
  graph = CPCM_graph_estimate(X, family_of_distributions = "Pareto")
  result3 = c(result3, graph[1,])
  result4 = c(result4, graph[2,])
  cat("Time remaining: ", number_of_repetitions-i, "\n")
}



result1 = result2[result1!='-']
result2 = result2[result2!='-']


result1=as.numeric(result1)
result2=as.numeric(result2)
result3=as.numeric(result3)
result4=as.numeric(result4)





binwidth=0.1

# Create the first pair of overlapping histograms
plot1 <- ggplot() +
  geom_histogram(aes(x = result1), binwidth = binwidth, fill = "blue", alpha = 0.2) +
  geom_histogram(aes(x = result2), binwidth = binwidth, fill = "red", alpha = 0.2) +
  labs(title = "Overlapping Histograms 1", x = "p-value", y = "Frequency") +  
  theme_minimal()

# Create the second pair of overlapping histograms
plot2 <- ggplot() +
  geom_histogram(aes(x = result3), binwidth = binwidth, fill = "green", alpha = 0.2) +
  geom_histogram(aes(x = result4), binwidth = binwidth, fill = "purple", alpha = 0.2) +
  labs(title = "Overlapping Histograms 2", x = "Value", y = "Frequency") +
  theme_minimal()

# Arrange the plots vertically
grid.arrange(plot1, plot2, ncol = 1)





# Combine data into data frames
df1 <- data.frame(value = c(result1, result2), direction = c(rep("X -> Y", 100), rep("Y -> X", 100)))
df2 <- data.frame(value = c(result3, result4), direction = c(rep("X -> Y", 100), rep("Y -> X", 100)))

binwidth <- 0.1

# Create the first pair of overlapping histograms with density lines
plot1 <- ggplot(df1, aes(x = value, fill = direction)) +
  geom_histogram(position = "identity", alpha = 0.1, binwidth = binwidth) +
  labs(title = "P-value distribution: case alpha = 0", x = "p-value", y = "Frequency") +
  scale_fill_manual(values = c("X -> Y" = "blue", "Y -> X" = "red")) +
  scale_color_manual(values = c("X -> Y" = "blue", "Y -> X" = "red")) +
  theme_minimal() +
  theme(legend.position = "top")

# Create the second pair of overlapping histograms with density lines
plot2 <- ggplot(df2, aes(x = value, fill = direction)) +
  geom_histogram(position = "identity", alpha = 0.3, binwidth = binwidth) +
  labs(title = "P-value distribution: case alpha = 2", x = "p-value", y = "Frequency") +
  scale_fill_manual(values = c("X -> Y" = "green", "Y -> X" = "purple")) +
  scale_color_manual(values = c("X -> Y" = "green", "Y -> X" = "purple")) +
  theme_minimal() +
  theme(legend.position = "top")

# Arrange the plots vertically
grid.arrange(plot1, plot2, ncol = 1)






######Graph number three
sample_sizes = c(50, 100,200,  300, 400, 500, 600, 800, 1000)
number_of_repetitions = 200

stupid_transmission <- function(x){ #return 1 if 1-->2, return 0 if 2-->1
  if(x=="Empty graph")return(0)else{return(2-as.numeric(substr(x[1],1,1)))}
}

alpha = 1
final_result = c()
for (n in sample_sizes) {
  result = c()
  for (i in 1:number_of_repetitions) {
    X = generate_our_data(n=n, alpha)
    graph = CPCM_graph_estimate(X, family_of_distributions = "Pareto")
    result = c(result, graph[3,])
    cat("Time remaining: n = ", n, '  and  number_of_repetitions = ', number_of_repetitions-i, "\n")
  }
  
  final_result = cbind(final_result, result)
  
}



results_1 = rep(0, length(final_result[1,]))
for (j in 1:length(final_result[1,])) {
  for (i in 1:length(final_result[,1])) {
    results_1[j] = results_1[j]+stupid_transmission(final_result[i,j])}
  results_1[j] = results_1[j]/length(final_result[,1])
}



alpha = 2 
final_result = c()
for (n in sample_sizes) {
  result = c()
  for (i in 1:number_of_repetitions) {
    X = generate_our_data(n=n, alpha)
    graph = CPCM_graph_estimate(X, family_of_distributions = "Pareto")
    result = c(result, graph[3,])
    cat("Time remaining: n = ", n, '  and  number_of_repetitions = ', number_of_repetitions-i, "\n")
  }
  
  final_result = cbind(final_result, result)
  
}

results_2 = rep(0, length(final_result[1,]))
for (j in 1:length(final_result[1,])) {
  for (i in 1:length(final_result[,1])) {
    results_2[j] = results_2[j]+stupid_transmission(final_result[i,j])}
  results_2[j] = results_2[j]/length(final_result[,1])
}






# Combine into a data frame
data <- data.frame(
  SampleSize = rep(sample_sizes, 2),
  CorrectEstimationFraction = c(results_1, results_2),
  Line = factor(rep(c("alpha = 1", "alpha = 2"), each = length(sample_sizes)))
)

# Plot the data
ggplot(data, aes(x = SampleSize, y = CorrectEstimationFraction, color = Line)) +
  geom_line(size = 1) +
  labs(
    title = "Score based algorithm estimate",
    x = "Size of a Sample n",
    y = "Correct Estimation Fraction"
  ) +
  theme_minimal() +
  ylim(0.3, 1)



