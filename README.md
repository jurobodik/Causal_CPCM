# CPCM

This repository contains the code for estimating the causal graph using the CPCM (Conditionally parametric causal models) method. The main function `CPCM_graph_estimate()` implements this method and allows for flexible specification of the family of distributions.
The `CPCM_graph_estimate` function is designed to estimate causal relationships among variables using the CPCM method, providing insights into the underlying causal structure based on the specified family of distributions.

## Function Overview

The function `CPCM_graph_estimate(X, family_of_distributions, force_estimate)` estimates the causal graph based on the provided data `X` and the specified family of distributions.

### Parameters

- **X**: The input data matrix or data frame where each column represents a variable.
- **family_of_distributions**: Specifies the family of distributions to use:
  - `1`: Uses a joint model CPCM(F1...Fk) with 1-parameter distributions.
  - `2`: Uses a joint model CPCM(F1...Fk) with 2-parameter distributions.
- **force_estimate**: Optional parameter (default `FALSE`). If set to `TRUE`, an empty graph is not allowed.

### Available Distributions

For `family_of_distributions = 1`:
- `Gaussian with fixed sigma`
- `Pareto`
- `Poisson`

For `family_of_distributions = 2`:
- `Gaussian`
- `Gamma`
- `Negative_binomial`


If you want to use CPCM(F) model with given $F$, the choices for 'family_of_distributions' are the following: `Gaussian`, `Gaussian with fixed sigma`, `Pareto`, `Gamma`, `Gamma with fixed scale`, `Gumbel`, `Gumbel with fixed scale`, `Negative_binomial`, `Poisson`

### Rule of Thumb

- Use `family_of_distributions = 1` if `n <= 1000`.
- Use `family_of_distributions = 2` if `n > 1000`, but this choice should depend on the dataset complexity.

### Output 
In the bivariate case, the output is in the following format:
-                                    Results
-p-value 1-->2                       0.909445
-p-value 2-->1                       0.000015
-Score-based graph estimate           1 --> 2
-Testing estimate                     1 --> 2
-Families used              Gaussian;Gaussian`

-p-value 1-->2  represents the p-value of the independence test between X1 and hat{epsilon_2}
-p-value 2-->1  represents the p-value of the independence test between X2 and hat{epsilon_1}
-Score-based graph estimate represents the final estimate of the graph
-Testing estimate represents the final estimate of the graph using Algorithm 1. Note that five different options can happen:  1) $X_1\indep X_2$, 2) $X_1\to X_2$, 3) $X_2\to X_1$, 4) ``unidentifiable setup'' (both directions appear to be plausible) and 5) ``Assumptions not fulfilled'' (neither direction appear to be plausible). 
-Families used: first is the family used for X1 and the second is the family used for X2 (that is, family with the same support). 

## Example

```r
# Generate example data
n <- 500
X1 <- rnorm(n)
X2 <- numeric()
for (i in 1:n) {
  X2 <- c(X2, rnorm(1, X1[i], X1[i]^2 + 1))
}
X <- data.frame(X1, X2)
plot(X)

# Estimate causal graph using CPCM
CPCM_graph_estimate(X, family_of_distributions = 2)
CPCM_graph_estimate(X, family_of_distributions = 'Gaussian') #The same results, but forces the use of Gaussian F
```
