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
  - `DistributionName`: Uses a model CPCM(F) where F=DistributionName. Choices are written below.
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
Good approach is to first try family1 and if all graphs are not plausible, use family2. However, if we do not want to compute everything twice, a good rule of thumb that works on the considered data is the following: 

- Use `family_of_distributions = 1` if `n <= 1000`.
- Use `family_of_distributions = 2` if `n > 1000`,
  but again, this choice should depend on the dataset complexity.

### Output 
Output format depends on the dimension of X. 

In the bivariate case, the output is in the following format:
                                     
- p-value 1-->2                       0.909445     (represents the p-value of the independence test between X1 and hat{epsilon_2})
- p-value 2-->1                       0.000015     (represents the p-value of the independence test between X2 and hat{epsilon_1})
- Score-based graph estimate           1 --> 2     (Score-based graph estimate represents the final estimate of the graph)
- Testing estimate                     1 --> 2     (Testing estimate represents the final estimate of the graph using Algorithm 1. Note that five different outputs are possible) 
- Families used              Gaussian;Gaussian     (first is the family used for X1 and the second is the family used for X2 (that is, family with the same support)

In the case d=3, the output is in the following format:
- $result                       Final score-based estimation of G; in particular, its adjacency matrix
- $plausible                    List of all plausible graphs
- $p_values                     p-values of the independence test between (hat{epsilon_1}, hat{epsilon_2}, hat{epsilon_3}) for all graphs
- $p_values_adjusted            scores of all graphs using score=-log(p_values)+ lambda*(number of edges)

In the case d>3, the output is the score-based estimation of G via greedy search: 
- DAG in the format of bnlearn environment (see https://cran.r-project.org/web/packages/bnlearn/index.html)


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

#>                                     Results
#>p-value 1-->2                       0.909445
#>p-value 2-->1                       0.000015
#>Score-based graph estimate           1 --> 2
#>Testing estimate                     1 --> 2
#>Families used              Gaussian;Gaussian
```
