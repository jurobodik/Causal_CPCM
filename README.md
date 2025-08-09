# CPCM

**Conditionally Parametric Causal Models (CPCM)** is a flexible framework for discovering causal structures when conditional distributions belong to known parametric families. This repository implements the `CPCM_graph_estimate()` function in R for inferring causal graphs under such assumptions.

---

##  Overview

CPCM models extend the identifiability of causal direction by allowing the effect’s conditional distribution to vary in mean, variance, tail behavior, etc., according to the parent variables. This removes restrictive assumptions typical in additive-noise models. For the theoretical foundations and empirical evaluation of CPCM, see:

- *Identifiability of causal graphs under nonadditive conditionally parametric causal models*, by Juraj Bodik & Valérie Chavez‑Demoulin (2023, arXiv) (https://arxiv.org/abs/2303.15376)

---

##  Installation

```r
install.packages(c("mgcv", "dHSIC", "bnlearn", "MASS", "gamlss", "stringr", "dplyr"))
git clone https://github.com/jurobodik/Causal_CPCM.git
```

---

##  Usage

```r
source("CPCM_graph_estimate.R")

result <- CPCM_graph_estimate(
  X, 
  family_of_distributions = 1,
  greedy_method = "RESIT_greedy",
  lambda = 1,
  quiet = TRUE
)
```

---

##  Parameters

| Parameter                | Description |
|--------------------------|-------------|
| **X**                    | Data frame or matrix; variables as columns. |
| **family_of_distributions** | Model to use: <br>**Joint-family models**:<br>`1` → {Gaussian with fixed sigma, Poisson, Exponential, Pareto}<br>`2` → {Gaussian, Negative_binomial, Gamma, Pareto2}<br>`"Sequential choice"` → Automatically tries family 1, then 2 if 1 is unplausible.<br>**Single-family models**: `"Gaussian"`, `"Gaussian with fixed sigma"`, `"Pareto"`, `"Pareto2"`, `"Exponential"`, `"Gamma"`, `"Gamma with fixed scale"`, `"Gumbel"`, `"Gumbel with fixed scale"`, `"Poisson"`, `"Negative_binomial"` |
| **greedy_method**         | Greedy search algorithm for `d > 2`: `"exact"`, `"RESIT_greedy"` (default), `"edge_greedy"`, `"RESIT"`. |
| **lambda**                | Complexity penalty (per edge) in score-based search. |
| **quiet**                 | Set `FALSE` to see progress updates. |

---

##  Output

- **Bivariate (d = 2)**:
  -  `p-value for empty graph`
  - `p-value 1→2`,
  - `p-value 2→1`
  - `Score-based graph estimate`
  - `Testing estimate`
  - `Forced estimate`
  - `Families used`

- **Multivariate (d > 2)**:
  - Returns estimated DAG in the **bnlearn** format ([bnlearn on CRAN](https://cran.r-project.org/web/packages/bnlearn)).

---

##  Extending the Method

To add new families or estimation logic:

```r
estimate_epsilon(Y, X, family = "YOUR_DISTRIBUTION")
```

The implementation leverages `gam()` from **mgcv** for smooth fitting.

---

##  Example

```r
library(mgcv)

n <- 1000
X1 <- rnorm(n)
X2 <- sapply(X1, function(x) rnorm(1, x^2, 1))
X <- data.frame(X1, X2)

# Automatic family selection
CPCM_graph_estimate(X, family_of_distributions = "Sequential choice")

# Force Gaussian
CPCM_graph_estimate(X, family_of_distributions = "Gaussian")
```

---

##  References

- Bodik, J. & Chavez-Demoulin, V. (2023): *Identifiability of causal graphs under non-additive conditionally parametric causal models*, arXiv:2303.15376 (https://arxiv.org/abs/2303.15376)

