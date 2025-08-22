# Cross-World Assumption and Refining Prediction Intervals for Individual Treatment Effects

[![R](https://img.shields.io/badge/R-4.2+-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://arxiv.org/abs/2507.12581)

We develop methods for constructing **valid and informative prediction intervals for individual treatment effects (ITE)** for new units not within the study, under assumptions on the **cross-world correlation parameter ρ(x) = corr(Y(1), Y(0) | X = x)**.

Assuming values or bounds on ρ(x) allows us to refine the uncertainty quantification around ITE estimates—balancing informativeness and validity.

![image](https://github.com/user-attachments/assets/8597020a-61f6-4fbc-97f5-46dfe4fdadd1)

![image](https://github.com/user-attachments/assets/ed7e234e-bc1f-46cd-a66c-2b742987f8e2)

---
---

## ⚙️ Core Function: `CW_rho_intervals`

The function `CW_rho_intervals()` constructs ITE prediction intervals under a user-specified cross-world correlation parameter ρ.

Supported options:
- CATE estimators: `T-learner`, `Causal Forest`
- Conformal methods: `CQR`, `Weighted CQR`, `PCS`, `CLEAR`, `Naive`
- Optional bootstrap-based confidence intervals

---

## 🧪 Example Usage (R or Python)

```r
source("CW_rho_intervals_function.R")

set.seed(0)
n <- 1000; d <- 5; rho <- 0.5
X <- as.data.frame(matrix(rnorm(n * d), ncol = d))
treatment <- rbinom(n, 1, 0.5)
epsilons <- MASS::mvrnorm(n, mu = c(0,0), Sigma = matrix(c(1, rho, rho, 1), 2))
Y0 <- X[,1] + epsilons[,1]
Y1 <- X[,1] + 1 + X[,2]^2 + epsilons[,2]
Y <- ifelse(treatment == 1, Y1, Y0)

new_points <- as.data.frame(matrix(rnorm(10 * d), ncol = d)) #at what points do you want to estimate ITE?

result <- CW_rho_intervals(X, Y, treatment, new_points, rho,
                          desired_coverage = 0.9, 
                          conformal = "CQR",
                          weighted_conformal = FALSE,
                          add_confidence_intervals = TRUE)

print(result$CATE)
print(result$lower)
print(result$upper)
```
We also implemented a translation of CW_rho_intervals into python:

```python
from CW_rho_intervals_function import CW_rho_intervals
import numpy as np

np.random.seed(0)
n, d = 1000, 5
rho = 0.5

X = np.random.randn(n, d)
treatment = np.random.binomial(1, 0.5, size=n)

cov = [[1, rho], [rho, 1]]
eps = np.random.multivariate_normal([0, 0], cov, size=n)
Y0 = X[:, 0] + eps[:, 0]
Y1 = X[:, 0] + 1 + X[:, 1]**2 + eps[:, 1]
Y = np.where(treatment == 1, Y1, Y0)

new_points = np.random.randn(10, d) #at what points do you want to estimate ITE?

result = CW_rho_intervals(X, Y, treatment, new_points, rho=rho,
                         desired_coverage = 0.9, 
                         conformal="CQR", weighted_conformal=False,
                         add_confidence_intervals=True)

print("CATE:", result["CATE"])
print("Lower:", result["lower"])
print("Upper:", result["upper"])
```
## How to choose ρ? 
Ask: “What proportion of the total noise variance can be attributed to hidden components affecting both potential outcomes similarly? ” In the additive model 

Y(1) = μ₁(X) + H + ε̃₁  
Y(0) = μ₀(X) + H + ε̃₀  

where:
- `H` is a hidden variable shared across potential outcomes,
- `ε̃₁` and `ε̃₀` are independent noise terms.

In this setup, ρ corresponds to the fraction of noise variance explained by the shared hidden component:

ρ = var(H) / [var(H) + var(ε̃)]

This gives a principled way to discuss plausible values of ρ based on how much of the total outcome noise is due to shared latent traits.

If no interpretation or domain knowledge is available, we recommend using ρ = 0 when coverage guarantees are of strong importance; otherwise, we recommend using ρ = 0.5 for the best performance across most practical scenarios. If treatment has little effect, do not hesitate to use ρ>=0.75. 

When in doubt, run sensitivity analysis over a range of ρ values.










