# Cross-World Assumptions and Their Role in Refining Prediction Intervals for Individual Treatment Effects

[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![R](https://img.shields.io/badge/R-4.2+-blue.svg)](https://www.r-project.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1234567.svg)](https://doi.org/10.5281/zenodo.1234567)


We develop methods for constructing **valid and informative prediction intervals for individual treatment effects (ITE)** under assumptions on the **cross-world correlation parameter**  
**ρ(x) = corr(Y(1), Y(0) | X = x)**.

Assuming values or bounds on ρ(x) allows us to refine the uncertainty quantification around ITE estimates—balancing informativeness and validity.
![image](https://github.com/user-attachments/assets/8597020a-61f6-4fbc-97f5-46dfe4fdadd1)

---
---

## ⚙️ Core Function: `D_rho_intervals`

The function `D_rho_intervals()` constructs ITE prediction intervals under a user-specified cross-world correlation parameter ρ.

Supported options:
- CATE estimators: `T-learner`, `Causal Forest`
- Conformal methods: `CQR`, `Weighted CQR`, `PCS`, `CLEAR`, `Naive`
- Optional bootstrap-based confidence intervals

---

## 🧪 Example Usage

```r
source("Main_function.R")

set.seed(0)
n <- 1000; d <- 5; rho <- 0.5
X <- as.data.frame(matrix(rnorm(n * d), ncol = d))
treatment <- rbinom(n, 1, 0.5)
epsilons <- MASS::mvrnorm(n, mu = c(0,0), Sigma = matrix(c(1, rho, rho, 1), 2))
Y0 <- X[,1] + epsilons[,1]
Y1 <- X[,1] + 1 + X[,2]^2 + epsilons[,2]
Y <- ifelse(treatment == 1, Y1, Y0)

new_points <- as.data.frame(matrix(rnorm(10 * d), ncol = d))

result <- D_rho_intervals(X, Y, treatment, new_points, rho,
                          conformal = "CQR", weighted_conformal = FALSE,
                          add_confidence_intervals = TRUE)

print(result$CATE)
print(result$lower)
print(result$upper)

```python
from Main_function import D_rho_intervals
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

new_points = np.random.randn(10, d)

result = D_rho_intervals(X, Y, treatment, new_points, rho=rho,
                         conformal="CQR", weighted_conformal=False,
                         add_confidence_intervals=True)

print("CATE:", result["CATE"])
print("Lower:", result["lower"])
print("Upper:", result["upper"])












