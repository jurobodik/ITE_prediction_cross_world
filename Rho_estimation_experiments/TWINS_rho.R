# ============================
# TWINS rho estimation
# oracle and hat{rho_L} and hat{rho_tilde_L} for binary outcomes 
# ============================
data <- upload_twins_dataset() #can be found in Helpers.R file
data$treatment <- sample(c(0, 1), size = length(data$Y_0), replace = TRUE)
head(data)

# X = observed covariates used for evaluation; Z = auxiliary covariates used only for rho estimation
x_cols <- c("X2", "X3") #choose any pair
z_cols <- paste0("X", c(1,4,5:51))


# ============================
library(ranger)
mu_method <- "ranger" # or 'lm' for linear model



# Fixed settings for bandwidths and whether to compute B_L and oracle rho
hA <- 0.5                     # bandwidth for weights in A(x)
hSigma <- 0.5                 # bandwidth for weights in sigma_t(x)
compute_B_L <- TRUE           # compute B_L(x) and rho_tilde_L(x)
estimate_oracle_rho <- TRUE   # estimating true rho using both (Y0, Y1)

# ----------------------------
# Build observed dataset
# ----------------------------
tw <- data
X_full <- as.data.frame(tw$X)
if (is.null(colnames(X_full))) colnames(X_full) <- paste0("X", seq_len(ncol(X_full)))

T_vec <- as.integer(tw$treatment)
stopifnot(all(T_vec %in% c(0L, 1L)))

# Observed outcome: pick the potential outcome corresponding to assigned treatment
Y1 <- as.numeric(tw$Y_1)
Y0 <- as.numeric(tw$Y_0)
Y_obs <- ifelse(T_vec == 1L, Y1, Y0)

# Combine into one data frame (Y, T, and all covariates X1..X51)
df <- cbind(
  Y = Y_obs,
  T = T_vec,
  X_full
)

# ----------------------------
# Top 3 most frequent (X1, X2)
# ----------------------------
# Evaluate rho(x) only on the most common discrete (X1, X2) combinations.
pairs <- df[, x_cols, drop = FALSE]

freq <- aggregate(
  rep(1, nrow(pairs)),
  by = list(X1 = pairs[[x_cols[1]]], X2 = pairs[[x_cols[2]]]),
  FUN = sum
)
names(freq)[3] <- "n"

freq <- freq[order(-freq$n, freq$X1, freq$X2), ]
top3 <- head(freq, 3)
x_eval <- as.matrix(top3[, c("X1", "X2"), drop = FALSE])

top3

# ----------------------------
# Prepare design matrices
# ----------------------------
# X used in kernel weights; Z used as "auxiliary" features in mu_t(x,z).
X <- as.matrix(df[, x_cols, drop = FALSE])
Z <- as.matrix(df[, z_cols, drop = FALSE])
Y <- as.numeric(df$Y)
T <- as.integer(df$T)

# Standardized feature names X1..Xp, Z1..Zq for modeling convenience
XZ <- as.data.frame(cbind(X, Z))
names(XZ) <- c(paste0("X", seq_len(ncol(X))), paste0("Z", seq_len(ncol(Z))))
dat <- cbind(Y = Y, T = T, XZ)

# Split by treatment arm, and drop T (arm-specific models)
d0 <- subset(dat, T == 0, select = -T)
d1 <- subset(dat, T == 1, select = -T)

# ----------------------------
# Fit mu models (probabilities)
# ----------------------------
# Goal: mu_t(x,z) = P(Y=1 | T=t, X=x, Z=z), t=0,1.
# - ranger: classification forest with probability=TRUE
# - lm: linear probability model (LPM), predictions clamped to [0,1]
if (mu_method == "ranger") {
  d0$Y <- factor(d0$Y, levels = c(0, 1))
  d1$Y <- factor(d1$Y, levels = c(0, 1))
  
  fit0 <- ranger::ranger(Y ~ ., data = d0, num.trees = 500,
                         probability = TRUE, respect.unordered.factors = TRUE)
  fit1 <- ranger::ranger(Y ~ ., data = d1, num.trees = 500,
                         probability = TRUE, respect.unordered.factors = TRUE)
  
} else if (mu_method == "lm") {
  fml <- as.formula(paste("Y ~", paste(colnames(XZ), collapse = " + ")))
  fit0 <- lm(fml, data = d0)
  fit1 <- lm(fml, data = d1)
  
} else {
  stop("mu_method must be 'ranger' or 'lm'.")
}

# ----------------------------
# Helpers (kept minimal)
# ----------------------------
eps <- 1e-12
gauss_kernel <- function(u) exp(-0.5 * u^2)  # unnormalized Gaussian kernel

# ----------------------------
# Loop over x_eval and compute rho
# ----------------------------
# A(x) = Cov(mu0(x,Z), mu1(x,Z) | X=x) via kernel weights in X
# sigma_t(x) = Var(Y | T=t, X=x) via local weighted variance within each arm
# rho_L(x) = A(x) / sqrt(sigma0(x)*sigma1(x))
# rho_tilde_L(x) = (A(x) - B_L(x)) / sqrt(sigma0(x)*sigma1(x))
nxe <- nrow(x_eval)

rho_L_hat <- numeric(nxe)
rho_tilde_L_hat <- if (compute_B_L) numeric(nxe) else NULL
B_L_hat <- if (compute_B_L) numeric(nxe) else NULL

A_hat <- numeric(nxe)
sigma0_hat <- numeric(nxe)
sigma1_hat <- numeric(nxe)

oracle_rho_hat <- if (estimate_oracle_rho) numeric(nxe) else NULL

for (j in seq_len(nxe)) {
  x0 <- as.numeric(x_eval[j, ])
  
  # Kernel weights w_i(x0) based on distance between X_i and x0
  d2_all <- rowSums((X - matrix(x0, nrow = nrow(X), ncol = length(x0), byrow = TRUE))^2)
  hA <- 0.5                     # bandwidth for weights in A(x)
  wA <- gauss_kernel(sqrt(d2_all) / hA) + 1e-12
  wsum <- sum(wA)
  
  # Build prediction dataset: X fixed to x0, Z varies across observed Z_i
  Z_df <- as.data.frame(Z); names(Z_df) <- paste0("Z", seq_len(ncol(Z)))
  X_df <- as.data.frame(matrix(x0, nrow = nrow(Z_df), ncol = length(x0), byrow = TRUE))
  names(X_df) <- paste0("X", seq_len(ncol(X)))
  newdat <- cbind(X_df, Z_df)
  
  # Predict mu0(x0,Z_i) and mu1(x0,Z_i)
  if (mu_method == "ranger") {
    pp0 <- predict(fit0, data = newdat)$predictions
    pp1 <- predict(fit1, data = newdat)$predictions
    c0 <- if ("1" %in% colnames(pp0)) "1" else colnames(pp0)[2]
    c1 <- if ("1" %in% colnames(pp1)) "1" else colnames(pp1)[2]
    mu0 <- as.numeric(pp0[, c0])
    mu1 <- as.numeric(pp1[, c1])
  } else {
    mu0 <- as.numeric(predict(fit0, newdata = newdat))
    mu1 <- as.numeric(predict(fit1, newdata = newdat))
    mu0 <- pmin(pmax(mu0, 0), 1)
    mu1 <- pmin(pmax(mu1, 0), 1)
  }
  
  # Weighted covariance between mu0 and mu1 (A_hat(x0))
  m0 <- sum(wA * mu0) / wsum
  m1 <- sum(wA * mu1) / wsum
  Aj <- sum(wA * (mu0 - m0) * (mu1 - m1)) / wsum
  A_hat[j] <- Aj
  
  # Estimate sigma_t(x0): local weighted variance of observed Y within each arm
  for (tval in 0:1) {
    idx <- which(T == tval)
    Xt <- X[idx, , drop = FALSE]
    Yt <- Y[idx]
    d2_t <- rowSums((Xt - matrix(x0, nrow = nrow(Xt), ncol = length(x0), byrow = TRUE))^2)
    wt <- gauss_kernel(sqrt(d2_t) / hSigma) + 1e-12
    wtsum <- sum(wt)
    my <- sum(wt * Yt) / wtsum
    vy <- sum(wt * (Yt - my)^2) / wtsum
    if (tval == 0) sigma0_hat[j] <- max(vy, eps)
    if (tval == 1) sigma1_hat[j] <- max(vy, eps)
  }
  
  denom <- max(sqrt(sigma0_hat[j] * sigma1_hat[j]), eps)
  
  # rho_L(x0)
  rho_L_hat[j] <- max(-1, min(1, Aj / denom))
  
  # B_L(x0) for binary: sigma_t(x0,z) = mu_t(x0,z)(1-mu_t(x0,z))
  if (compute_B_L) {
    s0_xz <- pmax(mu0 * (1 - mu0), eps)
    s1_xz <- pmax(mu1 * (1 - mu1), eps)
    BLj <- sum(wA * sqrt(s0_xz * s1_xz)) / wsum
    B_L_hat[j] <- BLj
    rho_tilde_L_hat[j] <- max(-1, min(1, (Aj - BLj) / denom))
  }
  
  # "True" rho(x0) based on potential outcomes (TWINS-only sanity check)
  if (estimate_oracle_rho) {
    mY0 <- sum(wA * Y0) / wsum
    mY1 <- sum(wA * Y1) / wsum
    c01 <- sum(wA * (Y0 - mY0) * (Y1 - mY1)) / wsum
    v0o <- sum(wA * (Y0 - mY0)^2) / wsum
    v1o <- sum(wA * (Y1 - mY1)^2) / wsum
    oracle_rho_hat[j] <- max(-1, min(1, c01 / max(sqrt(v0o * v1o), eps)))
  }
}

# ----------------------------
# Results table (top 3 only)
# ----------------------------
rho_tab <- top3[, c("X1", "X2", "n"), drop = FALSE]
rho_tab$rho_L <- round(rho_L_hat, 3)

if (compute_B_L) {rho_tab$rho_tilde_L <- round(rho_tilde_L_hat, 3)}

if (estimate_oracle_rho) {
  rho_tab$rho_oracle <- round(oracle_rho_hat, 3)
}

rho_tab

# Optional prettier printing if you have knitr installed
if (requireNamespace("knitr", quietly = TRUE)) {
  knitr::kable(rho_tab, caption = paste0("Rho estimates (mu_method = ", mu_method, ") on the 3 most frequent (X1, X2) pairs"))
}
