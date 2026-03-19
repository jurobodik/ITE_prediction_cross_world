## ============================================================
## Copy-paste script: gamma sweep vs gap, rho_eps in {-0.5,0,0.5}
## Baseline DGP, lm only, using your estimate_rho_L().
## Plot uses math-mode axis labels + 95% CI bands around lines.
## ============================================================

## -----------------------
## Baseline DGP 
## -----------------------
simulate_baseline <- function(n, rho_eps = 0.5, gamma = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  X1 <- rnorm(n)
  Z1 <- rnorm(n)
  
  e0 <- rnorm(n)
  e1 <- rho_eps * e0 + sqrt(pmax(1 - rho_eps^2, 0)) * rnorm(n)
  
  m0 <- runif(1, 0.5, 1.5) * X1 + gamma * Z1
  m1 <- runif(1, 0.5, 1.5) * X1 + gamma * Z1
  
  Y0 <- m0 + e0
  Y1 <- m1 + e1
  
  T <- rbinom(n, 1, 0.5)
  Y <- ifelse(T == 1, Y1, Y0)
  
  dat <- data.frame(Y = Y, T = T, X1 = X1, Z1 = Z1, Z2 = Z2)
  list(dat = dat, X1 = X1, Y0 = Y0, Y1 = Y1)
}


## -----------------------
## Kernel helpers (for "truth" rho(x))
## -----------------------
wmean <- function(x, w) sum(w * x) / sum(w)

wcov <- function(x, y, w) {
  mx <- wmean(x, w); my <- wmean(y, w)
  sum(w * (x - mx) * (y - my)) / sum(w)
}

wvar <- function(x, w) wcov(x, x, w)

gauss_kernel <- function(u) exp(-0.5 * u^2)

kernel_weights_X <- function(X_mat, x0, h) {
  X_mat <- as.matrix(X_mat); x0 <- as.numeric(x0)
  d2 <- rowSums((X_mat - matrix(x0, nrow = nrow(X_mat), ncol = length(x0), byrow = TRUE))^2)
  w <- gauss_kernel(sqrt(d2) / h)
  w + 1e-12
}

true_rho_at_x <- function(X, Y0, Y1, x_eval, h_truth = 0.5, eps = 1e-12) {
  X <- as.matrix(X); x_eval <- as.matrix(x_eval)
  out <- numeric(nrow(x_eval))
  for (j in seq_len(nrow(x_eval))) {
    w <- kernel_weights_X(X, x_eval[j, ], h_truth)
    c01 <- wcov(Y0, Y1, w)
    v0  <- wvar(Y0, w)
    v1  <- wvar(Y1, w)
    denom <- sqrt(pmax(v0, eps) * pmax(v1, eps))
    out[j] <- max(-1, min(1, c01 / denom))
  }
  out
}


## -----------------------
## One replication (lm only)
## -----------------------
run_one_rep_lm <- function(
    n = 1000,
    rho_eps = 0.5,
    gamma = 0.5,
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    x_eval = NULL,
    seed = NULL,
    compute_B_L = TRUE
) {
  sim <- simulate_baseline(n = n, rho_eps = rho_eps, gamma = gamma, seed = seed)
  dat <- sim$dat
  
  if (is.null(x_eval)) {
    x_eval <- matrix(
      seq(quantile(dat$X1, 0.05), quantile(dat$X1, 0.95), length.out = 25),
      ncol = 1
    )
    colnames(x_eval) <- "X1"
  } else {
    x_eval <- as.matrix(x_eval)
    if (ncol(x_eval) != 1) stop("This script assumes 1D X (X1).")
    colnames(x_eval) <- "X1"
  }
  
  rho_true_x <- true_rho_at_x(
    X = matrix(sim$X1, ncol = 1),
    Y0 = sim$Y0,
    Y1 = sim$Y1,
    x_eval = x_eval,
    h_truth = h_truth
  )
  
  out_hat <- estimate_rho_L(
    data = dat,
    y_col = "Y",
    t_col = "T",
    x_cols = "X1",
    z_cols = c("Z1"),
    x_eval = x_eval,
    mu_method = "lm",
    compute_B_L = compute_B_L,
    hA = hA,
    hSigma = hSigma
  )
  
  rhoL_x <- out_hat$rho_L_hat
  tilde_x <- if (compute_B_L) out_hat$rho_tilde_L_hat else rep(NA_real_, length(rhoL_x))
  
  list(
    gap_rho_L_mean = mean(rho_true_x - rhoL_x, na.rm = TRUE),
    gap_rho_tilde_mean = if (compute_B_L) mean(rho_true_x - tilde_x, na.rm = TRUE) else NA_real_,
    valid_rho_L = mean(rhoL_x <= rho_true_x + 1e-8, na.rm = TRUE),
    valid_rho_tilde = if (compute_B_L) mean(tilde_x <= rho_true_x + 1e-8, na.rm = TRUE) else NA_real_
  )
}

## -----------------------
## Monte Carlo loop (lm only)
## -----------------------
run_mc_lm <- function(
    R = 100,
    n = 1000,
    rho_eps = 0.5,
    gamma = 0.5,
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    seed = 1,
    compute_B_L = TRUE
) {
  out <- vector("list", R)
  for (r in seq_len(R)) {
    one <- run_one_rep_lm(
      n = n, rho_eps = rho_eps, gamma = gamma,
      hA = hA, hSigma = hSigma, h_truth = h_truth,
      seed = seed + r,
      compute_B_L = compute_B_L
    )
    out[[r]] <- data.frame(
      rep = r,
      gap_rho_L_mean = one$gap_rho_L_mean,
      gap_rho_tilde_mean = one$gap_rho_tilde_mean,
      valid_rho_L = one$valid_rho_L,
      valid_rho_tilde = one$valid_rho_tilde
    )
  }
  do.call(rbind, out)
}

## -----------------------
## Sweep gamma for multiple rho_eps values with 95% CI
## -----------------------
sweep_gamma_by_rhoeps_lm <- function(
    gamma_grid,
    rho_eps_grid = c(-0.5, 0, 0.5),
    R = 100,
    n = 1000,
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    seed = 1,
    compute_B_L = TRUE,
    verbose = TRUE
) {
  res <- list()
  k <- 1
  
  for (re in rho_eps_grid) {
    if (verbose) cat("rho_eps =", re, "\n")
    for (g in seq_along(gamma_grid)) {
      gam <- gamma_grid[g]
      if (verbose) cat("  gamma =", gam, "\n")
      
      mc <- run_mc_lm(
        R = R, n = n,
        rho_eps = re, gamma = gam,
        hA = hA, hSigma = hSigma, h_truth = h_truth,
        seed = seed + 100000 * which(rho_eps_grid == re) + 1000 * g,
        compute_B_L = compute_B_L
      )
      
      ## Means + standard errors across Monte Carlo reps
      qL <- quantile(mc$gap_rho_L_mean, probs = c(0.05, 0.95), na.rm = TRUE)
      qt <- quantile(mc$gap_rho_tilde_mean, probs = c(0.05, 0.95), na.rm = TRUE)
      
      m_gapL <- mean(mc$gap_rho_L_mean, na.rm = TRUE)
      m_gapt <- mean(mc$gap_rho_tilde_mean, na.rm = TRUE)
      
      res[[k]] <- data.frame(
        gamma = gam,
        rho_eps = re,
        gap_rho_L_mean = m_gapL,
        gap_rho_L_lo = as.numeric(qL[1]),
        gap_rho_L_hi = as.numeric(qL[2]),
        gap_rho_tilde_mean = m_gapt,
        gap_rho_tilde_lo = as.numeric(qt[1]),
        gap_rho_tilde_hi = as.numeric(qt[2]),
        valid_rho_L_mean = mean(mc$valid_rho_L, na.rm = TRUE),
        valid_rho_tilde_mean = mean(mc$valid_rho_tilde, na.rm = TRUE)
      )
      k <- k + 1
    }
  }
  
  do.call(rbind, res)
}

## ============================================================
## RUN IT
## ============================================================

gamma_grid <- seq(0, 2, by = 0.25)
rho_eps_grid <- c(-1, -0.5, 0, 0.5, 1)

res <- sweep_gamma_by_rhoeps_lm(
  gamma_grid = gamma_grid,
  rho_eps_grid = rho_eps_grid,
  R = 100,         ## increase if you want smoother CI
  n = 1000,
  hA = 0.5,
  hSigma = 0.5,
  h_truth = 0.5,
  seed = 123,
  compute_B_L = TRUE,
  verbose = TRUE
)

print(res)




















library(ggplot2)

## ============================================================
## Long data for two panels
## ============================================================

df_tilde <- data.frame(
  gamma = res$gamma,
  rho_eps = res$rho_eps,
  mean = res$gap_rho_tilde_mean,
  lo   = res$gap_rho_tilde_lo,
  hi   = res$gap_rho_tilde_hi,
  estimator = "tilde"
)

df_rho <- data.frame(
  gamma = res$gamma,
  rho_eps = res$rho_eps,
  mean = res$gap_rho_L_mean,
  lo   = res$gap_rho_L_lo,
  hi   = res$gap_rho_L_hi,
  estimator = "rho"
)

df <- rbind(df_tilde, df_rho)

## ordered legend
rho_levels <- sort(unique(df$rho_eps))
df$rho_eps <- factor(df$rho_eps, levels = rho_levels)

## ribbon alpha per panel (optional)
df$alpha_ribbon <- ifelse(df$estimator == "tilde", 0.18, 0.10)

## ---- parsed facet titles: include "- rho" because it's a GAP ----
df$estimator_plotmath <- ifelse(df$estimator == "tilde",
                                "rho - hat(tilde(rho))",
                                "rho - hat(rho)")
df$estimator_plotmath <- factor(df$estimator_plotmath,
                                levels = c("rho - hat(tilde(rho))",
                                           "rho - hat(rho)"))

## ============================================================
## Plot
## ============================================================

p <- ggplot(df, aes(x = gamma, y = mean, color = rho_eps, group = rho_eps)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi, fill = rho_eps, alpha = alpha_ribbon),
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_wrap(~ estimator_plotmath, nrow = 1, labeller = label_parsed) +
  scale_alpha_identity() +
  labs(
    x = expression(gamma),
    y = 'Average gap',
    color = expression(rho[epsilon])
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title.align = 0.5,
    legend.background = element_rect(color = "black", fill = "white", linewidth = 0.4),
    strip.background = element_blank(),
    strip.text = element_text(size = 13),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))

p

#save
ggsave(
  filename = "rho_estimate_gap.pdf",
  plot = p,
  device = cairo_pdf,   # nicer text embedding; if not available, remove this line
  width = 7.0, height = 3.6, units = "in"
)










