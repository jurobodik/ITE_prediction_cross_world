## ============================================================
## FULL SCRIPT: gamma sweep vs gap, rho_eps grid
## Baseline DGP, lm only, using your estimate_rho_L()
## Higher-dimensional Z (aux covariates) + Monte Carlo CI bands
## Plot: 2 x 3 grid (rows: rho/tilde; cols: dim(Z)=pZ)
## ============================================================

## -----------------------
## Baseline DGP with high-dim Z
## -----------------------
simulate_baseline_hdZ <- function(
    n,
    rho_eps = 0.5,
    gamma = 0.5,
    pZ = 10,
    s = NULL,          # number of Z coords that truly matter; default = pZ
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  X1 <- rnorm(n)
  
  Z <- matrix(rnorm(n * pZ), nrow = n, ncol = pZ)
  colnames(Z) <- paste0("Z", seq_len(pZ))
  
  t <- if (is.null(s)) pZ else s
  if (t > pZ) stop("s must be <= pZ")
  
  a <- numeric(pZ)
  a[seq_len(t)] <- rnorm(t)
  a <- a / sqrt(sum(a^2) + 1e-12)
  H <- as.numeric(Z %*% a)
  
  e0 <- rnorm(n)
  e1 <- rho_eps * e0 + sqrt(pmax(1 - rho_eps^2, 0)) * rnorm(n)
  
  b0 <- runif(1, 0.5, 1.5)
  b1 <- runif(1, 0.5, 1.5)
  
  m0 <- b0 * X1 + gamma * H
  m1 <- b1 * X1 + gamma * H
  
  Y0 <- m0 + e0
  Y1 <- m1 + e1
  
  Tobs <- rbinom(n, 1, 0.5)
  Yobs <- ifelse(Tobs == 1, Y1, Y0)
  
  dat <- data.frame(Y = Yobs, T = Tobs, X1 = X1)
  dat <- cbind(dat, as.data.frame(Z))
  
  list(dat = dat, X1 = X1, Z = Z, Y0 = Y0, Y1 = Y1)
}

## -----------------------
## Kernel helpers (truth rho(x) = cor(Y0,Y1 | X=x))
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
run_one_rep_lm_hdZ <- function(
    n = 1000,
    rho_eps = 0.5,
    gamma = 0.5,
    pZ = 10,
    s = NULL,
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    x_eval = NULL,
    seed = NULL,
    compute_B_L = TRUE
) {
  sim <- simulate_baseline_hdZ(
    n = n, rho_eps = rho_eps, gamma = gamma, pZ = pZ, s = s, seed = seed
  )
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
  
  z_cols <- paste0("Z", seq_len(pZ))
  
  out_hat <- estimate_rho_L(
    data = dat,
    y_col = "Y",
    t_col = "T",
    x_cols = "X1",
    z_cols = z_cols,
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
    gap_rho_tilde_mean = if (compute_B_L) mean(rho_true_x - tilde_x, na.rm = TRUE) else NA_real_
  )
}

## -----------------------
## Monte Carlo loop: MEAN + quantile CI across reps
## -----------------------
run_mc_lm_hdZ <- function(
    R = 100,
    n = 1000,
    rho_eps = 0.5,
    gamma = 0.5,
    pZ = 10,
    s = NULL,
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    seed = 1,
    compute_B_L = TRUE,
    ci_level = 0.95
) {
  gapL <- numeric(R)
  gapt <- numeric(R)
  
  for (r in seq_len(R)) {
    one <- run_one_rep_lm_hdZ(
      n = n, rho_eps = rho_eps, gamma = gamma,
      pZ = pZ, s = s,
      hA = hA, hSigma = hSigma, h_truth = h_truth,
      seed = seed + r,
      compute_B_L = compute_B_L
    )
    gapL[r] <- one$gap_rho_L_mean
    gapt[r] <- one$gap_rho_tilde_mean
  }
  
  alpha <- 1 - ci_level
  q <- c(alpha / 2, 1 - alpha / 2)
  
  list(
    gap_rho_L_mean = mean(gapL, na.rm = TRUE),
    gap_rho_L_lo   = as.numeric(quantile(gapL, probs = q[1], na.rm = TRUE)),
    gap_rho_L_hi   = as.numeric(quantile(gapL, probs = q[2], na.rm = TRUE)),
    
    gap_rho_tilde_mean = mean(gapt, na.rm = TRUE),
    gap_rho_tilde_lo   = as.numeric(quantile(gapt, probs = q[1], na.rm = TRUE)),
    gap_rho_tilde_hi   = as.numeric(quantile(gapt, probs = q[2], na.rm = TRUE))
  )
}

## -----------------------
## Sweep gamma for multiple rho_eps values (WITH CI)
## -----------------------
sweep_gamma_by_rhoeps_lm_hdZ <- function(
    gamma_grid,
    rho_eps_grid = c(-0.5, 0, 0.5),
    R = 100,
    n = 1000,
    pZ = 10,
    s = NULL,
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    seed = 1,
    compute_B_L = TRUE,
    ci_level = 0.95,
    verbose = TRUE
) {
  out <- vector("list", length(rho_eps_grid) * length(gamma_grid))
  k <- 1
  
  for (re in rho_eps_grid) {
    if (verbose) cat("rho_eps =", re, "\n")
    for (g in seq_along(gamma_grid)) {
      gam <- gamma_grid[g]
      if (verbose) cat("  gamma =", gam, "\n")
      
      mc <- run_mc_lm_hdZ(
        R = R, n = n,
        rho_eps = re, gamma = gam,
        pZ = pZ, s = s,
        hA = hA, hSigma = hSigma, h_truth = h_truth,
        seed = seed + 100000 * which(rho_eps_grid == re) + 1000 * g,
        compute_B_L = compute_B_L,
        ci_level = ci_level
      )
      
      out[[k]] <- data.frame(
        gamma = gam,
        rho_eps = re,
        pZ = pZ,
        s = if (is.null(s)) pZ else s,
        
        gap_rho_L_mean = mc$gap_rho_L_mean,
        gap_rho_L_lo   = mc$gap_rho_L_lo,
        gap_rho_L_hi   = mc$gap_rho_L_hi,
        
        gap_rho_tilde_mean = mc$gap_rho_tilde_mean,
        gap_rho_tilde_lo   = mc$gap_rho_tilde_lo,
        gap_rho_tilde_hi   = mc$gap_rho_tilde_hi
      )
      k <- k + 1
    }
  }
  
  do.call(rbind, out)
}

## ============================================================
## RUN IT
## ============================================================

gamma_grid <- seq(0, 2, by = 0.25)
rho_eps_grid <- c(-1, -0.5, 0, 0.5, 1)

pZ_list <- c(5, 20, 100)

all_res <- list()
for (pZ in pZ_list) {
  cat("\n=== Running pZ =", pZ, "===\n")
  all_res[[as.character(pZ)]] <- sweep_gamma_by_rhoeps_lm_hdZ(
    gamma_grid = gamma_grid,
    rho_eps_grid = rho_eps_grid,
    R = 100,
    n = 1000,
    pZ = pZ,
    s = min(5, pZ),
    hA = 0.5,
    hSigma = 0.5,
    h_truth = 0.5,
    seed = 123,
    compute_B_L = TRUE,
    ci_level = 0.95,
    verbose = TRUE
  )
}

res <- do.call(rbind, all_res)
print(res)
#save res to file if desired, e.g.:
#saveRDS(res, file = "gamma_sweep_rhoeps_grid_hdZ_results.rds")

## ============================================================
## PLOT: same style as yours (with ribbons)
## ============================================================

library(ggplot2)

res$pZ <- factor(res$pZ, levels = sort(unique(res$pZ)))
pZ_labs <- setNames(paste0("dim(Z) = ", levels(res$pZ)), levels(res$pZ))

df_tilde <- data.frame(
  gamma = res$gamma,
  rho_eps = res$rho_eps,
  pZ = res$pZ,
  mean = res$gap_rho_tilde_mean,
  lo   = res$gap_rho_tilde_lo,
  hi   = res$gap_rho_tilde_hi,
  estimator = "tilde"
)

df_rho <- data.frame(
  gamma = res$gamma,
  rho_eps = res$rho_eps,
  pZ = res$pZ,
  mean = res$gap_rho_L_mean,
  lo   = res$gap_rho_L_lo,
  hi   = res$gap_rho_L_hi,
  estimator = "rho"
)

df <- rbind(df_tilde, df_rho)

rho_levels <- sort(unique(df$rho_eps))
df$rho_eps <- factor(df$rho_eps, levels = rho_levels)

df$alpha_ribbon <- 0.20

df$estimator_plotmath <- ifelse(df$estimator == "tilde",
                                "rho - hat(tilde(rho))",
                                "rho - hat(rho)")
df$estimator_plotmath <- factor(df$estimator_plotmath,
                                levels = c("rho - hat(tilde(rho))",
                                           "rho - hat(rho)"))

y_breaks_fun <- function(lims) sort(unique(c(pretty(lims), 2)))

p <- ggplot(df, aes(x = gamma, y = mean, color = rho_eps, group = rho_eps)) +
  geom_hline(yintercept = 0, color = "red", linewidth = 0.5) +
  geom_ribbon(
    aes(ymin = lo, ymax = hi, fill = rho_eps, alpha = alpha_ribbon),
    color = NA,
    show.legend = FALSE
  ) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2) +
  facet_grid(
    estimator_plotmath ~ pZ,
    labeller = labeller(
      estimator_plotmath = label_parsed,
      pZ = pZ_labs
    )
  ) +
  scale_y_continuous(breaks = y_breaks_fun) +
  scale_alpha_identity() +
  labs(
    x = expression(gamma),
    y = "Average gap",
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

## Optional save
 ggsave(
   filename = "rho_estimate_high_dim.pdf",
   plot = p,
   device = cairo_pdf,
   width = 7.0, height = 4, units = "in"
 )
