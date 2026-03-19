# =============================================================================
# estimate_rho_L()
# =============================================================================
# Estimate the lower bound rho_L(x) for the cross-world correlation rho(x)
# using auxiliary covariates Z
#
# The function can also compute the more conservative bound rho_tilde_L(x),
# and optional bootstrap-based confidence intervals.
#
# -----------------------------------------------------------------------------
# Example
# -----------------------------------------------------------------------------
# # Simulate simple data
# # n <- 1000
# # rho_epsilon <- 0
# # gamma_Z <- 1.5
# #
# # e0 <- rnorm(n)
# # e1 <- rho_epsilon * e0 + sqrt(1 - rho_epsilon^2) * rnorm(n) # cor(e0, e1) = rho_epsilon
# #
# # X1 <- rnorm(n)
# # Z1 <- rnorm(n)
# #
# # tilde_e0 <- gamma_Z * Z1 + e0
# # tilde_e1 <- gamma_Z * Z1 + e1
# #
# # Y0 <- 0.6 * X1 + tilde_e0
# # Y1 <- 0.8 * X1 + tilde_e1
# # rho_true <- cor(tilde_e0, tilde_e1)
# #
# # T <- rbinom(n, 1, 0.5)
# # Y <- ifelse(T == 1, Y1, Y0)
# # dat <- data.frame(Y = Y, T = T, X1 = X1, Z1 = Z1)
# #
# # out <- estimate_rho_L(
# #   data = dat,
# #   y_col = "Y",
# #   t_col = "T",
# #   x_cols = "X1",
# #   z_cols = "Z1",
# #   mu_method = "lm",
# #   ci = FALSE
# # )
# #
# # out$rho_L_mean
# # out$rho_tilde_L_mean
# # rho_true
# =============================================================================

# library(ranger) #only if ranger is used
estimate_rho_L <- function(
    data,
    y_col = "Y",
    t_col = "T",
    x_cols = c("X1"),
    z_cols = c("Z1"),
    x_eval = NULL,
    mu_method = c("ranger", "lm"),
    compute_B_L = TRUE,
    hA = 0.5,
    hSigma = 0.5,
    ranger_num_trees = 500,
    ci = FALSE,
    ci_method = c("bootstrap_full", "bootstrap_fast"),
    n_boot = 50,
    conf_level = 0.95,
    seed = 1,
    store_boot = FALSE,
    ranger_seed = TRUE
) {
  mu_method <- match.arg(mu_method)
  if (ci) ci_method <- match.arg(ci_method)

  # ---------------------------------------------------------------------------
  # 1. Input checks
  # ---------------------------------------------------------------------------
  required_cols <- unique(c(y_col, t_col, x_cols, z_cols))
  if (!all(required_cols %in% names(data))) {
    stop("Some of y_col, t_col, x_cols, or z_cols are not present in `data`.")
  }

  Y_raw <- data[[y_col]]
  T <- data[[t_col]]

  if (is.factor(T)) T <- as.integer(as.character(T))
  T <- as.integer(T)
  if (!all(T %in% c(0L, 1L))) stop("`T` must be coded as 0/1.")

  # Detect whether outcome is binary
  is_binary <- FALSE
  if (is.logical(Y_raw)) {
    is_binary <- TRUE
  } else if (is.factor(Y_raw) && nlevels(Y_raw) == 2) {
    is_binary <- TRUE
  } else {
    uy <- sort(unique(na.omit(as.numeric(Y_raw))))
    is_binary <- (length(uy) <= 2 && all(uy %in% c(0, 1)))
  }

  eps <- 1e-12

  # Convert Y to numeric
  if (is_binary) {
    if (is.logical(Y_raw)) {
      Y <- as.numeric(Y_raw)
    } else if (is.factor(Y_raw) && nlevels(Y_raw) == 2) {
      Y <- as.numeric(Y_raw == levels(Y_raw)[2])
    } else {
      Y <- as.numeric(Y_raw)
    }

    uy <- sort(unique(na.omit(Y)))
    if (!(length(uy) <= 2 && all(uy %in% c(0, 1)))) {
      stop("Binary outcome detected, but `Y` is not coded as 0/1.")
    }
  } else {
    Y <- as.numeric(Y_raw)
  }

  # ---------------------------------------------------------------------------
  # 2. Prepare data objects
  # ---------------------------------------------------------------------------
  X_mat <- as.matrix(data[, x_cols, drop = FALSE])
  storage.mode(X_mat) <- "double"

  X_df_model <- as.data.frame(data[, x_cols, drop = FALSE])
  Z_df_model <- as.data.frame(data[, z_cols, drop = FALSE])

  names(X_df_model) <- paste0("X", seq_len(ncol(X_df_model)))
  names(Z_df_model) <- paste0("Z", seq_len(ncol(Z_df_model)))

  XZ <- cbind(X_df_model, Z_df_model)
  dat <- cbind(Y = Y, T = T, XZ)

  # ---------------------------------------------------------------------------
  # 3. Default evaluation points x_eval
  # ---------------------------------------------------------------------------
  if (is.null(x_eval)) {
    X_df_raw <- data[, x_cols, drop = FALSE]
    n <- nrow(X_df_raw)
    p <- ncol(X_df_raw)

    uniq_counts <- vapply(
      X_df_raw,
      function(col) length(unique(col[!is.na(col)])),
      integer(1)
    )
    is_non_numeric <- vapply(
      X_df_raw,
      function(col) !(is.numeric(col) || is.integer(col)),
      logical(1)
    )

    is_discrete <- all(is_non_numeric | (uniq_counts <= 10) | (uniq_counts / n <= 0.05))

    if (is_discrete) {
      freq <- aggregate(rep(1, n), by = as.list(X_df_raw), FUN = sum)
      names(freq)[ncol(freq)] <- "n"
      freq <- freq[order(-freq$n), , drop = FALSE]
      k <- min(3, nrow(freq))
      x_eval <- freq[seq_len(k), x_cols, drop = FALSE]
    } else {
      if (p > 5) {
        x_eval <- as.data.frame(lapply(X_df_raw, function(col) {
          if (is.numeric(col) || is.integer(col)) {
            median(col, na.rm = TRUE)
          } else {
            u <- unique(col[!is.na(col)])
            u[which.max(tabulate(match(col, u)))]
          }
        }))
        names(x_eval) <- x_cols
      } else {
        qs <- lapply(
          X_df_raw,
          function(col) quantile(col, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, type = 7)
        )
        x_eval <- expand.grid(
          lapply(qs, function(v) as.numeric(unique(v))),
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        )
        names(x_eval) <- x_cols
      }
    }
  } else {
    x_eval <- as.data.frame(x_eval)
    if (!all(x_cols %in% names(x_eval))) {
      stop("`x_eval` must contain all columns listed in `x_cols`.")
    }
    x_eval <- x_eval[, x_cols, drop = FALSE]
  }

  x_eval_mat <- as.matrix(x_eval)
  storage.mode(x_eval_mat) <- "double"

  # ---------------------------------------------------------------------------
  # 4. Helper functions
  # ---------------------------------------------------------------------------
  wmean <- function(x, w) sum(w * x) / sum(w)

  wcov <- function(x, y, w) {
    mx <- wmean(x, w)
    my <- wmean(y, w)
    sum(w * (x - mx) * (y - my)) / sum(w)
  }

  wvar <- function(x, w) wcov(x, x, w)

  gauss_kernel <- function(u) exp(-0.5 * u^2)

  kernel_weights_X <- function(X, x0, h) {
    d2 <- rowSums((X - matrix(x0, nrow = nrow(X), ncol = length(x0), byrow = TRUE))^2)
    gauss_kernel(sqrt(d2) / h) + eps
  }

  # ---------------------------------------------------------------------------
  # 5. Fit nuisance models mu_t(x, z)
  # ---------------------------------------------------------------------------
  d0 <- subset(dat, T == 0, select = -T)
  d1 <- subset(dat, T == 1, select = -T)

  fit0 <- NULL
  fit1 <- NULL
  pred_mu0_train <- NULL
  pred_mu1_train <- NULL

  base_rf_seed <- if (ranger_seed && !is.null(seed)) as.integer(seed) else NULL

  if (mu_method == "ranger") {
    if (!requireNamespace("ranger", quietly = TRUE)) {
      stop("Package `ranger` is required when mu_method = 'ranger'.")
    }

    if (is_binary) {
      d0$Y <- factor(d0$Y, levels = c(0, 1))
      d1$Y <- factor(d1$Y, levels = c(0, 1))

      args0 <- list(
        formula = Y ~ ., data = d0, num.trees = ranger_num_trees,
        probability = TRUE, respect.unordered.factors = TRUE
      )
      args1 <- list(
        formula = Y ~ ., data = d1, num.trees = ranger_num_trees,
        probability = TRUE, respect.unordered.factors = TRUE
      )

      if (!is.null(base_rf_seed)) {
        args0$seed <- base_rf_seed + 101L
        args1$seed <- base_rf_seed + 102L
      }

      fit0 <- do.call(ranger::ranger, args0)
      fit1 <- do.call(ranger::ranger, args1)

      p0_tr <- predict(fit0, data = d0)$predictions
      p1_tr <- predict(fit1, data = d1)$predictions
      col0 <- if ("1" %in% colnames(p0_tr)) "1" else colnames(p0_tr)[2]
      col1 <- if ("1" %in% colnames(p1_tr)) "1" else colnames(p1_tr)[2]

      pred_mu0_train <- as.numeric(p0_tr[, col0])
      pred_mu1_train <- as.numeric(p1_tr[, col1])

    } else {
      args0 <- list(
        formula = Y ~ ., data = d0, num.trees = ranger_num_trees,
        respect.unordered.factors = TRUE
      )
      args1 <- list(
        formula = Y ~ ., data = d1, num.trees = ranger_num_trees,
        respect.unordered.factors = TRUE
      )

      if (!is.null(base_rf_seed)) {
        args0$seed <- base_rf_seed + 111L
        args1$seed <- base_rf_seed + 112L
      }

      fit0 <- do.call(ranger::ranger, args0)
      fit1 <- do.call(ranger::ranger, args1)

      pred_mu0_train <- as.numeric(predict(fit0, data = d0)$predictions)
      pred_mu1_train <- as.numeric(predict(fit1, data = d1)$predictions)
    }

  } else if (mu_method == "lm") {
    fml <- as.formula(paste("Y ~", paste(names(XZ), collapse = " + ")))

    if (is_binary) {
      fit0 <- glm(fml, data = d0, family = binomial())
      fit1 <- glm(fml, data = d1, family = binomial())
      pred_mu0_train <- as.numeric(predict(fit0, newdata = d0, type = "response"))
      pred_mu1_train <- as.numeric(predict(fit1, newdata = d1, type = "response"))
    } else {
      fit0 <- lm(fml, data = d0)
      fit1 <- lm(fml, data = d1)
      pred_mu0_train <- as.numeric(predict(fit0, newdata = d0))
      pred_mu1_train <- as.numeric(predict(fit1, newdata = d1))
    }
  }

  predict_mu_at_x <- function(x0_vec) {
    newX <- as.data.frame(
      matrix(x0_vec, nrow = nrow(Z_df_model), ncol = ncol(X_df_model), byrow = TRUE)
    )
    names(newX) <- names(X_df_model)
    newdat <- cbind(newX, Z_df_model)

    if (mu_method == "ranger") {
      if (is_binary) {
        pp0 <- predict(fit0, data = newdat)$predictions
        pp1 <- predict(fit1, data = newdat)$predictions
        c0 <- if ("1" %in% colnames(pp0)) "1" else colnames(pp0)[2]
        c1 <- if ("1" %in% colnames(pp1)) "1" else colnames(pp1)[2]
        mu0 <- as.numeric(pp0[, c0])
        mu1 <- as.numeric(pp1[, c1])
      } else {
        mu0 <- as.numeric(predict(fit0, data = newdat)$predictions)
        mu1 <- as.numeric(predict(fit1, data = newdat)$predictions)
      }
    } else {
      if (is_binary) {
        mu0 <- as.numeric(predict(fit0, newdata = newdat, type = "response"))
        mu1 <- as.numeric(predict(fit1, newdata = newdat, type = "response"))
      } else {
        mu0 <- as.numeric(predict(fit0, newdata = newdat))
        mu1 <- as.numeric(predict(fit1, newdata = newdat))
      }
    }

    if (is_binary) {
      mu0 <- pmin(pmax(mu0, eps), 1 - eps)
      mu1 <- pmin(pmax(mu1, eps), 1 - eps)
    }

    list(mu0 = mu0, mu1 = mu1)
  }

  # ---------------------------------------------------------------------------
  # 6. Optional models for sigma_t(x, z) used in B_L
  # ---------------------------------------------------------------------------
  predict_sigma_at_x <- NULL

  if (compute_B_L && !is_binary) {
    r0_2 <- (d0$Y - pred_mu0_train)^2
    r1_2 <- (d1$Y - pred_mu1_train)^2

    v0_dat <- cbind(r2 = r0_2, d0[, setdiff(names(d0), "Y"), drop = FALSE])
    v1_dat <- cbind(r2 = r1_2, d1[, setdiff(names(d1), "Y"), drop = FALSE])

    if (mu_method == "ranger") {
      argsv0 <- list(
        formula = r2 ~ ., data = v0_dat, num.trees = ranger_num_trees,
        respect.unordered.factors = TRUE
      )
      argsv1 <- list(
        formula = r2 ~ ., data = v1_dat, num.trees = ranger_num_trees,
        respect.unordered.factors = TRUE
      )

      if (!is.null(base_rf_seed)) {
        argsv0$seed <- base_rf_seed + 201L
        argsv1$seed <- base_rf_seed + 202L
      }

      vfit0 <- do.call(ranger::ranger, argsv0)
      vfit1 <- do.call(ranger::ranger, argsv1)

      predict_sigma_at_x <- function(x0_vec) {
        newX <- as.data.frame(
          matrix(x0_vec, nrow = nrow(Z_df_model), ncol = ncol(X_df_model), byrow = TRUE)
        )
        names(newX) <- names(X_df_model)
        newdat <- cbind(newX, Z_df_model)

        s0 <- pmax(as.numeric(predict(vfit0, data = newdat)$predictions), eps)
        s1 <- pmax(as.numeric(predict(vfit1, data = newdat)$predictions), eps)
        list(s0 = s0, s1 = s1)
      }

    } else {
      vfit0 <- lm(r2 ~ ., data = v0_dat)
      vfit1 <- lm(r2 ~ ., data = v1_dat)

      predict_sigma_at_x <- function(x0_vec) {
        newX <- as.data.frame(
          matrix(x0_vec, nrow = nrow(Z_df_model), ncol = ncol(X_df_model), byrow = TRUE)
        )
        names(newX) <- names(X_df_model)
        newdat <- cbind(newX, Z_df_model)

        s0 <- pmax(as.numeric(predict(vfit0, newdata = newdat)), eps)
        s1 <- pmax(as.numeric(predict(vfit1, newdata = newdat)), eps)
        list(s0 = s0, s1 = s1)
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 7. Estimate sigma_t(x) via local weighted variance in each arm
  # ---------------------------------------------------------------------------
  estimate_sigma_t_at_x <- function(x0_vec, mult = NULL) {
    out <- numeric(2)

    for (tval in 0:1) {
      idx <- which(T == tval)
      w <- kernel_weights_X(X_mat[idx, , drop = FALSE], x0_vec, hSigma)
      if (!is.null(mult)) w <- w * mult[idx]
      out[tval + 1] <- wvar(Y[idx], w)
    }

    names(out) <- c("sigma0", "sigma1")
    pmax(out, eps)
  }

  # ---------------------------------------------------------------------------
  # 8. Main estimation loop over x_eval
  # ---------------------------------------------------------------------------
  nxe <- nrow(x_eval_mat)
  rho_L_hat <- numeric(nxe)

  B_L_hat <- if (compute_B_L) numeric(nxe) else NULL
  rho_tilde_L_hat <- if (compute_B_L) numeric(nxe) else NULL

  for (j in seq_len(nxe)) {
    x0 <- x_eval_mat[j, ]
    wA <- kernel_weights_X(X_mat, x0, hA)

    mu_hat <- predict_mu_at_x(x0)
    Aj <- wcov(mu_hat$mu0, mu_hat$mu1, wA)

    sigx <- estimate_sigma_t_at_x(x0)
    denom <- max(sqrt(sigx["sigma0"] * sigx["sigma1"]), eps)

    rho_L_hat[j] <- max(-1, min(1, Aj / denom))

    if (compute_B_L) {
      if (is_binary) {
        s0 <- pmax(mu_hat$mu0 * (1 - mu_hat$mu0), eps)
        s1 <- pmax(mu_hat$mu1 * (1 - mu_hat$mu1), eps)
      } else {
        sig_xz <- predict_sigma_at_x(x0)
        s0 <- sig_xz$s0
        s1 <- sig_xz$s1
      }

      BLj <- wmean(sqrt(s0 * s1), wA)
      B_L_hat[j] <- BLj
      rho_tilde_L_hat[j] <- max(-1, min(1, (Aj - BLj) / denom))
    }
  }

  # ---------------------------------------------------------------------------
  # 9. Output
  # ---------------------------------------------------------------------------
  out <- list(
    x_eval = x_eval,
    mu_method = mu_method,
    rho_L_hat = rho_L_hat,
    rho_L_mean = mean(rho_L_hat, na.rm = TRUE),
    rho_L_var = var(rho_L_hat, na.rm = TRUE)
  )

  if (compute_B_L) {
    out$B_L_hat <- B_L_hat
    out$rho_tilde_L_hat <- rho_tilde_L_hat
    out$rho_tilde_L_mean <- mean(rho_tilde_L_hat, na.rm = TRUE)
    out$rho_tilde_L_var <- var(rho_tilde_L_hat, na.rm = TRUE)
  }

  # ---------------------------------------------------------------------------
  # 10. Bootstrap confidence intervals
  # ---------------------------------------------------------------------------
  if (ci) {
    alpha <- (1 - conf_level) / 2
    fisher_z <- function(r) atanh(pmin(pmax(r, -0.999999), 0.999999))
    inv_fisher_z <- function(z) tanh(z)

    if (!is.null(seed)) set.seed(seed)

    if (ci_method == "bootstrap_full") {
      idx0 <- which(T == 0L)
      idx1 <- which(T == 1L)

      rho_boot <- matrix(NA_real_, nrow = n_boot, ncol = nxe)
      rho_tilde_boot <- if (compute_B_L) matrix(NA_real_, nrow = n_boot, ncol = nxe) else NULL

      for (b in seq_len(n_boot)) {
        b0 <- sample(idx0, length(idx0), replace = TRUE)
        b1 <- sample(idx1, length(idx1), replace = TRUE)
        dat_b <- data[c(b0, b1), , drop = FALSE]

        est_b <- estimate_rho_L(
          data = dat_b,
          y_col = y_col,
          t_col = t_col,
          x_cols = x_cols,
          z_cols = z_cols,
          x_eval = x_eval,
          mu_method = mu_method,
          compute_B_L = compute_B_L,
          hA = hA,
          hSigma = hSigma,
          ranger_num_trees = ranger_num_trees,
          ci = FALSE,
          seed = if (is.null(seed)) NULL else (as.integer(seed) + 10000L + b),
          ranger_seed = ranger_seed
        )

        rho_boot[b, ] <- est_b$rho_L_hat
        if (compute_B_L) rho_tilde_boot[b, ] <- est_b$rho_tilde_L_hat
      }

      zmat <- apply(rho_boot, 2, fisher_z)
      z_lo <- apply(zmat, 2, quantile, probs = alpha, na.rm = TRUE, names = FALSE)
      z_hi <- apply(zmat, 2, quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)
      out$rho_L_ci <- cbind(lower = inv_fisher_z(z_lo), upper = inv_fisher_z(z_hi))

      if (compute_B_L) {
        zmat2 <- apply(rho_tilde_boot, 2, fisher_z)
        z2_lo <- apply(zmat2, 2, quantile, probs = alpha, na.rm = TRUE, names = FALSE)
        z2_hi <- apply(zmat2, 2, quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)
        out$rho_tilde_L_ci <- cbind(lower = inv_fisher_z(z2_lo), upper = inv_fisher_z(z2_hi))
      }

      if (store_boot) {
        out$rho_L_boot <- rho_boot
        if (compute_B_L) out$rho_tilde_L_boot <- rho_tilde_boot
      }
    }

    if (ci_method == "bootstrap_fast") {
      n <- nrow(data)
      rho_boot <- matrix(NA_real_, nrow = n_boot, ncol = nxe)
      rho_tilde_boot <- if (compute_B_L) matrix(NA_real_, nrow = n_boot, ncol = nxe) else NULL

      for (b in seq_len(n_boot)) {
        g <- rexp(n, rate = 1)
        g <- g / mean(g)

        for (j in seq_len(nxe)) {
          x0 <- x_eval_mat[j, ]
          wA <- kernel_weights_X(X_mat, x0, hA) * g

          mu_hat <- predict_mu_at_x(x0)
          Aj <- wcov(mu_hat$mu0, mu_hat$mu1, wA)

          sigx <- estimate_sigma_t_at_x(x0, mult = g)
          denom <- max(sqrt(sigx["sigma0"] * sigx["sigma1"]), eps)

          rho_boot[b, j] <- max(-1, min(1, Aj / denom))

          if (compute_B_L) {
            if (is_binary) {
              s0 <- pmax(mu_hat$mu0 * (1 - mu_hat$mu0), eps)
              s1 <- pmax(mu_hat$mu1 * (1 - mu_hat$mu1), eps)
            } else {
              sig_xz <- predict_sigma_at_x(x0)
              s0 <- sig_xz$s0
              s1 <- sig_xz$s1
            }

            BLj <- wmean(sqrt(s0 * s1), wA)
            rho_tilde_boot[b, j] <- max(-1, min(1, (Aj - BLj) / denom))
          }
        }
      }

      zmat <- apply(rho_boot, 2, fisher_z)
      z_lo <- apply(zmat, 2, quantile, probs = alpha, na.rm = TRUE, names = FALSE)
      z_hi <- apply(zmat, 2, quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)
      out$rho_L_ci <- cbind(lower = inv_fisher_z(z_lo), upper = inv_fisher_z(z_hi))

      if (compute_B_L) {
        zmat2 <- apply(rho_tilde_boot, 2, fisher_z)
        z2_lo <- apply(zmat2, 2, quantile, probs = alpha, na.rm = TRUE, names = FALSE)
        z2_hi <- apply(zmat2, 2, quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)
        out$rho_tilde_L_ci <- cbind(lower = inv_fisher_z(z2_lo), upper = inv_fisher_z(z2_hi))
      }

      if (store_boot) {
        out$rho_L_boot <- rho_boot
        if (compute_B_L) out$rho_tilde_L_boot <- rho_tilde_boot
      }
    }
  }

  return(out)
}
