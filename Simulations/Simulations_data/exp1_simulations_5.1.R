###############################################################################################
# Experiment 1 (Section 5.1 and Appendix C): synthetic settings (A) d = 1, (B) d = 5, (C) d = 10.
# Reproduces Figure fig_sim1 (sim_main_compact.pdf), the Appendix C figures sim_coverage_boxplots.pdf,
# sim_cov_by_propensity.pdf and sim_cov_width_loss.pdf, the additive-model check, and the numbers quoted in the text.
# Methods: CW(rho), CW(rho - 0.25), CW^{+CI}_K(rho) with K = K_hat (Appendix 'Choice of K'), CMC, DR conformal,
# the inexact and naive intervals of Lei and Candes, and the oracle interval.
# Self-contained: needs no other file. Packages: mgcv, qgam, randomForest, quantregForest, gbm, MASS, ggplot2, patchwork.
# Runtime on one core: about 1 h for 3 x 20 replications, plus about 25 min for the additive check.
# Usage: Rscript exp1_simulations_5.1.R      (resumable; outputs in OUT_DIR)
###############################################################################################
OUT_DIR      <- "results_exp1"                          # cache, results (CSV), figures
RUN_ADDITIVE <- TRUE                                    # Appendix C: setting (B) with purely additive GAMs
suppressPackageStartupMessages({
  library(mgcv); library(qgam); library(randomForest); library(quantregForest); library(gbm)
  library(MASS); library(ggplot2); library(patchwork)
})

ALPHA    <- 0.1
RHO_GRID <- seq(-1, 1, by = 0.25)
N_TRAIN  <- 3000; N_TEST <- 2000; N_REPS <- 20; N_BOOT <- 50
DIMS     <- c(A = 1, B = 5, C = 10)
# GAM model structure used when d <= 5 (te(X1, X2) in (B), where tau is not additive)
GAM_RHS  <- c(A = "s(X1)", B = "te(X1, X2) + s(X3) + s(X4) + s(X5)", C = NA)
N_REPS <- as.integer(Sys.getenv("N_REPS", N_REPS))         # 20 in the paper; e.g. N_REPS=1 for a quick test

# ------------------------------------------------------------------------------------------------
# CW^{+CI}_K(rho): D_{rho,K}(a,b) = max{D_rho(a,b), lambda a, lambda b},
#   lambda^2 = 1 - 1/K^2 - (2/K)(rho - 1/K)_+  (lambda = 1 if K = Inf, 0 if K = 1),  c = c*(x), hence
#   U = tau + r_u + D_{rho,K}((l_0 - r_u/2)_+, (u_1 - r_u/2)_+),  L = tau - r_l - D_{rho,K}((l_1 - r_l/2)_+, (u_0 - r_l/2)_+)
# ------------------------------------------------------------------------------------------------
lambda_rhoK <- function(rho, K) ifelse(is.infinite(K), 1, sqrt(pmax(1 - 1 / K^2 - (2 / K) * pmax(rho - 1 / K, 0), 0)))
D_rhoK      <- function(a, b, rho, K) pmax(rho_distance(a, b, rho), lambda_rhoK(rho, K) * pmax(a, b))
ci_half_K   <- function(a, b, r, rho, K) r + D_rhoK(pmax(a - r / 2, 0), pmax(b - r / 2, 0), rho, K)
cw_ci_interval <- function(o, rho, K = o$Khat)
  cbind(o$cate - ci_half_K(o$l1, o$u0, o$rl, rho, K), o$cate + ci_half_K(o$l0, o$u1, o$ru, rho, K))


# ------------------------------------------------------------------------------------------------
# Data-generating processes
# ------------------------------------------------------------------------------------------------
varsigma <- function(u) 1 + 1 / (1 + exp(-20 * (u - 1/3)))

make_data <- function(setting, n_train = N_TRAIN, n_test = N_TEST) {
  N <- n_train + n_test; d <- DIMS[[setting]]
  if (d == 1) {
    X <- matrix(runif(N), ncol = 1)
  } else {
    Sig <- matrix(0.25, d, d); diag(Sig) <- 1
    X <- pnorm(MASS::mvrnorm(N, rep(0, d), Sig))
  }
  colnames(X) <- paste0("X", seq_len(d))
  x1 <- X[, 1]; x2 <- if (d > 1) X[, 2] else X[, 1]
  tau   <- varsigma(x1) * varsigma(x2)
  beta  <- runif(d, -1, 1)
  kappa <- runif(1, 1, 5)
  mu0 <- as.numeric(X %*% beta); s0 <- rep(1, N)
  if (setting %in% c("A", "B")) {
    p1 <- (3 - x1) / 4                                       # P(T = 1 | X) in [0.5, 0.75]
  } else {
    p1  <- pmin(0.9, pmax(0.1, sin(pi * x1 * x2)))           # Nie & Wager (2021), Setup A
    mu0 <- mu0 + 2 * sin(pi * x1 * x2) + 4 * (X[, 3] - 0.5)^2
    s0  <- s0 * (0.5 + X[, 3])
  }
  s1 <- kappa * s0; mu1 <- mu0 + tau
  T <- rbinom(N, 1, p1); z <- rnorm(N)
  Y <- ifelse(T == 1, mu1 + s1 * z, mu0 + s0 * z)            # observed data: law does not depend on rho
  tr <- seq_len(n_train); te <- n_train + seq_len(n_test)
  list(X = X[tr, , drop = FALSE], T = T[tr], Y = Y[tr], Xte = X[te, , drop = FALSE],
       tau = tau[te], s0 = s0[te], s1 = s1[te], p1 = p1[te],
       Z0 = rnorm(n_test), Z1 = rnorm(n_test), kappa = kappa)
}

ite_given_rho <- function(dat, rho)
  dat$tau + dat$s1 * (rho * dat$Z0 + sqrt(max(1 - rho^2, 0)) * dat$Z1) - dat$s0 * dat$Z0

# ------------------------------------------------------------------------------------------------
# Learners: GAMs (d <= 5) and random forests (d > 5)
# ------------------------------------------------------------------------------------------------
fit_quantiles <- function(Xtr, ytr, X_eval, setting, qs = c(ALPHA / 2, 0.5, 1 - ALPHA / 2)) {
  d <- ncol(Xtr)
  if (d > 5) {
    qrf <- quantregForest(x = Xtr, y = ytr, ntree = 300, nodesize = 5, mtry = max(floor(d / 3), 1))
    # matrix() keeps one column per quantile: predict() returns a plain vector when length(qs) == 1
    return(lapply(X_eval, function(Xe) matrix(predict(qrf, newdata = Xe, what = qs), ncol = length(qs))))
  }
  df <- data.frame(y = ytr, Xtr); form <- as.formula(paste("y ~", GAM_RHS[[setting]]))
  fits <- lapply(qs, function(q) qgam(form, data = df, qu = q))
  lapply(X_eval, function(Xe) sapply(fits, function(f) as.numeric(predict(f, newdata = data.frame(Xe)))))
}

fit_mean <- function(Xtr, ytr, setting, ntree = 300) {            # classical mean regression
  d <- ncol(Xtr)
  if (d > 5) {
    rf <- randomForest(x = Xtr, y = ytr, ntree = ntree, nodesize = 5, mtry = max(floor(d / 3), 1))
    return(function(Xn) as.numeric(predict(rf, newdata = Xn)))
  }
  g <- gam(as.formula(paste("y ~", GAM_RHS[[setting]])), data = data.frame(y = ytr, Xtr), method = "GCV.Cp")
  function(Xn) as.numeric(predict(g, newdata = data.frame(Xn)))
}

gbm_propensity <- function(X, T) {
  fit <- gbm(T ~ ., data = data.frame(T = T, X), distribution = "bernoulli", n.trees = 300,
             interaction.depth = 3, shrinkage = 0.05, verbose = FALSE)
  function(Xn) predict(fit, newdata = data.frame(Xn), n.trees = 300, type = "response")
}

wquantile <- function(x, w, p) {
  o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w)
  x[min(which(cw >= p)[1], length(x), na.rm = TRUE)]
}

# ------------------------------------------------------------------------------------------------
# CW(rho) and CW+CI(rho): building blocks fitted once, combined for every rho
# ------------------------------------------------------------------------------------------------
rho_distance <- function(x, y, rho) sqrt(pmax(x^2 + y^2 - 2 * rho * x * y, 0))
CI_addition  <- function(r, a, b, rho) {                          # exact c* r (as in CW_rho_intervals)
  v_a <- pmax(a - r / 2, 0); v_b <- pmax(b - r / 2, 0)
  ifelse(r == 0, 0, pmax(0, rho_distance(v_a, v_b, rho) + (a - v_a) + (b - v_b) - rho_distance(a, b, rho)))
}


cw_fit <- function(X, T, Y, Xte, setting, seed, B = N_BOOT, clip = 0.2) {
  idx <- list(which(T == 0), which(T == 1))
  tr  <- lapply(idx, function(i) i[seq_len(floor(0.8 * length(i)))])
  ca  <- lapply(seq_along(idx), function(k) setdiff(idx[[k]], tr[[k]]))
  ps  <- gbm_propensity(X[c(tr[[1]], tr[[2]]), , drop = FALSE], T[c(tr[[1]], tr[[2]])])
  out <- list(); mfun <- list()
  for (t in 0:1) {
    k <- t + 1
    q  <- fit_quantiles(X[tr[[k]], , drop = FALSE], Y[tr[[k]]],
                        list(X[ca[[k]], , drop = FALSE], Xte), setting)
    fix <- function(m) { m[, 3] <- pmax(m[, 3], m[, 2]); m[, 1] <- pmin(m[, 1], m[, 2]); m }
    qc <- fix(q[[1]]); qt <- fix(q[[2]]); yc <- Y[ca[[k]]]
    sc <- pmax((qc[, 2] - yc) / pmax(qc[, 2] - qc[, 1], 1e-12),          # multiplicative CQR score
               (yc - qc[, 2]) / pmax(qc[, 3] - qc[, 2], 1e-12))
    e  <- pmin(pmax(ps(X[ca[[k]], , drop = FALSE]), clip), 1 - clip)
    g  <- wquantile(sc, if (t == 1) 1 / e else 1 / (1 - e), 1 - ALPHA)
    lower <- qt[, 2] - g * (qt[, 2] - qt[, 1]); upper <- qt[, 2] + g * (qt[, 3] - qt[, 2])
    mfun[[k]] <- fit_mean(X[tr[[k]], , drop = FALSE], Y[tr[[k]]], setting); f <- mfun[[k]](Xte)
    out[[paste0("f", t)]] <- f
    out[[paste0("l", t)]] <- pmax(f - lower, 0); out[[paste0("u", t)]] <- pmax(upper - f, 0)
  }
  ## K_hat (Appendix 'Choice of K'): sigma_t = regression of |Y - mu_hat_t(X)| on X in the calibration split of
  ## arm t (learner of the mean model, floored at 0.05 mean|residual|); K_hat = 90% quantile over the test points
  ## of max{kappa, 1/kappa}, kappa = sigma_1/sigma_0. Own random-number stream: the main stream is left unchanged.
  old_seed <- get(".Random.seed", envir = globalenv()); set.seed(seed + 777)
  sig <- lapply(1:2, function(k) { Xc <- X[ca[[k]], , drop = FALSE]; r <- abs(Y[ca[[k]]] - mfun[[k]](Xc))
                                   pmax(fit_mean(Xc, r, setting)(Xte), 0.05 * mean(r)) })
  assign(".Random.seed", old_seed, envir = globalenv())
  kap <- sig[[2]] / sig[[1]]; out$Khat <- max(1, quantile(pmax(kap, 1 / kap), 0.9, names = FALSE))
  boot <- list(matrix(NA, B, nrow(Xte)), matrix(NA, B, nrow(Xte))); n <- length(Y)
  for (b in seq_len(B)) {                                           # bootstrap CI for mu_0 and mu_1
    bi <- sample.int(n, n, replace = TRUE)
    for (t in 0:1) {
      rows <- bi[T[bi] == t]; rows <- rows[seq_len(floor(0.8 * length(rows)))]
      boot[[t + 1]][b, ] <- fit_mean(X[rows, , drop = FALSE], Y[rows], setting, ntree = 50)(Xte)
    }
  }
  ci <- lapply(boot, function(M) apply(M, 2, quantile, probs = c(0.05, 0.95)))
  out$rl   <- 2 * pmax(out$f1 - ci[[2]][1, ], ci[[1]][2, ] - out$f0, 0)
  out$ru   <- 2 * pmax(out$f0 - ci[[1]][1, ], ci[[2]][2, ] - out$f1, 0)
  out$cate <- out$f1 - out$f0
  out
}


cw_interval <- function(o, rho, ci = FALSE) {
  lo <- o$cate - rho_distance(o$l1, o$u0, rho); hi <- o$cate + rho_distance(o$l0, o$u1, rho)
  if (ci) {
    lo <- lo - CI_addition(o$rl, o$l1, o$u0, rho); hi <- hi + CI_addition(o$ru, o$l0, o$u1, rho)
  }
  cbind(lo, hi)
}

# ------------------------------------------------------------------------------------------------
# Baselines, all with the SAME nuisance estimators as CW: fit_mean / fit_quantiles (GAMs when d <= 5,
# random forests when d > 5) and the gbm propensity score
# ------------------------------------------------------------------------------------------------
wcp_threshold <- function(scores, w_cal, w_test, level) {           # weighted split conformal quantile
  o <- order(scores); s <- scores[o]; cw <- cumsum(w_cal[o])
  k <- findInterval(level * (cw[length(cw)] + w_test), cw, left.open = TRUE) + 1
  s[pmin(k, length(s))]                                             # cfcausal useInf = FALSE
}
clipw <- function(w) pmin(pmax(w, 0.05), 20)

lei_naive <- function(X, T, Y, Xte, setting) {                       # Lei & Candes (2021), naive (Bonferroni)
  n <- length(Y); perm <- sample.int(n); ntr <- floor(0.75 * n)
  itr <- perm[seq_len(ntr)]; ica <- perm[-seq_len(ntr)]
  ps <- gbm_propensity(X[itr, , drop = FALSE], T[itr]); ec <- ps(X[ica, , drop = FALSE]); et <- ps(Xte)
  b <- list()
  for (t in 0:1) {
    atr <- itr[T[itr] == t]; cal <- ica[T[ica] == t]; e_c <- ec[T[ica] == t]
    q <- fit_quantiles(X[atr, , drop = FALSE], Y[atr], list(X[cal, , drop = FALSE], Xte), setting, qs = c(0.05, 0.95))
    sc  <- pmax(q[[1]][, 1] - Y[cal], Y[cal] - q[[1]][, 2])           # CQR score
    thr <- wcp_threshold(sc, clipw(if (t == 1) 1 / e_c else 1 / (1 - e_c)),
                         clipw(if (t == 1) 1 / et else 1 / (1 - et)), 1 - ALPHA / 2)
    b[[t + 1]] <- cbind(q[[2]][, 1] - thr, q[[2]][, 2] + thr)
  }
  cbind(b[[2]][, 1] - b[[1]][, 2], b[[2]][, 2] - b[[1]][, 1])
}

lei_inexact <- function(X, T, Y, Xte, setting) {                     # inexact nested (cfcausal defaults)
  f2 <- unlist(lapply(0:1, function(t) { i <- which(T == t); i[sample.int(length(i), ceiling(length(i) / 2))] }))
  f1 <- sample(setdiff(seq_along(Y), f2)); ntr <- floor(0.75 * length(f1))
  itr <- f1[seq_len(ntr)]; ica <- f1[-seq_len(ntr)]
  ps <- gbm_propensity(X[itr, , drop = FALSE], T[itr]); ec <- ps(X[ica, , drop = FALSE]); e2 <- ps(X[f2, , drop = FALSE])
  lo2 <- hi2 <- numeric(length(f2))
  for (t in 0:1) {                                                  # interval for the missing Y(t)
    atr <- itr[T[itr] == t]; cal <- ica[T[ica] == t]; tgt <- which(T[f2] == 1 - t)
    q <- fit_quantiles(X[atr, , drop = FALSE], Y[atr], list(X[cal, , drop = FALSE], X[f2[tgt], , drop = FALSE]),
                       setting, qs = c(0.05, 0.95))
    sc <- pmax(q[[1]][, 1] - Y[cal], Y[cal] - q[[1]][, 2])
    r  <- if (t == 0) function(e) e / (1 - e) else function(e) (1 - e) / e
    thr <- wcp_threshold(sc, clipw(r(ec[T[ica] == t])), clipw(r(e2[tgt])), 1 - ALPHA)
    cl <- q[[2]][, 1] - thr; cu <- q[[2]][, 2] + thr; y2 <- Y[f2[tgt]]
    if (t == 0) { lo2[tgt] <- y2 - cu; hi2[tgt] <- y2 - cl } else { lo2[tgt] <- cl - y2; hi2[tgt] <- cu - y2 }
  }
  Xf2 <- X[f2, , drop = FALSE]
  cbind(fit_quantiles(Xf2, lo2, list(Xte), setting, qs = 0.4)[[1]][, 1],
        fit_quantiles(Xf2, hi2, list(Xte), setting, qs = 0.6)[[1]][, 1])
}

cps_draw <- function(yhat, r, w, w_test, u) {                       # weighted conformal predictive distribution
  o <- order(r); r <- r[o]; cw <- cumsum(w[o])
  k <- findInterval(u * (cw[length(cw)] + w_test), cw, left.open = TRUE) + 1
  yhat + r[pmin(k, length(r))]
}

cmc_intervals <- function(X, T, Y, Xte, p, setting, mc = 100) {     # CMC T-learner with pseudo Monte Carlo
  n <- length(Y)                                                    # (as in the reference implementation)
  i_nu <- sample.int(n, floor(0.75 * n)); i_cal <- setdiff(seq_len(n), i_nu)
  i_pt <- sample(i_nu, floor(0.75 * length(i_nu))); i_cn <- setdiff(i_nu, i_pt)
  m0 <- fit_mean(X[i_pt[T[i_pt] == 0], , drop = FALSE], Y[i_pt[T[i_pt] == 0]], setting)
  m1 <- fit_mean(X[i_pt[T[i_pt] == 1], , drop = FALSE], Y[i_pt[T[i_pt] == 1]], setting)
  c0 <- i_cn[T[i_cn] == 0]; c1 <- i_cn[T[i_cn] == 1]
  r0 <- Y[c0] - m0(X[c0, , drop = FALSE]); w0 <- p[c0] / (1 - p[c0])
  r1 <- Y[c1] - m1(X[c1, , drop = FALSE]); w1 <- rep(1, length(c1))
  Xc <- X[i_cal, , drop = FALSE]; y0c <- m0(Xc); y1c <- m1(Xc)
  ite <- t(sapply(seq_along(i_cal), function(j) {                   # pseudo-ITEs of the calibration units
    i <- i_cal[j]; u <- runif(mc)
    if (T[i] == 0) cps_draw(y1c[j], r1, w1, (1 - p[i]) / p[i], u) - Y[i]
    else           Y[i] - cps_draw(y0c[j], r0, w0, p[i] / (1 - p[i]), u)
  }))
  q <- quantile(as.vector(ite - (y1c - y0c)), c(0.05, 0.95))        # predictive distribution of the ITE
  tau_te <- m1(Xte) - m0(Xte)
  cbind(tau_te + q[1], tau_te + q[2])
}

dr_pseudo <- function(W, p, Y, y0, y1) (W - p) / (p * (1 - p)) * (Y - ifelse(W == 1, y1, y0)) + (y1 - y0)

dr_intervals <- function(X, T, Y, Xte, p, setting, n_folds = 5) {   # conformal DR meta-learner
  n <- length(Y); perm <- sample.int(n); ntr <- floor(0.75 * n)
  itr <- perm[seq_len(ntr)]; ica <- perm[-seq_len(ntr)]
  folds <- split(itr, rep(seq_len(n_folds), length.out = length(itr)))
  m0s <- m1s <- mms <- vector("list", n_folds)
  for (k in seq_len(n_folds)) {                                     # cross-fitting
    te_k <- folds[[k]]; tr_k <- setdiff(itr, te_k)
    a0 <- tr_k[T[tr_k] == 0]; a1 <- tr_k[T[tr_k] == 1]
    m0s[[k]] <- fit_mean(X[a0, , drop = FALSE], Y[a0], setting)
    m1s[[k]] <- fit_mean(X[a1, , drop = FALSE], Y[a1], setting)
    Xk <- X[te_k, , drop = FALSE]
    mms[[k]] <- fit_mean(Xk, dr_pseudo(T[te_k], p[te_k], Y[te_k], m0s[[k]](Xk), m1s[[k]](Xk)), setting)
  }
  avg <- function(ms, Z) Reduce(`+`, lapply(ms, function(m) m(Z))) / length(ms)
  Xc <- X[ica, , drop = FALSE]
  yc <- dr_pseudo(T[ica], p[ica], Y[ica], avg(m0s, Xc), avg(m1s, Xc))
  res <- sort(abs(avg(mms, Xc) - yc)); nc <- length(res)
  off <- res[min(ceiling((1 - ALPHA) * (1 + 1 / nc) * nc), nc)]     # held-out calibration
  th <- avg(mms, Xte); cbind(th - off, th + off)
}

# ------------------------------------------------------------------------------------------------
# Metrics and one replication
# ------------------------------------------------------------------------------------------------
interval_score <- function(lo, hi, y)                               # Gneiting & Raftery (2007)
  mean((hi - lo) + (2 / ALPHA) * pmax(lo - y, 0) + (2 / ALPHA) * pmax(y - hi, 0))


run_replication <- function(setting, rep) {
  seed <- 10000 * match(setting, names(DIMS)) + rep; set.seed(seed)
  dat <- make_data(setting); X <- dat$X; T <- dat$T; Y <- dat$Y; Xte <- dat$Xte
  o <- cw_fit(X, T, Y, Xte, setting, seed)
  p_shared <- pmin(pmax(gbm_propensity(X, T)(X), 0.05), 0.95)
  fits <- list(Lei_naive   = tryCatch(lei_naive(X, T, Y, Xte, setting),   error = function(e) NULL),
               Lei_inexact = tryCatch(lei_inexact(X, T, Y, Xte, setting), error = function(e) NULL),
               CMC = tryCatch(cmc_intervals(X, T, Y, Xte, p_shared, setting), error = function(e) NULL),
               DR  = tryCatch(dr_intervals(X, T, Y, Xte, p_shared, setting), error = function(e) NULL))
  fits <- fits[!vapply(fits, is.null, logical(1))]
  ## building blocks: everything needed to (re)evaluate any combination rule without refitting
  list(setting = setting, rep = rep, seed = seed, o = o, fits = fits,
       test = dat[c("tau", "s0", "s1", "p1", "Z0", "Z1", "kappa")])
}

## CW^{+CI} from the cached blocks; version "new" = revised Section section_CI, "old" = c* r rule

evaluate_replication <- function(B, label = B$setting) {     # metrics of every method at every rho
  o <- B$o; dat <- B$test
  bins <- findInterval(dat$p1, quantile(dat$p1, c(0.2, 0.4, 0.6, 0.8)))   # propensity quintiles 0..4
  zq <- qnorm(1 - ALPHA / 2); rows <- list()
  for (rho in RHO_GRID) {
    ite <- ite_given_rho(dat, rho)
    ivs <- c(list(CW = cw_interval(o, rho), CW_misspec = cw_interval(o, max(rho - 0.25, -1)),
                  CW_CI = cw_ci_interval(o, rho)), B$fits)
    ow <- 2 * zq * rho_distance(dat$s0, dat$s1, rho)
    ivs$Oracle <- cbind(dat$tau - ow / 2, dat$tau + ow / 2)          # shortest 90% interval (Gaussian)
    for (m in names(ivs)) {
      lo <- ivs[[m]][, 1]; hi <- ivs[[m]][, 2]; cov <- (lo < ite) & (ite < hi)
      q <- sapply(0:4, function(k) mean(cov[bins == k]))
      rows[[length(rows) + 1]] <- data.frame(
        setting = label, rep = B$rep, rho = rho, method = m, kappa = dat$kappa, coverage = mean(cov),
        width = median(hi - lo), interval_score = interval_score(lo, hi, ite),
        cov_q1 = q[1], cov_q2 = q[2], cov_q3 = q[3], cov_q4 = q[4], cov_q5 = q[5], Khat = o$Khat)
    }
  }
  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------------------------
# Run (resumable: each replication is cached in OUT_DIR/cache)
# ------------------------------------------------------------------------------------------------
dir.create(file.path(OUT_DIR, "cache"), recursive = TRUE, showWarnings = FALSE)
FIG <- file.path(OUT_DIR, "figures"); dir.create(FIG, showWarnings = FALSE)
cache_file <- function(label, rep) file.path(OUT_DIR, "cache", sprintf("%s_rep%02d.rds", label, rep))
run_cached <- function(label, setting, rep) {
  f <- cache_file(label, rep)
  if (!file.exists(f)) {
    t0 <- Sys.time(); saveRDS(run_replication(setting, rep), f)
    cat(sprintf("%-4s rep %2d  %.1f s\n", label, rep, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  }
  invisible(f)
}
load_res <- function(label) do.call(rbind, lapply(seq_len(N_REPS), function(rep) evaluate_replication(readRDS(cache_file(label, rep)), label)))
for (rep in seq_len(N_REPS)) for (s in c("A", "B", "C")) run_cached(s, s, rep)
if (RUN_ADDITIVE) {                                      # same data and seeds as (B); additive GAMs only
  rhs_B <- GAM_RHS[["B"]]; GAM_RHS[["B"]] <- "s(X1) + s(X2) + s(X3) + s(X4) + s(X5)"
  for (rep in seq_len(N_REPS)) run_cached("Badd", "B", rep)
  GAM_RHS[["B"]] <- rhs_B
}
res <- do.call(rbind, lapply(c("A", "B", "C"), load_res))
write.csv(res, file.path(OUT_DIR, "results_A_B_C.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------------------------------------
METHODS <- c("CW", "CW_misspec", "CW_CI", "CMC", "DR", "Lei_inexact", "Lei_naive")
LABELS  <- list(CW = quote(CW(rho)), CW_misspec = quote(CW(rho - 0.25)), CW_CI = quote(CW^{"+CI"}*(rho)),
                CMC = "CMC", DR = "DR", Lei_inexact = "Lei (inexact)", Lei_naive = "Lei (naive)",
                CW_CI_old = quote(CW^{"+CI"}*(rho)~"previous rule"), Oracle = "oracle")
COLORS  <- c(CW = "#e31a1c", CW_misspec = "#fd8d3c", CW_CI = "#800026", CMC = "#969696", DR = "#9ecae1",
             Lei_inexact = "#a1d99b", Lei_naive = "#6a51a3", CW_CI_old = "#800026", Oracle = "black")
LTYPES  <- c(CW = "solid", CW_misspec = "dashed", CW_CI = "solid", CMC = "solid", DR = "solid",
             Lei_inexact = "solid", Lei_naive = "dotted", CW_CI_old = "22", Oracle = "dashed")
SHAPES  <- c(CW = 16, CW_misspec = 15, CW_CI = 18, CMC = 17, DR = 16, Lei_inexact = 16, Lei_naive = 3,
             CW_CI_old = 5, Oracle = NA)
LWD     <- c(CW = 0.9, CW_misspec = 0.9, CW_CI = 0.9, CMC = 0.45, DR = 0.45, Lei_inexact = 0.45,
             Lei_naive = 0.6, CW_CI_old = 0.7, Oracle = 0.5)
PSIZE   <- c(CW = 1.3, CW_misspec = 1.2, CW_CI = 1.8, CMC = 1, DR = 0.9, Lei_inexact = 0.9, Lei_naive = 1.2,
             CW_CI_old = 1.3, Oracle = 0)
SET_LAB <- c(A = "(A)~~italic(d)==1", B = "(B)~~italic(d)==5",
             C = "atop((C)~~italic(d)==10*','~complex~pi(x)*',',nonlinear*','~heteroskedastic)")
lab_expr <- function(ms) do.call(c, lapply(ms, function(m) as.expression(LABELS[[m]])))
scales_for <- function(ms) list(
  scale_colour_manual(values = COLORS[ms], labels = lab_expr(ms), breaks = ms, limits = ms, name = NULL),
  scale_linetype_manual(values = LTYPES[ms], labels = lab_expr(ms), breaks = ms, limits = ms, name = NULL),
  scale_shape_manual(values = SHAPES[ms], labels = lab_expr(ms), breaks = ms, limits = ms, name = NULL, na.value = NA),
  scale_discrete_manual("linewidth", values = LWD[ms], labels = lab_expr(ms), breaks = ms, limits = ms, name = NULL),
  scale_size_manual(values = PSIZE[ms], labels = lab_expr(ms), breaks = ms, limits = ms, name = NULL))
## lines and points: size mapped on points only; the oracle is drawn as a line without points
geoms_lp <- function() list(geom_line(aes(linewidth = method)),
                            geom_point(aes(size = method), data = function(d) d[d$method != "Oracle", ]))
AES <- aes(colour = method, linetype = method, shape = method)
one_row <- function(n = 1) guides(colour = guide_legend(nrow = n, byrow = TRUE), linetype = guide_legend(nrow = n, byrow = TRUE),
  shape = guide_legend(nrow = n, byrow = TRUE), linewidth = guide_legend(nrow = n, byrow = TRUE), size = guide_legend(nrow = n, byrow = TRUE))
TH <- theme_bw(base_size = 8.5) + theme(legend.position = "bottom", legend.key.width = unit(1.4, "lines"),
  panel.grid.minor = element_blank(), strip.background = element_blank(), legend.margin = margin(0, 0, 0, 0),
  strip.text = element_text(size = 8))
xrho <- scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("-1", "-.5", "0", ".5", "1"))
XLAB <- expression(paste("cross-world correlation ", rho))

## ratios of averages over replications, relative to the oracle (as in the paper)
summarise_rel <- function(r) {
  S <- aggregate(cbind(coverage, width, interval_score) ~ setting + method + rho, data = r, FUN = mean)
  O <- S[S$method == "Oracle", c("setting", "rho", "width", "interval_score")]
  names(O)[3:4] <- c("width_or", "is_or")
  S <- merge(S[S$method != "Oracle", ], O, by = c("setting", "rho"))
  S$width_rel <- S$width / S$width_or; S$is_rel <- S$interval_score / S$is_or
  S
}
S_all <- summarise_rel(res)

line_panel <- function(S, y, ylab, ref, ms, log = FALSE, ylim = NULL, free = FALSE) {
  S <- S[S$method %in% ms, ]; S$method <- factor(S$method, levels = rev(ms))
  p <- ggplot(S, aes(rho, .data[[y]])) + AES +
    geom_hline(yintercept = ref, colour = "black", linewidth = 0.35, linetype = if (log) "dashed" else "solid") +
    geoms_lp() + facet_wrap(~setting, nrow = 1, labeller = as_labeller(SET_LAB, label_parsed), scales = if (free) "free_y" else "fixed") +
    scales_for(ms) + xrho + labs(x = XLAB, y = ylab) + TH + one_row(ceiling(length(ms) / 7))
  if (log) p <- p + scale_y_continuous(trans = "log2", breaks = c(0.5, 0.75, 1, 1.25, 1.5, 1.75, 2))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}
cov_lo <- function(S, ms) min(0.65, floor(20 * min(S$coverage[S$method %in% ms])) / 20)


## ---------------------------------------------------------------- Figure fig_sim1 (main text)
rel_lo <- min(0.95, 0.95 * min(S_all$width_rel[S_all$method %in% METHODS]))
f1 <- (line_panel(S_all, "is_rel", "interval score", 1, METHODS, log = TRUE, ylim = c(0.95, 2)) +
         theme(axis.title.x = element_blank(), axis.text.x = element_blank())) /
      (line_panel(S_all, "coverage", "coverage", 1 - ALPHA, METHODS, ylim = c(cov_lo(S_all, METHODS), 1)) +
         theme(axis.title.x = element_blank(), axis.text.x = element_blank(), strip.text = element_blank())) /
      (line_panel(S_all, "width_rel", "median width", 1, METHODS, log = TRUE, ylim = c(rel_lo, 2)) +
         theme(strip.text = element_blank())) +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(FIG, "sim_main_compact.pdf"), f1, width = 6.3, height = 5.9)

## ---------------------------------------------------------------- Figure fig_sim_boxplots
B <- res[res$method %in% METHODS & res$rho %in% c(-1, -0.5, 0, 0.5, 1), ]
B$method <- factor(B$method, levels = rev(METHODS))
B$rho_lab <- factor(sprintf("rho == %s", B$rho), levels = sprintf("rho == %s", c(-1, -0.5, 0, 0.5, 1)))
B$set_lab <- factor(B$setting, levels = c("A", "B", "C"), labels = c("(A)~d==1", "(B)~d==5", "(C)~d==10"))
f2 <- ggplot(B, aes(coverage, method, fill = method)) +
  geom_vline(xintercept = 1 - ALPHA, linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.5, linewidth = 0.25, width = 0.7) +
  facet_grid(set_lab ~ rho_lab, labeller = label_parsed) +
  scale_fill_manual(values = COLORS, guide = "none") +
  scale_y_discrete(labels = function(b) lab_expr(b)) +
  scale_x_continuous(breaks = c(0.6, 0.8, 1), labels = c("0.6", "0.8", "1.0")) +
  labs(x = "coverage across replications", y = NULL) + TH
ggsave(file.path(FIG, "sim_coverage_boxplots.pdf"), f2, width = 6.5, height = 3.2 + 0.2 * length(METHODS))

## ---------------------------------------------------------------- Figure fig_sim_propensity (setting C)
MP <- c("CW", "CW_CI", "CMC", "DR", "Lei_inexact")
P <- res[res$setting == "C" & res$rho %in% c(-0.5, 0, 0.5) & res$method %in% MP, ]
P <- do.call(rbind, lapply(1:5, function(k) data.frame(P[, c("rho", "method")], q = k, cov = P[[paste0("cov_q", k)]])))
P <- aggregate(cov ~ rho + method + q, data = P, FUN = mean)
P$method <- factor(P$method, levels = rev(MP))
P$rho_lab <- factor(sprintf("(C)~~rho == %s", P$rho), levels = sprintf("(C)~~rho == %s", c(-0.5, 0, 0.5)))
f3 <- ggplot(P, aes(q, cov)) + AES +
  geom_hline(yintercept = 1 - ALPHA, linewidth = 0.35) + geoms_lp() +
  facet_wrap(~rho_lab, nrow = 1, labeller = label_parsed) + scales_for(MP) +
  labs(x = expression(paste("quintile of ", pi(x))), y = "coverage") + TH
ggsave(file.path(FIG, "sim_cov_by_propensity.pdf"), f3 + one_row(ceiling(length(MP) / 6)), width = 6.3, height = 2.7)

## ---------------------------------------------------------------- Figure coverage-width loss
L <- res[res$method %in% METHODS, ]
key <- interaction(L$setting, L$rep, L$rho, drop = TRUE)
wmin <- ave(L$width, key, FUN = min); wmax <- ave(L$width, key, FUN = max)
L$loss <- (L$width - wmin) / pmax(wmax - wmin, 1e-12) + (2 / ALPHA) * abs(L$coverage - (1 - ALPHA))
LS <- aggregate(loss ~ method + rho, data = L, FUN = mean)
LS$method <- factor(LS$method, levels = rev(METHODS))
f4 <- ggplot(LS, aes(rho, loss)) + AES + geoms_lp() + scales_for(METHODS) + xrho + labs(x = XLAB, y = "coverage-width loss") + TH +
  theme(legend.position = "right")
ggsave(file.path(FIG, "sim_cov_width_loss.pdf"), f4, width = 6.3, height = 2.6)


# ------------------------------------------------------------------------------------------------
# Numbers quoted in Section 5.1 and Appendix C
# ------------------------------------------------------------------------------------------------
sink(file.path(OUT_DIR, "numbers.txt"), split = TRUE)
g <- function(s, m, col) { d <- S_all[S_all$setting == s & S_all$method == m, ]; d[order(d$rho), col] }
cat(sprintf("rho grid: %s\n", paste(RHO_GRID, collapse = " ")))
for (s in c("A", "B", "C")) for (m in c("CW", "CW_CI")) cat(sprintf("%s %-6s coverage %s | width/oracle %s\n", s, m,
  paste(sprintf("%.3f", g(s, m, "coverage")), collapse = " "), paste(sprintf("%.2f", g(s, m, "width_rel")), collapse = " ")))
W <- reshape(S_all[S_all$method %in% c("CW", "CW_CI", "CMC", "Lei_naive"), c("setting", "rho", "method", "width")],
             idvar = c("setting", "rho"), timevar = "method", direction = "wide")
W1 <- W[W$rho == 1, ]; W1 <- W1[order(W1$setting), ]
cat(sprintf("CW+CI: minimum coverage over settings and rho %.3f; width relative to CW, mean over settings and rho %.3f\n",
            min(S_all$coverage[S_all$method == "CW_CI"]), mean(W$width.CW_CI / W$width.CW)))
cat(sprintf("rho = 1 (A/B/C): naive / CW+CI %s; CMC / CW+CI %s; coverage of CW %s\n",
            paste(sprintf("%.2f", W1$width.Lei_naive / W1$width.CW_CI), collapse = "/"),
            paste(sprintf("%.2f", W1$width.CMC / W1$width.CW_CI), collapse = "/"),
            paste(sprintf("%.3f", sapply(c("A", "B", "C"), function(s) tail(g(s, "CW", "coverage"), 1))), collapse = "/")))
kh <- unique(res[, c("setting", "rep", "Khat", "kappa")]); kh$true <- pmax(kh$kappa, 1 / kh$kappa)
for (s in c("A", "B", "C")) { d <- kh[kh$setting == s, ]
  cat(sprintf("K_hat %s: median %.2f [%.2f, %.2f]; K_hat >= true ratio in %d of %d replications\n", s, median(d$Khat),
              min(d$Khat), max(d$Khat), sum(d$Khat >= d$true - 1e-9), nrow(d))) }
for (m in c("CW_CI", "Lei_inexact")) cat(sprintf("setting (C), propensity quintiles, rho in {-0.5, 0, 0.5}: %-11s %.2f-%.2f\n",
  m, min(P$cov[P$method == m]), max(P$cov[P$method == m])))
best <- do.call(rbind, lapply(split(LS, LS$rho), function(d) d[which.min(d$loss), ]))
cat("smallest coverage-width loss:", paste(sprintf("rho=%s %s", best$rho, best$method), collapse = "; "), "\n")
if (RUN_ADDITIVE) {
  rmse <- function(label) mean(sapply(seq_len(N_REPS), function(rep) { B <- readRDS(cache_file(label, rep)); sqrt(mean((B$o$cate - B$test$tau)^2)) }))
  isr <- function(d) { S <- aggregate(interval_score ~ method + rho, d, mean); O <- S[S$method == "Oracle", c("rho", "interval_score")]
    S <- merge(S, O, by = "rho", suffixes = c("", "_or")); S <- S[S$method == "CW_CI", ]; mean(S$interval_score / S$interval_score_or) }
  ra <- load_res("Badd"); rb <- res[res$setting == "B", ]
  cat(sprintf("additive check (additive | tensor product): RMSE of tau_hat %.3f | %.3f; coverage of CW at rho = 1 %.3f | %.3f; interval score of CW+CI / oracle, mean over rho %.3f | %.3f\n",
              rmse("Badd"), rmse("B"), mean(ra$coverage[ra$method == "CW" & ra$rho == 1]), mean(rb$coverage[rb$method == "CW" & rb$rho == 1]), isr(ra), isr(rb)))
}
sink()
