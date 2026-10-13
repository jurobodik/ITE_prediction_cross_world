## Real-covariate case study
## =====================================================================================

suppressPackageStartupMessages({
  library(splines)
  library(quantreg)
  library(ranger)
})

ALPHA <- 0.10        # target miscoverage of the ITE prediction intervals
BETA_CI <- 0.10      # CI level 1 - beta for the CW^{+CI} term (each arm one-sided at 1 - beta/2)
DELTA <- 0.05        # rho^boot_{L,delta}: 95% lower confidence bound for rho_L(x)
N_TEST <- 1002       # new children; the remaining 3,800 form the trial
FRAC_K <- 0.08       # bandwidth = distance to the k-th nearest neighbour, k = FRAC_K * n ...
HMIN_SD <- 0.25      # ... but at least 0.25 * sd(X)
W_TOL <- 1e-4        # units with relative kernel weight below W_TOL are dropped (speed only)
N_CI_BOOT <- 100     # bootstrap draws for the CI of the CATE center
N_RHO_BOOT <- 50     # bootstrap draws (with refitting) for rho^boot_{L,delta}
DELTA_G <- 0.01      # calibrated variant: 99% Gaussian lower confidence bound for rho_L ...
BETA_G <- 0.01       # ... a 99% CI for the CATE term ...
LEVEL_G <- round((1 - ALPHA) / ((1 - BETA_G) * (1 - DELTA_G)), 6)  # ... marginal level, (1-a)(1-b)(1-d) = 0.9
NB_MAX <- 5          # number of covariate bins for conditional coverage

CPP_NAMES <- c(
  "mom_age", "mar_status", "mom_cigs_per_day", "mom_years_smoked", "mom_height",
  "mom_weight_prior", "mom_num_cardio_cond", "mom_num_pulm_cond", "mom_num_hema_cond",
  "mom_num_endocrine_cond", "mom_num_veneral_cond", "mom_num_urin_cond", "mom_num_gyne_cond",
  "mom_num_neur_cond", "mom_num_obst_compl", "mom_num_infect_dis", "mom_work_status",
  "mom_years_educ", "family_income", "housing_density", "mom_birth_place", "consanguinity",
  "socio_eco", "mom_race", "age_menarche", "dias_blood_pres", "mom_weight_birth", "dad_age",
  "dad_years_educ", "num_premes", "num_abortions", "num_prior_pregs", "num_stillbirths",
  "bayley_mental", "bayley_motor", "placental_weight", "cord_length", "sex", "apgar_1m_total",
  "apgar_5m_total", "bottle_feed_days", "breast_feed_days", "child_bilirubin",
  "child_hematocrit", "child_hemoglobin", "child_num_neur_abn", "child_num_cns_cond",
  "child_num_muscoskel", "child_num_resp_abn", "child_num_cardio_abn", "child_num_liver_abn",
  "child_num_hemo_cond", "child_num_infect", "child_num_synd", "child_num_endo_dis",
  "child_num_proc", "head_size_1yr", "gest_delivery")
CAT_COLS <- c("mar_status", "mom_birth_place", "mom_race")
MAIN_X <- "gest_delivery"
ROBUST_X <- c("gest_delivery", "mom_age", "mom_weight_prior", "bayley_mental",
              "placental_weight", "child_hematocrit")
CAUSALLIB_URL <- paste0("https://raw.githubusercontent.com/BiomedSciAI/causallib/master/",
                        "causallib/datasets/data/acic_challenge_2016")

## ------------------------------------------------------------------------------ data
load_acic <- function(data_dir = "acic_data") {
  get <- function(f) {
    p <- file.path(data_dir, f)
    if (!file.exists(p)) {
      dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
      utils::download.file(paste0(CAUSALLIB_URL, "/", f), p, quiet = TRUE, mode = "wb")
    }
    utils::read.csv(p, stringsAsFactors = FALSE)
  }
  X <- get("x.csv")
  stopifnot(ncol(X) == 58)
  names(X) <- CPP_NAMES
  for (cc in CAT_COLS) X[[cc]] <- factor(X[[cc]])
  dgps <- lapply(1:10, function(k) {
    z <- get(sprintf("zymu_%d.csv", k))
    res <- c(z$y0 - z$mu0, z$y1 - z$mu1)
    list(id = k, mu0 = z$mu0, mu1 = z$mu1, sigma = sd_pop(res))
  })
  list(X = X, dgps = dgps)
}

## ------------------------------------------------------------------------------ helpers
sd_pop <- function(v) sqrt(mean((v - mean(v))^2))
var_pop <- function(v) mean((v - mean(v))^2)

bandwidth <- function(x, x0) {
  n <- length(x)
  k <- max(floor(FRAC_K * n), 10)
  d <- sort(abs(x - x0))
  max(HMIN_SD * sd_pop(x), d[min(k, n - 1) + 1], 1e-8)
}
kern <- function(x, x0, h) exp(-0.5 * ((x - x0) / h)^2) + 1e-12
kern_w <- function(x, x0) kern(x, x0, bandwidth(x, x0))

wcov <- function(a, b, w) {
  w <- w / sum(w)
  ma <- sum(w * a)
  mb <- sum(w * b)
  sum(w * (a - ma) * (b - mb))
}
wcor <- function(a, b, w) wcov(a, b, w) / sqrt(max(wcov(a, a, w) * wcov(b, b, w), 1e-12))

make_grid <- function(x, n_max = 20) {
  q <- stats::quantile(x, c(0.01, 0.99), names = FALSE)
  xs <- x[x >= q[1] & x <= q[2]]
  u <- sort(unique(xs))
  if (length(u) <= n_max) return(as.numeric(u))
  sort(unique(stats::quantile(xs, seq(0, 1, length.out = n_max), names = FALSE)))
}

make_bins <- function(xvar, x_pop, x_new) {
  if (xvar == "gest_delivery") {          # clinical categories
    edges <- c(36, 38, 40, 41)
    labels <- c("<=36", "37-38", "39-40", "41", ">=42")
  } else {
    edges <- sort(unique(stats::quantile(x_pop, c(0.2, 0.4, 0.6, 0.8), names = FALSE)))
    labels <- paste0("Q", seq_len(length(edges) + 1))
  }
  list(bins = findInterval(x_new, edges, left.open = TRUE), labels = labels)
}

## cubic B-spline basis with knots at sample quantiles; constant extrapolation (clamping)
spline_basis <- function(x, n_knots = 6) {
  q <- sort(unique(stats::quantile(x, seq(0, 1, length.out = n_knots), names = FALSE)))
  if (length(q) < 3) q <- seq(min(x), max(x), length.out = 4)
  lo <- q[1]
  hi <- q[length(q)]
  inner <- q[-c(1, length(q))]
  function(v) splines::bs(pmin(pmax(v, lo), hi), knots = inner, degree = 3,
                          Boundary.knots = c(lo, hi), intercept = TRUE)
}

## exact linear quantile regression on a given basis (no extra intercept)
qr_coef <- function(B, y, tau) {
  fit <- tryCatch(suppressWarnings(quantreg::rq.fit(B, y, tau = tau, method = "br")),
                  error = function(e) suppressWarnings(quantreg::rq.fit(B, y, tau = tau, method = "fn")))
  fit$coefficients
}

## correlation-adjusted Euclidean distance D_rho(a, b)
Dfun <- function(a, b, rho) sqrt(pmax(a * a + b * b - 2 * rho * a * b, 0))

## c*(x) of Theorem `theorem_CI`, eq. (c_star), for one interval end:
## c* = {D(v_a, v_b) + a - v_a + b - v_b - D(a, b)} / r,  v = max(width - r/2, 0);  c* = 0 if r = 0.
## upper end: (a, b, r) = (l_0, u_1, r_u);  lower end: (l_1, u_0, r_l)
c_star <- function(a, b, r, rho) {
  va <- pmax(a - r / 2, 0)
  vb <- pmax(b - r / 2, 0)
  num <- Dfun(va, vb, rho) + a - va + b - vb - Dfun(a, b, rho)
  cc <- ifelse(r > 0, num / ifelse(r > 0, r, 1), 0)
  pmin(pmax(cc, 0), 1)
}

## ------------------------------------------------------------------------------ nuisance learners
gbm_fit <- function(d, n.trees, train.fraction) {
  suppressWarnings(gbm::gbm(.y ~ ., data = d, distribution = "gaussian", n.trees = n.trees,
                            interaction.depth = 4, shrinkage = 0.05, n.minobsinnode = 20,
                            bag.fraction = 0.8, train.fraction = train.fraction, keep.data = FALSE,
                            verbose = FALSE))
}

## gbm: number of trees selected on a validation split (unless supplied), then the model is
## refitted on all rows with that number held fixed (the same refit is used in every bootstrap)
fit_learner <- function(df, y, learner, seed, n_trees = NULL) {
  set.seed(seed)
  d <- df
  d$.y <- y
  if (learner == "gbm") {
    if (is.null(n_trees)) {
      sel <- gbm_fit(d[sample.int(nrow(d)), , drop = FALSE], 600, 0.85)
      n_trees <- max(1, suppressWarnings(gbm::gbm.perf(sel, method = "test", plot.it = FALSE)))
    }
    list(type = "gbm", model = gbm_fit(d, n_trees, 1), n.trees = n_trees)
  } else {
    list(type = "ranger", n.trees = NA,
         model = ranger::ranger(dependent.variable.name = ".y", data = d, num.trees = 500,
                                num.threads = 1, seed = seed))
  }
}

predict_learner <- function(f, newdf) {
  if (f$type == "gbm") {
    suppressWarnings(stats::predict(f$model, newdata = newdf, n.trees = f$n.trees))
  } else {
    stats::predict(f$model, data = newdf, num.threads = 1)$predictions
  }
}

## arm-specific regressions mu_t(x, z) on ALL 58 covariates (trial data only)
fit_nuisance <- function(Wtr, T, Y, learner, seed) {
  i0 <- which(T == 0)
  i1 <- which(T == 1)
  m0 <- fit_learner(Wtr[i0, , drop = FALSE], Y[i0], learner, seed + 1)
  m1 <- fit_learner(Wtr[i1, , drop = FALSE], Y[i1], learner, seed + 2)
  ## honest (out-of-fold) own-arm residuals, used only for the assumption-free bound
  set.seed(seed + 3)
  res2 <- numeric(length(Y))
  for (idx in list(i0, i1)) {
    p <- idx[sample.int(length(idx))]
    h <- floor(length(p) / 2)
    a <- p[seq_len(h)]
    b <- p[(h + 1):length(p)]
    ma <- fit_learner(Wtr[a, , drop = FALSE], Y[a], learner, seed + 4)
    mb <- fit_learner(Wtr[b, , drop = FALSE], Y[b], learner, seed + 5)
    res2[b] <- (Y[b] - predict_learner(ma, Wtr[b, , drop = FALSE]))^2
    res2[a] <- (Y[a] - predict_learner(mb, Wtr[a, , drop = FALSE]))^2
  }
  ## imputation pseudo-pairs: the counterfactual model was fit on the other arm
  Yt0 <- Y
  Yt1 <- Y
  Yt0[i1] <- predict_learner(m0, Wtr[i1, , drop = FALSE])
  Yt1[i0] <- predict_learner(m1, Wtr[i0, , drop = FALSE])
  list(m0 = m0, m1 = m1, res2 = res2, Yt0 = Yt0, Yt1 = Yt1)
}

## ------------------------------------------------------------------------------ rho estimation
## untruncated rho_L(x) = A(x) / {sigma_0(x) sigma_1(x)} on the grid.  Wu, Tu, Yu: distinct units,
## mult: their multiplicities (bootstrap counts), x_bw: the full X sample used for the bandwidth.
rhoL_ratio <- function(m0, m1, Wu, xv, Tu, Yu, grid, mult, x_bw) {
  x <- Wu[[xv]]
  i0 <- Tu == 0
  i1 <- Tu == 1
  parts <- keeps <- ws <- vector("list", length(grid))
  for (g in seq_along(grid)) {
    w <- kern(x, grid[g], bandwidth(x_bw, grid[g])) * mult
    keep <- which(w >= W_TOL * max(w))
    nd <- Wu[keep, , drop = FALSE]
    nd[[xv]] <- grid[g]
    parts[[g]] <- nd
    keeps[[g]] <- keep
    ws[[g]] <- w
  }
  nd <- do.call(rbind, parts)
  P0 <- predict_learner(m0, nd)
  P1 <- predict_learner(m1, nd)
  out <- numeric(length(grid))
  pos <- 0
  for (g in seq_along(grid)) {
    k <- keeps[[g]]
    sl <- pos + seq_along(k)
    pos <- pos + length(k)
    w <- ws[[g]]
    A <- wcov(P0[sl], P1[sl], w[k])
    s0 <- wcov(Yu[i0], Yu[i0], w[i0])
    s1 <- wcov(Yu[i1], Yu[i1], w[i1])
    out[g] <- A / sqrt(max(s0 * s1, 1e-12))
  }
  list(ratio = out, w = ws)
}

## rho_L (untruncated and clipped), assumption-free rho_tilde_L (B_L^+ = sqrt(E[r0^2|x]E[r1^2|x]),
## a Cauchy-Schwarz upper bound of B_L) and the imputation estimator rho_imp on the grid
rho_on_grid <- function(nuis, Wtr, xv, T, Y, grid) {
  n <- length(Y)
  i0 <- T == 0
  i1 <- T == 1
  rr <- rhoL_ratio(nuis$m0, nuis$m1, Wtr, xv, T, Y, grid, rep(1, n), Wtr[[xv]])
  tl <- imp <- numeric(length(grid))
  for (g in seq_along(grid)) {
    w <- rr$w[[g]]
    s0 <- wcov(Y[i0], Y[i0], w[i0])
    s1 <- wcov(Y[i1], Y[i1], w[i1])
    v0 <- sum(w[i0] * nuis$res2[i0]) / sum(w[i0])
    v1 <- sum(w[i1] * nuis$res2[i1]) / sum(w[i1])
    tl[g] <- rr$ratio[g] - sqrt(v0 * v1) / sqrt(max(s0 * s1, 1e-12))
    imp[g] <- wcor(nuis$Yt0, nuis$Yt1, w)
  }
  cl <- function(v) pmin(pmax(v, -1), 1)
  list(rhoL = cl(rr$ratio), rhoL_raw = rr$ratio, rhoTL = cl(tl), rhoImp = cl(imp))
}

## centered bootstrap lower confidence bound (Appendix appendix_rho_confidence):
## resample trial units, refit both outcome models (tuning fixed), recompute the untruncated ratio,
## bound = clip(rho_L - q*_{1-delta}(rho*_L - rho_L))
rho_boot_lcb <- function(nuis, Wtr, xv, T, Y, grid, raw, learner, B = N_RHO_BOOT, delta = DELTA) {
  n <- length(Y)
  reps <- matrix(NA_real_, B, length(grid))
  for (b in seq_len(B)) {
    s <- sample.int(n, n, replace = TRUE)
    cnt <- tabulate(s, nbins = n)
    u <- which(cnt > 0)
    Ws <- Wtr[s, , drop = FALSE]
    Ts <- T[s]
    Ys <- Y[s]
    sb <- sample.int(.Machine$integer.max - 1, 1)
    m0 <- fit_learner(Ws[Ts == 0, , drop = FALSE], Ys[Ts == 0], learner, sb, n_trees = nuis$m0$n.trees)
    m1 <- fit_learner(Ws[Ts == 1, , drop = FALSE], Ys[Ts == 1], learner, sb + 1, n_trees = nuis$m1$n.trees)
    reps[b, ] <- rhoL_ratio(m0, m1, Wtr[u, , drop = FALSE], xv, T[u], Y[u], grid, cnt[u], Wtr[[xv]][s])$ratio
  }
  q <- apply(sweep(reps, 2, raw), 2, stats::quantile, probs = 1 - delta, names = FALSE)
  list(bound = pmin(pmax(raw - q, -1), 1), reps = reps)
}

## true rho(x) and true rho_L(x) in the population (known surfaces + noise law), same kernel
oracle_rho <- function(x_pop, mu0, mu1, sigma, rho_eps, grid) {
  r <- rl <- numeric(length(grid))
  for (g in seq_along(grid)) {
    w <- kern_w(x_pop, grid[g])
    A <- wcov(mu0, mu1, w)
    v0 <- wcov(mu0, mu0, w) + sigma^2
    v1 <- wcov(mu1, mu1, w) + sigma^2
    r[g] <- (A + rho_eps * sigma^2) / sqrt(v0 * v1)
    rl[g] <- A / sqrt(v0 * v1)
  }
  list(rho = pmin(pmax(r, -1), 1), rhoL = pmin(pmax(rl, -1), 1))
}

## reference only: equal-tailed 90% quantiles of the TRUE ITE distribution given X = x0
oracle_ite_interval <- function(x_pop, mu0, mu1, sigma, rho_eps, grid, seed, M = 20) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))  # keep the main stream
  set.seed(seed)
  tau <- mu1 - mu0
  s <- rep(tau, times = M) + sigma * sqrt(2 - 2 * rho_eps) * stats::rnorm(M * length(tau))
  o <- order(s)
  ss <- s[o]
  lo <- hi <- numeric(length(grid))
  for (g in seq_along(grid)) {
    w <- rep(kern_w(x_pop, grid[g]), times = M)[o]
    cw_ <- cumsum(w) / sum(w)
    lo[g] <- ss[min(findInterval(ALPHA / 2, cw_) + 1, length(ss))]
    hi[g] <- ss[min(findInterval(1 - ALPHA / 2, cw_) + 1, length(ss))]
  }
  list(lo = lo, hi = hi)
}

## ------------------------------------------------------------------------------ intervals (X only)
## least squares on the spline basis: the conditional-MEAN fit used as the CW center
ols_coef <- function(B, y) drop(solve(crossprod(B) + diag(1e-8, ncol(B)), crossprod(B, y)))

## split conformal C_t(x) = [m(x) - l(x), m(x) + u(x)] for Y(t) given X.  center = "mean": m is the spline
## least-squares estimate of E[Y(t) | X = x] (as in the theory); center = "median" (Lei-Candes nested only).
## Widths from spline quantile regression relative to m, calibrated with max{(m - y)/l0, (y - m)/u0} and
## floored at 5% of the quantile band.  n_boot > 0: the SAME center fit is bootstrapped -> m_star
marginal_cqr <- function(x, y, x_new, levels = c(0.90, 0.95), trainprop = 0.8, n_boot = 0, center = "mean") {
  n <- length(y)
  p <- sample.int(n)
  ntr <- floor(trainprop * n)
  tr <- p[seq_len(ntr)]
  ca <- p[(ntr + 1):n]
  sb <- spline_basis(x[tr])
  Btr <- sb(x[tr])
  Bca <- sb(x[ca])
  Bn <- sb(x_new)
  qs <- sort(unique(round(c((1 - levels) / 2, (1 + levels) / 2), 6)))
  co <- lapply(qs, function(q) qr_coef(Btr, y[tr], q))
  names(co) <- as.character(qs)
  pr <- function(B, q) drop(B %*% co[[as.character(round(q, 6))]])
  fit_c <- if (center == "mean") ols_coef else function(B, yy) qr_coef(B, yy, 0.5)
  bc <- fit_c(Btr, y[tr])
  c_c <- drop(Bca %*% bc)
  c_n <- drop(Bn %*% bc)
  widths <- function(cc, lo, hi) {
    fl <- pmax(0.05 * abs(hi - lo), 1e-6)
    list(l = pmax(cc - lo, fl), u = pmax(hi - cc, fl))
  }
  out <- list(m = c_n)
  for (l in levels) {
    a <- round((1 - l) / 2, 6)
    b <- round((1 + l) / 2, 6)
    wc <- widths(c_c, pr(Bca, a), pr(Bca, b))
    sc <- pmax((c_c - y[ca]) / wc$l, (y[ca] - c_c) / wc$u)
    k <- ceiling((length(ca) + 1) * l)
    gam <- max(sort(sc)[min(k, length(sc))], 0)
    wn <- widths(c_n, pr(Bn, a), pr(Bn, b))
    out[[as.character(l)]] <- list(l = gam * wn$l, u = gam * wn$u)
  }
  if (n_boot > 0) {
    ytr <- y[tr]
    out$m_star <- t(vapply(seq_len(n_boot), function(bb) {
      s <- sample.int(ntr, ntr, replace = TRUE)
      drop(Bn %*% fit_c(Btr[s, , drop = FALSE], ytr[s]))
    }, numeric(length(x_new))))
  }
  out
}

## half-widths (r_l, r_u) of CI(x) = [tau - r_l, tau + r_u] satisfying (conditional_beta):
## one-sided centered-bootstrap bounds at level 1 - beta/2 for each arm's center
ci_from_boot <- function(c0, c1, beta = BETA_CI) {
  d0 <- sweep(c0$m_star, 2, c0$m)
  d1 <- sweep(c1$m_star, 2, c1$m)
  q0 <- apply(d0, 2, stats::quantile, probs = c(beta / 2, 1 - beta / 2), names = FALSE)
  q1 <- apply(d1, 2, stats::quantile, probs = c(beta / 2, 1 - beta / 2), names = FALSE)
  list(rl = 2 * pmax(pmax(q1[2, ], -q0[1, ]), 0),     # mu_1 >= m_1 - r_l/2, mu_0 <= m_0 + r_l/2
       ru = 2 * pmax(pmax(-q1[1, ], q0[2, ]), 0))     # mu_1 <= m_1 + r_u/2, mu_0 >= m_0 - r_u/2
}

## CW(rho) of Definition 2: [tau - D(l1,u0), tau + D(l0,u1)];  CW^{+CI} adds c* r_l and c* r_u
## (Theorem theorem_CI, one constant per end).  swapped = TRUE: width assignment of the old repo code
## normal-approximation version of (conditional_beta) for small beta: r = 2 z_{1-beta/2} max_t se_t(x),
## se_t = bootstrap standard error of the same median fit (percentiles at beta/2 need far more draws)
ci_gauss <- function(c0, c1, beta) {
  r <- 2 * stats::qnorm(1 - beta / 2) * pmax(apply(c0$m_star, 2, stats::sd), apply(c1$m_star, 2, stats::sd))
  list(rl = r, ru = r)
}

cw <- function(tau, c0, c1, rho, ci = NULL, swapped = FALSE) {
  if (swapped) {
    al <- c0$l; bl <- c1$u; au <- c0$u; bu <- c1$l
  } else {
    al <- c1$l; bl <- c0$u; au <- c0$l; bu <- c1$u
  }
  L <- tau - Dfun(al, bl, rho)
  U <- tau + Dfun(au, bu, rho)
  if (!is.null(ci)) {
    L <- L - c_star(al, bl, ci$rl, rho) * ci$rl
    U <- U + c_star(au, bu, ci$ru, rho) * ci$ru
  }
  list(L = L, U = U)
}

## Lei & Candes (2021) nested approach, inexact (cfcausal defaults: cfprop = 0.5, counterfactual
## CQR at level 1 - alpha with trainprop = 0.75, loquantile = 0.4, upquantile = 0.6)
lei_nested_inexact <- function(x, T, y, x_new, cfprop = 0.5, loq = 0.4, upq = 0.6) {
  i1 <- which(T == 1)
  i0 <- which(T == 0)
  cf1 <- i1[sample.int(length(i1), ceiling(length(i1) * cfprop))]
  cf0 <- i0[sample.int(length(i0), ceiling(length(i0) * cfprop))]
  rest <- setdiff(seq_along(y), c(cf1, cf0))
  t1 <- rest[T[rest] == 1]
  t0 <- rest[T[rest] == 0]
  lv <- 1 - ALPHA
  k <- as.character(lv)
  c1 <- marginal_cqr(x[t1], y[t1], x[cf0], levels = lv, trainprop = 0.75, center = "median")
  c0 <- marginal_cqr(x[t0], y[t0], x[cf1], levels = lv, trainprop = 0.75, center = "median")
  lo1 <- c1$m - c1[[k]]$l
  hi1 <- c1$m + c1[[k]]$u
  lo0 <- c0$m - c0[[k]]$l
  hi0 <- c0$m + c0[[k]]$u
  Lcf <- c(y[cf1] - hi0, lo1 - y[cf0])
  Ucf <- c(y[cf1] - lo0, hi1 - y[cf0])
  xcf <- c(x[cf1], x[cf0])
  sb <- spline_basis(xcf)
  Bcf <- sb(xcf)
  Bn <- sb(x_new)
  list(L = drop(Bn %*% qr_coef(Bcf, Lcf, loq)), U = drop(Bn %*% qr_coef(Bcf, Ucf, upq)))
}

rf1d <- function(x, y) {
  ranger::ranger(y ~ x, data = data.frame(x = x, y = y), num.trees = 100, min.node.size = 1,
                 num.threads = 1, seed = sample.int(.Machine$integer.max, 1))
}
pr1d <- function(m, x) stats::predict(m, data = data.frame(x = x), num.threads = 1)$predictions

## Jonkers et al. (2024) pseudo-MC CMC T-learner, non-normalised, RF(100) base learners
cmc <- function(x, T, y, x_new, mc = 100) {
  n <- length(y)
  p <- sample.int(n)
  ncal <- round(0.25 * n)
  cal <- p[seq_len(ncal)]
  nu <- p[(ncal + 1):n]
  nnc <- round(0.25 * length(nu))
  nu_cal <- nu[seq_len(nnc)]
  nu_tr <- nu[(nnc + 1):length(nu)]
  mods <- vector("list", 2)
  res <- vector("list", 2)
  for (t in 0:1) {
    a <- nu_tr[T[nu_tr] == t]
    b <- nu_cal[T[nu_cal] == t]
    mods[[t + 1]] <- rf1d(x[a], y[a])
    res[[t + 1]] <- y[b] - pr1d(mods[[t + 1]], x[b])
  }
  mu0 <- pr1d(mods[[1]], x[cal])
  mu1 <- pr1d(mods[[2]], x[cal])
  yc <- y[cal]
  nc <- length(cal)
  y1cf <- pmin(pmax(mu1 + matrix(sample(res[[2]], nc * mc, replace = TRUE), nc, mc), min(y)), max(y))
  y0cf <- pmin(pmax(mu0 + matrix(sample(res[[1]], nc * mc, replace = TRUE), nc, mc), min(y)), max(y))
  Tm <- matrix(T[cal], nc, mc)
  ite <- ifelse(Tm == 0, y1cf - yc, yc - y0cf)
  r <- sort(as.vector(ite - (mu1 - mu0)))
  m <- length(r)
  lo <- r[max(floor((m + 1) * ALPHA / 2), 1)]
  hi <- r[min(ceiling((m + 1) * (1 - ALPHA / 2)), m)]
  tau <- pr1d(mods[[2]], x_new) - pr1d(mods[[1]], x_new)
  list(L = tau + lo, U = tau + hi)
}

## Alaa et al. (2023) conformal DR meta-learner (RF, known propensity), held-out 25% calibration
dr_conformal <- function(x, T, y, x_new, nfold = 5, pi = 0.5) {
  n <- length(y)
  p <- sample.int(n)
  ncal <- round(0.25 * n)
  cal <- p[seq_len(ncal)]
  fit <- p[(ncal + 1):n]
  fit <- fit[sample.int(length(fit))]
  folds <- split(fit, cut(seq_along(fit), nfold, labels = FALSE))
  pseudo <- function(m, xs, Ts, ys) {
    m0 <- pr1d(m[[1]], xs)
    m1 <- pr1d(m[[2]], xs)
    (m1 - m0) + (Ts - pi) / (pi * (1 - pi)) * (ys - ifelse(Ts == 1, m1, m0))
  }
  models <- lapply(seq_len(nfold), function(k) {
    te <- folds[[k]]
    tr <- unlist(folds[-k], use.names = FALSE)
    m <- lapply(0:1, function(t) {
      a <- tr[T[tr] == t]
      rf1d(x[a], y[a])
    })
    list(m = m, tm = rf1d(x[te], pseudo(m, x[te], T[te], y[te])))
  })
  res <- rowMeans(sapply(models, function(o) abs(pseudo(o$m, x[cal], T[cal], y[cal]) - pr1d(o$tm, x[cal]))))
  off <- sort(res)[min(ceiling((1 - ALPHA) * (ncal + 1)), ncal)]
  tau <- rowMeans(sapply(models, function(o) pr1d(o$tm, x_new)))
  list(L = tau - off, U = tau + off)
}

## ------------------------------------------------------------------------------ evaluation
evaluate <- function(L, U, ite, bins, nb) {
  covd <- ite >= L & ite <= U
  width <- U - L
  is <- width + (2 / ALPHA) * (pmax(L - ite, 0) + pmax(ite - U, 0))
  out <- list(coverage = mean(covd), width = mean(width), median_width = stats::median(width),
              iscore = mean(is))
  bc <- rep(NA_real_, NB_MAX)
  for (b in 0:(NB_MAX - 1)) {
    mask <- bins == b
    v <- if (b < nb && sum(mask) >= 15) mean(covd[mask]) else NA_real_
    out[[paste0("cov_bin", b)]] <- v
    out[[paste0("n_bin", b)]] <- sum(mask)
    bc[b + 1] <- v
  }
  out$worst_bin_cov <- if (all(is.na(bc))) NA_real_ else min(bc, na.rm = TRUE)
  out
}

## ------------------------------------------------------------------------------ one replication
one_replication <- function(W, dgp, xvars, noise, rho_eps, seed, full, learner, cache, lcb_B = 0) {
  set.seed(seed)
  N <- nrow(W)
  sigma <- if (noise == "acic") dgp$sigma else sd_pop(dgp$mu0)
  e0 <- stats::rnorm(N)
  e1 <- rho_eps * e0 + sqrt(1 - rho_eps^2) * stats::rnorm(N)
  Y0 <- dgp$mu0 + sigma * e0
  Y1 <- dgp$mu1 + sigma * e1
  perm <- sample.int(N)
  te <- perm[seq_len(N_TEST)]
  tr <- perm[(N_TEST + 1):N]
  T <- stats::rbinom(length(tr), 1, 0.5)
  Y <- ifelse(T == 1, Y1[tr], Y0[tr])
  ite <- Y1[te] - Y0[te]
  Wtr <- W[tr, , drop = FALSE]
  Wte <- W[te, , drop = FALSE]
  
  t0 <- proc.time()[["elapsed"]]
  nuis <- fit_nuisance(Wtr, T, Y, learner, seed)
  t_nuis <- proc.time()[["elapsed"]] - t0
  r2 <- c(1 - mean((predict_learner(nuis$m0, Wte) - dgp$mu0[te])^2) / var_pop(dgp$mu0[te]),
          1 - mean((predict_learner(nuis$m1, Wte) - dgp$mu1[te])^2) / var_pop(dgp$mu1[te]))
  
  rows <- list()
  curves <- list()
  for (xv in xvars) {
    x_pop <- W[[xv]]
    x_tr <- Wtr[[xv]]
    x_te <- Wte[[xv]]
    grid <- make_grid(x_pop)
    est <- rho_on_grid(nuis, Wtr, xv, T, Y, grid)
    t_lcb <- 0
    if (lcb_B > 0 && xv == xvars[1]) {
      t1_ <- proc.time()[["elapsed"]]
      bl <- rho_boot_lcb(nuis, Wtr, xv, T, Y, grid, est$rhoL_raw, learner, B = lcb_B)
      est$rhoL_boot <- bl$bound
      ## Gaussian construction (Appendix appendix_rho_confidence) with the refitting-bootstrap SE
      est$rhoL_g01 <- pmin(pmax(est$rhoL_raw - stats::qnorm(1 - DELTA_G) * apply(bl$reps, 2, stats::sd), -1), 1)
      t_lcb <- proc.time()[["elapsed"]] - t1_
    }
    key <- paste(dgp$id, xv, signif(sigma, 8), rho_eps)
    if (is.null(cache[[key]])) {
      cache[[key]] <- list(or = oracle_rho(x_pop, dgp$mu0, dgp$mu1, sigma, rho_eps, grid),
                           q = oracle_ite_interval(x_pop, dgp$mu0, dgp$mu1, sigma, rho_eps, grid,
                                                   seed = dgp$id))
    }
    est$oracle <- cache[[key]]$or$rho
    est$oracle_rhoL <- cache[[key]]$or$rhoL
    q_lo <- stats::approx(grid, cache[[key]]$q$lo, xout = x_te, rule = 2)$y
    q_hi <- stats::approx(grid, cache[[key]]$q$hi, xout = x_te, rule = 2)$y
    rho_te <- lapply(est[names(est) != "rhoL_raw"], function(v) stats::approx(grid, v, xout = x_te, rule = 2)$y)
    
    ## marginal intervals, CATE center and its bootstrap CI, all from X only
    set.seed(7000000 + seed + 100000 * match(xv, CPP_NAMES))   # intervals: own stream (independent of the rho bootstrap)
    c0 <- marginal_cqr(x_tr[T == 0], Y[T == 0], x_te, levels = c(0.90, 0.95, LEVEL_G), n_boot = N_CI_BOOT)
    c1 <- marginal_cqr(x_tr[T == 1], Y[T == 1], x_te, levels = c(0.90, 0.95, LEVEL_G), n_boot = N_CI_BOOT)
    tau <- c1$m - c0$m
    a0 <- c0[["0.9"]]
    a1 <- c1[["0.9"]]
    n0 <- c0[["0.95"]]
    n1 <- c1[["0.95"]]
    ci <- ci_from_boot(c0, c1)
    
    ints <- list()
    ints[["Oracle ITE quantiles"]] <- list(L = q_lo, U = q_hi)          # unattainable reference
    ints[["Lei naive"]] <- list(L = tau - (n1$l + n0$u), U = tau + (n0$l + n1$u))
    ints[["CW+CI(0)"]] <- cw(tau, a0, a1, 0, ci)
    ints[["CW+CI(0.5)"]] <- cw(tau, a0, a1, 0.5, ci)
    ints[["CW(rhoL)"]] <- cw(tau, a0, a1, rho_te$rhoL)
    ints[["CW+CI(rhoL)"]] <- cw(tau, a0, a1, rho_te$rhoL, ci)
    if (!is.null(rho_te$rhoL_boot)) ints[["CW+CI(rhoL_boot)"]] <- cw(tau, a0, a1, rho_te$rhoL_boot, ci)
    if (!is.null(rho_te$rhoL_g01))     # calibrated variant: marginals at LEVEL_G, 99% CI, 99% bound
      ints[["CW+CI(rhoL_g01,cal)"]] <- cw(tau, c0[[as.character(LEVEL_G)]], c1[[as.character(LEVEL_G)]],
                                          rho_te$rhoL_g01, ci_gauss(c0, c1, BETA_G))
    ints[["CW+CI(rhoTL)"]] <- cw(tau, a0, a1, rho_te$rhoTL, ci)
    ints[["CW+CI(rhoImp)"]] <- cw(tau, a0, a1, rho_te$rhoImp, ci)
    ints[["CW(oracle)"]] <- cw(tau, a0, a1, rho_te$oracle)
    ints[["CW+CI(oracle)"]] <- cw(tau, a0, a1, rho_te$oracle, ci)
    ints[["CW+CI(rhoL) swapped"]] <- cw(tau, a0, a1, rho_te$rhoL, ci, swapped = TRUE)
    if (full) {
      set.seed(8000000 + seed + 100000 * match(xv, CPP_NAMES))
      ints[["Lei nested"]] <- lei_nested_inexact(x_tr, T, Y, x_te)
      ints[["CMC"]] <- cmc(x_tr, T, Y, x_te)
      ints[["DR-conformal"]] <- dr_conformal(x_tr, T, Y, x_te)
    }
    
    bb <- make_bins(xv, x_pop, x_te)
    rho_used <- list("CW+CI(0)" = 0, "CW+CI(0.5)" = 0.5, "CW(rhoL)" = rho_te$rhoL,
                     "CW+CI(rhoL)" = rho_te$rhoL, "CW+CI(rhoL_boot)" = rho_te$rhoL_boot, "CW+CI(rhoL_g01,cal)" = rho_te$rhoL_g01,
                     "CW+CI(rhoTL)" = rho_te$rhoTL, "CW+CI(rhoImp)" = rho_te$rhoImp,
                     "CW(oracle)" = rho_te$oracle, "CW+CI(oracle)" = rho_te$oracle)
    has_b <- !is.null(est$rhoL_boot)
    base <- list(dgp = dgp$id, xvar = xv, noise = noise, rho_eps = rho_eps, sigma = sigma,
                 seed = seed, r2_mu0 = r2[1], r2_mu1 = r2[2], t_nuis = t_nuis, t_lcb = t_lcb,
                 rho_oracle_mean = mean(rho_te$oracle), rhoL_true_mean = mean(rho_te$oracle_rhoL),
                 rhoL_mean = mean(rho_te$rhoL), rhoTL_mean = mean(rho_te$rhoTL),
                 rhoImp_mean = mean(rho_te$rhoImp),
                 rhoL_boot_mean = if (has_b) mean(rho_te$rhoL_boot) else NA_real_,
                 frac_rhoL_below_oracle = mean(est$rhoL <= est$oracle + 1e-9),
                 frac_rhoLboot_below_rhoL = if (has_b) mean(est$rhoL_boot <= est$oracle_rhoL + 1e-9) else NA_real_,
                 rhoL_g01_mean = if (has_b) mean(rho_te$rhoL_g01) else NA_real_,
                 frac_rhoLg01_below_rhoL = if (has_b) mean(est$rhoL_g01 <= est$oracle_rhoL + 1e-9) else NA_real_,
                 ci_mean = mean(ci$rl + ci$ru), bin_labels = paste(bb$labels, collapse = "|"))
    for (nm in names(ints)) {
      ev <- evaluate(ints[[nm]]$L, ints[[nm]]$U, ite, bb$bins, length(bb$labels))
      ru <- if (is.null(rho_used[[nm]])) NA_real_ else mean(rho_used[[nm]])
      rows[[length(rows) + 1]] <- as.data.frame(c(base, list(method = nm), ev, list(rho_used_mean = ru)),
                                                stringsAsFactors = FALSE, check.names = FALSE)
    }
    cv <- data.frame(dgp = dgp$id, xvar = xv, noise = noise, rho_eps = rho_eps, seed = seed, x0 = grid,
                     rhoL = est$rhoL, rhoL_raw = est$rhoL_raw, rhoTL = est$rhoTL, rhoImp = est$rhoImp,
                     rhoL_boot = if (has_b) est$rhoL_boot else NA_real_,
                     rhoL_g01 = if (has_b) est$rhoL_g01 else NA_real_,
                     oracle = est$oracle, oracle_rhoL = est$oracle_rhoL, stringsAsFactors = FALSE)
    curves[[length(curves) + 1]] <- cv
  }
  list(rows = do.call(rbind, rows), curves = do.call(rbind, curves))
}

## ------------------------------------------------------------------------------ driver
collect <- function(out, study) {
  d <- file.path(out, study)
  for (pp in list(c("rows", "results"), c("curves", "rho_curves"))) {
    fs <- sort(list.files(d, pattern = paste0("^", pp[1], "_.*\\.csv$"), full.names = TRUE))
    if (length(fs)) {
      dfs <- lapply(fs, utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE)
      cols <- unique(unlist(lapply(dfs, names)))
      dfs <- lapply(dfs, function(z) { for (cc in setdiff(cols, names(z))) z[[cc]] <- NA; z[cols] })
      utils::write.csv(do.call(rbind, dfs), file.path(out, sprintf("%s_%s.csv", pp[2], study)), row.names = FALSE)
    }
  }
}

parse_args <- function(args) {
  opt <- list(study = "main", reps = "15", dgps = "1,2,3,4,5,6,7,8,9,10", out = "results_R",
              data_dir = "acic_data", learner = "gbm", lcb_reps = NA, lcb_B = as.character(N_RHO_BOOT),
              collect_only = FALSE)
  i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[i])
    if (key == "collect_only") {
      opt$collect_only <- TRUE
      i <- i + 1
    } else if (grepl("=", key, fixed = TRUE)) {
      kv <- strsplit(key, "=", fixed = TRUE)[[1]]
      opt[[kv[1]]] <- kv[2]
      i <- i + 1
    } else {
      opt[[key]] <- args[i + 1]
      i <- i + 2
    }
  }
  opt$reps <- as.integer(opt$reps)
  opt$lcb_B <- as.integer(opt$lcb_B)
  opt$lcb_reps <- if (is.na(opt$lcb_reps)) (if (opt$study == "main") 3L else 0L) else as.integer(opt$lcb_reps)
  opt
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  opt <- parse_args(args)
  if (opt$collect_only) {
    collect(opt$out, opt$study)
    return(invisible(NULL))
  }
  if (opt$learner == "gbm" && !requireNamespace("gbm", quietly = TRUE))
    stop("package 'gbm' not installed: install.packages('gbm') or use --learner ranger")
  dat <- load_acic(opt$data_dir)
  W <- dat$X
  keep <- as.integer(strsplit(opt$dgps, ",", fixed = TRUE)[[1]])
  dgps <- Filter(function(d) d$id %in% keep, dat$dgps)
  cfg <- switch(opt$study,
                main = list(configs = list(c("acic", "0")), xvars = MAIN_X, full = TRUE, code = 1),
                xrobust = list(configs = list(c("acic", "0")), xvars = ROBUST_X, full = FALSE, code = 2),
                noise = list(configs = list(c("r2half", "-0.5"), c("r2half", "0"), c("r2half", "0.5")),
                             xvars = MAIN_X, full = FALSE, code = 3),
                stop("--study must be main, xrobust or noise"))
  rep_dir <- file.path(opt$out, opt$study)
  dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
  cache <- new.env()
  t_start <- proc.time()[["elapsed"]]
  for (ci in seq_along(cfg$configs)) {
    noise <- cfg$configs[[ci]][1]
    rho_eps <- as.numeric(cfg$configs[[ci]][2])
    for (d in dgps) {
      for (rep in seq_len(opt$reps) - 1) {
        seed <- 1000000 * cfg$code + 10000 * (ci - 1) + 100 * d$id + rep
        fr <- file.path(rep_dir, sprintf("rows_%d.csv", seed))
        fc <- file.path(rep_dir, sprintf("curves_%d.csv", seed))
        if (file.exists(fr) && file.exists(fc)) next
        t0 <- proc.time()[["elapsed"]]
        o <- one_replication(W, d, cfg$xvars, noise, rho_eps, seed, cfg$full, opt$learner, cache,
                             lcb_B = if (rep < opt$lcb_reps) opt$lcb_B else 0)
        utils::write.csv(o$curves, paste0(fc, ".tmp"), row.names = FALSE)
        file.rename(paste0(fc, ".tmp"), fc)
        utils::write.csv(o$rows, paste0(fr, ".tmp"), row.names = FALSE)
        file.rename(paste0(fr, ".tmp"), fr)
        m <- o$rows[o$rows$xvar == cfg$xvars[1], ]
        g <- function(meth, col) m[m$method == meth, col][1]
        cat(sprintf(paste0("[%s] noise=%s rho_eps=%+.1f dgp=%2d rep=%2d %5.1fs  rho_or=%.2f rhoL=%.2f  ",
                           "cov(CW+CI rhoL)=%.3f w=%.2f | naive w=%.2f%s  [total %6.0fs]\n"),
                    opt$study, noise, rho_eps, d$id, rep, proc.time()[["elapsed"]] - t0,
                    g("CW+CI(oracle)", "rho_oracle_mean"), g("CW+CI(rhoL)", "rhoL_mean"),
                    g("CW+CI(rhoL)", "coverage"), g("CW+CI(rhoL)", "width"), g("Lei naive", "width"),
                    if (any(m$method == "CW+CI(rhoL_boot)"))
                      sprintf(" boot=%.2f cov=%.3f", g("CW+CI(rhoL_boot)", "rhoL_boot_mean"),
                              g("CW+CI(rhoL_boot)", "coverage")) else "",
                    proc.time()[["elapsed"]] - t_start))
        flush.console()
      }
    }
  }
  collect(opt$out, opt$study)
}

if (sys.nframe() == 0L) main()
