## =====================================================================================================
##  Experiment 2 (Section 5.2 and Appendix 'Beyond elliptical distributions', label appendix_non_elliptical):
##  robustness of CW(rho) and CW^{+CI}(rho) to non-Gaussian and non-elliptical joint distributions.
##  CW^{+CI}_K(rho) in the finite-sample study E2 uses K = K_hat (Appendix 'Choice of K'), computed from the
##  calibration splits of the marginal CQR (see marginal_cqr() and K_hat() below).
##  Self-contained: needs no other file.
##
##  Outputs (in OUT_DIR, default cw_nonelliptical_output; set CW_OUT_DIR to change it):
##     figures/robust_E1_copulas_sym.pdf   population-level coverage of CW(rho), sigma1 = 2 sigma0
##     figures/robust_E1_cop_sym.pdf       CW(rho) vs copula-based intervals (Gaussian family, union)
##     figures/robust_E2_summary_sym.pdf   finite samples (baseline design, n = 1000), incl. CW^{+CI}
##     tables/table_E1cop_sym_body.tex     Table tab_E1cop body
##     tables/table_E2_sym_body.tex        Table tab_E2 body
##     text_numbers.txt                    every number quoted in the text of the section
##
##  Required packages: quantreg, splines, ggplot2, and an LP solver: 'highs' (CRAN) or 'Rglpk'.
##
##  Usage:   Rscript exp2_nonelliptical_5.2.R                    (full run, about 2 hours on one core:
##                                                               E1 ~35 min, E2 ~1 h, E3 ~2 min; resumable per stage)
##           CW_FAST=TRUE Rscript exp2_nonelliptical_5.2.R       (quick test, ~20 min, coarse grids,
##                                                               4 replications; numbers are approximate)
## =====================================================================================================

suppressPackageStartupMessages({
  library(quantreg)
  library(splines)
  library(ggplot2)
})

## ----------------------------------------------------------------------------------- settings
FAST    <- as.logical(Sys.getenv("CW_FAST", "FALSE"))   # TRUE = quick test run (coarse grids, few replications)
OUT_DIR <- Sys.getenv("CW_OUT_DIR", "cw_nonelliptical_output")
ALPHA   <- 0.10
RHOS    <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
SRAT    <- c(1, 2)               # sigma1 / sigma0 (sigma0 = 1)
N_TRAIN <- 1000; N_TEST <- 5000
R_REPS  <- if (FAST) 4 else 50   # replications per configuration (finite samples)
M_DRAWS <- if (FAST) 5000 else 20000   # copula draws for the copula-based intervals
NP      <- if (FAST) 401 else 801      # grid size for the Pearson-correlation quadrature (error ~1e-4)
NZ      <- if (FAST) 2401 else 4801    # grid size for the ITE-distribution quadrature
M_LP    <- if (FAST) 50 else 120       # patchwork grid for the worst-case linear programs
set.seed(20240601)
dir.create(file.path(OUT_DIR, "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"), showWarnings = FALSE)
cache_file <- function(name) file.path(OUT_DIR, "cache", paste0(name, if (FAST) "_fast" else "", ".rds"))
cached <- function(name, expr) {            # evaluate expr once and store the result
  f <- cache_file(name)
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time(); val <- expr; saveRDS(val, f)
  message(sprintf("[%s] computed in %.1f min", name, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  val
}

## ----------------------------------------------------------------------------------- helpers
ZC <- 37                                   # Phi(-37) ~ 1e-300
clipz  <- function(z) pmin(pmax(z, -ZC), ZC)
lndtr  <- function(z) pnorm(z, log.p = TRUE)                     # log Phi(z), accurate in both tails
score_from_logs <- function(lp, lq) clipz(ifelse(lp < log(0.5), qnorm(lp, log.p = TRUE), -qnorm(lq, log.p = TRUE)))
score_from_cdf_sf <- function(p, q) score_from_logs(log(pmax(p, 1e-320)), log(pmax(q, 1e-320)))
softplus <- function(r) ifelse(r > 0, r + log1p(exp(-r)), log1p(exp(r)))
logaddexp <- function(a, b) { m <- pmax(a, b); m + log1p(exp(-abs(a - b))) }
log_expm1 <- function(x) ifelse(x > 30, x, log(expm1(pmin(x, 30))))
Drho <- function(a, b, r) sqrt(pmax(a^2 + b^2 - 2 * r * a * b, 0))

## ----------------------------------------------------------------------------------- marginals
## Every standardised marginal law (mean 0, variance 1) is given through its normal-score quantile
## g(z) = G^{-1}(Phi(z)) and its normal score zscore(y) = Phi^{-1}(G(y)).
new_marg <- function(name, label, g, zscore) {
  list(name = name, label = label, g = g, zscore = zscore,
       cdf = function(y) pnorm(zscore(y)), sf = function(y) pnorm(-zscore(y)),
       rvs = function(n) g(rnorm(n)))
}
tq_fun <- function(nu) {                   # Phi(z) -> t_nu quantile (unscaled), via accurate log-probabilities
  function(z) ifelse(z < 0, qt(lndtr(z), nu, log.p = TRUE), qt(lndtr(-z), nu, lower.tail = FALSE, log.p = TRUE))
}
tz_fun <- function(nu) {                   # t_nu value -> normal score
  function(x) score_from_logs(pt(x, nu, log.p = TRUE), pt(x, nu, lower.tail = FALSE, log.p = TRUE))
}
m_gauss <- new_marg("gauss", "Gaussian", function(z) z, function(y) clipz(y))
m_t <- function(nu) {
  s <- sqrt(nu / (nu - 2)); tq <- tq_fun(nu); tz <- tz_fun(nu)
  new_marg(paste0("t", nu), paste0("Student-t", nu), function(z) tq(z) / s, function(y) tz(y * s))
}
m_laplace <- local({
  b <- 1 / sqrt(2)
  new_marg("laplace", "Laplace",
           function(z) ifelse(z < 0, b * (log(2) + lndtr(z)), -b * (log(2) + lndtr(-z))),
           function(y) { lt <- log(0.5) - abs(y) / b
           clipz(ifelse(y < 0, qnorm(lt, log.p = TRUE), -qnorm(lt, log.p = TRUE))) })
})
m_unif <- local({
  a <- sqrt(3)
  new_marg("unif", "Uniform",
           function(z) ifelse(z < 0, -a + 2 * a * pnorm(z), a - 2 * a * pnorm(-z)),
           function(y) score_from_cdf_sf(pmin(pmax((y + a) / (2 * a), 0), 1), pmin(pmax((a - y) / (2 * a), 0), 1)))
})
m_bimodal <- local({
  a <- 0.9; s <- sqrt(1 - a^2)
  zs <- function(y) score_from_cdf_sf(0.5 * (pnorm((y + a) / s) + pnorm((y - a) / s)),
                                      0.5 * (pnorm(-(y + a) / s) + pnorm(-(y - a) / s)))
  yg <- seq(-a - 40 * s, a + 40 * s, length.out = 400001); zg <- zs(yg)
  keep <- abs(zg) < ZC - 0.01
  new_marg("bimodal", "Bimodal", function(z) approx(zg[keep], yg[keep], z, rule = 2)$y, zs)
})
m_lognorm <- local({
  sg <- 0.5; m <- exp(sg^2 / 2); sd <- sqrt((exp(sg^2) - 1) * exp(sg^2))
  new_marg("lognorm", "Log-normal", function(z) (exp(sg * z) - m) / sd,
           function(y) { x <- y * sd + m; clipz(ifelse(x > 0, log(pmax(x, 1e-300)) / sg, -Inf)) })
})
m_exp <- local({
  k <- 1
  new_marg("exp", "Exponential",
           function(z) (ifelse(z < 0, qgamma(lndtr(z), k, log.p = TRUE),
                               qgamma(lndtr(-z), k, lower.tail = FALSE, log.p = TRUE)) - k) / sqrt(k),
           function(y) { x <- k + sqrt(k) * y
           score_from_logs(pgamma(x, k, log.p = TRUE), pgamma(x, k, lower.tail = FALSE, log.p = TRUE)) })
})
negated <- function(m) new_marg(paste0("neg_", m$name), paste("neg.", m$label),
                                function(z) -m$g(-z), function(y) -m$zscore(-y))
MARG <- list(gauss = m_gauss, t3 = m_t(3), laplace = m_laplace, unif = m_unif, bimodal = m_bimodal,
             lognorm = m_lognorm, exp = m_exp)
MARG_CONFIGS <- list(              # (Y(0) noise, Y(1) noise, label)
  "gauss/gauss"     = c("gauss", "gauss", "Gaussian"),
  "t3/t3"           = c("t3", "t3", "Student-t3"),
  "laplace/laplace" = c("laplace", "laplace", "Laplace"),
  "unif/unif"       = c("unif", "unif", "Uniform"),
  "bimodal/bimodal" = c("bimodal", "bimodal", "Bimodal"),
  "lognorm/lognorm" = c("lognorm", "lognorm", "Log-normal"),
  "exp/exp"         = c("exp", "exp", "Exponential"),
  "gauss/lognorm"   = c("gauss", "lognorm", "Gauss. -> log-normal"),
  "gauss/t3"        = c("gauss", "t3", "Gauss. -> t3"))

## ----------------------------------------------------------------------------------- copulas
## hz(zv, zu, th) = P(Z_V <= zv | Z_U = zu) on the normal-score scale (zu, zv of equal length).
## sample_z(n, th) returns normal scores (zu, zv).
bisect_inv <- function(hz, w, zu, th, iters = 45) {
  lo <- rep(-ZC, length(w)); hi <- rep(ZC, length(w))
  for (i in seq_len(iters)) { mid <- (lo + hi) / 2; up <- hz(mid, zu, th) < w
  lo <- ifelse(up, mid, lo); hi <- ifelse(up, hi, mid) }
  (lo + hi) / 2
}
cop_gauss <- list(name = "gauss", label = "Gaussian", sing = FALSE,
                  hz = function(zv, zu, th) pnorm((zv - th * zu) / sqrt(1 - th^2)),
                  sample_z = function(n, th) { zu <- rnorm(n); list(zu = zu, zv = th * zu + sqrt(1 - th^2) * rnorm(n)) })
cop_t <- local({
  nu <- 4; tz <- tz_fun(nu); zg <- seq(-ZC, ZC, by = 5e-4); tg <- tq_fun(nu)(zg)
  tq <- function(z) approx(zg, tg, z, rule = 2)$y
  hz <- function(zv, zu, th) { x <- tq(zu); y <- tq(zv)
  pt((y - th * x) / sqrt((nu + x^2) * (1 - th^2) / (nu + 1)), nu + 1) }
  list(name = "t4", label = "Student-t4", sing = FALSE, hz = hz,
       sample_z = function(n, th) { zu <- rnorm(n); x <- tq(zu)
       y <- th * x + sqrt((nu + x^2) * (1 - th^2) / (nu + 1)) * qt(runif(n), nu + 1)
       list(zu = zu, zv = tz(y)) })
})
cop_clayton <- list(name = "clayton", label = "Clayton", sing = FALSE,
                    hz = function(zv, zu, th) { lr <- log_expm1(-th * lndtr(zv)) + th * lndtr(zu)
                    exp(-(1 + 1 / th) * softplus(lr)) },
                    sample_z = function(n, th) { zu <- rnorm(n); w <- runif(n)
                    lr <- log(expm1(-log(w) * th / (1 + th))) - th * lndtr(zu)
                    list(zu = zu, zv = clipz(qnorm(-softplus(lr) / th, log.p = TRUE))) })
cop_gumbel <- list(name = "gumbel", label = "Gumbel", sing = FALSE,
                   hz = function(zv, zu, th) { a <- -lndtr(zu); b <- -lndtr(zv); la <- log(a); lb <- log(b)
                   lw <- logaddexp(th * la, th * lb) / th
                   exp(pmin(-exp(lw) + (th - 1) * (la - lw) + a, 0)) })
cop_joe <- list(name = "joe", label = "Joe", sing = FALSE,
                hz = function(zv, zu, th) { lub <- lndtr(-zu); lvb <- lndtr(-zv); lP <- th * lub; lQ <- th * lvb
                P <- exp(lP); Q <- exp(lQ)
                lS <- logaddexp(lP, lQ + log1p(-pmin(P, 1 - 1e-16)))
                exp(pmin((th - 1) * lub + log1p(-pmin(Q, 1)) + (1 / th - 1) * lS, 0)) })
frank_hpos <- function(zv, zu, th) {
  u <- pnorm(zu); v <- pnorm(zv); vb <- pnorm(-zv)
  l1 <- -th * u + log(-expm1(-th * v)); l2 <- -th * v + log(-expm1(-th * vb))
  1 / (1 + exp(pmin(pmax(l2 - l1, -700), 700)))
}
cop_frank <- list(name = "frank", label = "Frank", sing = FALSE,
                  hz = function(zv, zu, th) {
                    if (abs(th) < 1e-6) return(pnorm(zv) + 0 * zu)
                    if (th > 0) frank_hpos(zv, zu, th) else 1 - frank_hpos(-zv, zu, -th) })
cop_plackett <- list(name = "plackett", label = "Plackett", sing = FALSE,
                     hz = function(zv, zu, th) { u <- pnorm(zu); v <- pnorm(zv); S <- 1 + (th - 1) * (u + v)
                     disc <- S^2 - 4 * th * (th - 1) * u * v
                     pmin(pmax(0.5 - (S - 2 * th * v) / (2 * sqrt(pmax(disc, 1e-300))), 0), 1) })
cop_gmix <- list(name = "gmix", label = "Gauss. mixture (+-0.9)", sing = FALSE,
                 hz = function(zv, zu, th) th * pnorm((zv - 0.9 * zu) / sqrt(1 - 0.81)) + (1 - th) * pnorm((zv + 0.9 * zu) / sqrt(1 - 0.81)),
                 sample_z = function(n, th) { zu <- rnorm(n); r <- ifelse(runif(n) < th, 0.9, -0.9)
                 list(zu = zu, zv = r * zu + sqrt(1 - r^2) * rnorm(n)) })
cop_frechet <- list(name = "frechet", label = "Frechet M/W (singular)", sing = TRUE,
                    sample_z = function(n, th) { zu <- rnorm(n); list(zu = zu, zv = ifelse(runif(n) < th, zu, -zu)) })
add_sampler <- function(cop) {             # generic conditional-inversion sampler
  if (is.null(cop$sample_z)) cop$sample_z <- function(n, th) {
    zu <- rnorm(n); list(zu = zu, zv = bisect_inv(cop$hz, runif(n), zu, th)) }
  cop
}
rotate <- function(cop) {                  # U -> 1 - U: negative dependence
  base <- cop
  r <- list(name = paste0(cop$name, "_rot"), label = paste(cop$label, "(rot.)"), sing = FALSE,
            hz = function(zv, zu, th) base$hz(zv, -zu, th))
  r$sample_z <- function(n, th) { s <- base$sample_z(n, th); list(zu = -s$zu, zv = s$zv) }
  r
}
COPS <- lapply(list(gauss = cop_gauss, t4 = cop_t, frank = cop_frank, plackett = cop_plackett,
                    clayton = cop_clayton, gumbel = cop_gumbel, joe = cop_joe, gmix = cop_gmix,
                    frechet = cop_frechet), add_sampler)
ROT <- lapply(list(clayton = COPS$clayton, gumbel = COPS$gumbel, joe = COPS$joe), rotate)
REG <- c("gauss", "t4", "frank", "plackett", "clayton", "gumbel", "joe")   # the seven parametric families
INDEP <- c("frank", "plackett", "clayton", "gumbel", "joe")                # independence at rho = 0
FAM_LABEL <- sapply(COPS, function(c) c$label)
cop_for <- function(fam, rho) if (rho < 0 && fam %in% names(ROT)) ROT[[fam]] else COPS[[fam]]

## ----------------------------------------------------------------------------------- quadrature
trap_w <- function(z) { d <- z[2] - z[1]; w <- rep(d, length(z)); w[c(1, length(z))] <- d / 2; w * dnorm(z) }
GZ  <- seq(-12, 12, length.out = NZ); GW <- trap_w(GZ)       # outer grid for ITE distribution
GPo <- seq(-12, 12, length.out = NP); GPw <- trap_w(GPo)     # Pearson grid (nodes = edges)
GPm <- (GPo[-1] + GPo[-NP]) / 2

## Pearson correlation of (eps0, eps1) under copula `cop` with parameter th and marginals m0, m1
pearson <- function(cop, th, m0, m1) {
  if (cop$sing) {                          # Frechet: linear in p
    rM <- sum(GW * m0$g(GZ) * m1$g(GZ)); rW <- sum(GW * m0$g(GZ) * m1$g(-GZ))
    return(th * rM + (1 - th) * rW)
  }
  ZU <- matrix(GPo, NP, NP); ZV <- matrix(GPo, NP, NP, byrow = TRUE)
  H  <- matrix(cop$hz(as.vector(ZV), as.vector(ZU), th), NP, NP)
  dH <- H[, -1, drop = FALSE] - H[, -NP, drop = FALSE]
  g1m <- m1$g(GPm); cm <- as.vector(dH %*% g1m); cm2 <- as.vector(dH %*% g1m^2)
  g0 <- m0$g(GPo); w <- GPw
  E0 <- sum(w * g0); E1 <- sum(w * cm)
  (sum(w * g0 * cm) - E0 * E1) / sqrt((sum(w * g0^2) - E0^2) * (sum(w * cm2) - E1^2))
}

## copula parameter giving Pearson correlation rho (NA if not attainable)
theta_for_rho <- function(fam, rho, m0, m1) {
  if (rho == 0 && fam %in% INDEP) return(0)          # independence (handled by Frank with theta = 0)
  cop <- cop_for(fam, rho)
  if (fam %in% c("gmix", "frechet")) {
    a <- pearson(cop, 0, m0, m1); b <- pearson(cop, 1, m0, m1); p <- (rho - a) / (b - a)
    return(if (p >= -1e-4 && p <= 1 + 1e-4) min(max(p, 0), 1) else NA)
  }
  br <- switch(fam, gauss = c(-0.99995, 0.99995), t4 = c(-0.99995, 0.99995), frank = c(-80, 80),
               plackett = c(log(1e-5), log(1e5)), clayton = c(1e-6, 60), gumbel = c(1 + 1e-9, 40), joe = c(1 + 1e-9, 60))
  tr <- if (fam == "plackett") exp else identity
  f <- function(x) pearson(cop, tr(x), m0, m1) - rho
  fl <- f(br[1]); fh <- f(br[2])
  if (fl * fh > 0) return(NA)
  tr(uniroot(f, br, f.lower = fl, f.upper = fh, tol = 1e-9)$root)
}
true_cop <- function(fam, rho) if (rho == 0 && fam %in% INDEP) COPS$frank else cop_for(fam, rho)

## P(s1 eps1 - s0 eps0 <= t) for a vector t
ite_cdf <- function(t, cop, th, m0, m1, s0, s1) {
  g0 <- m0$g(GZ)
  if (cop$sing) {
    dM <- s1 * m1$g(GZ) - s0 * g0; dW <- s1 * m1$g(-GZ) - s0 * g0
    return(sapply(t, function(tt) th * sum(GW[dM <= tt]) + (1 - th) * sum(GW[dW <= tt])))
  }
  sapply(t, function(tt) sum(GW * cop$hz(m1$zscore((tt + s0 * g0) / s1), GZ, th)))
}
find_quantile <- function(Fun, p) {
  lo <- -2; hi <- 2
  while (Fun(lo) > p) lo <- 2 * lo
  while (Fun(hi) < p) hi <- 2 * hi
  uniroot(function(t) Fun(t) - p, c(lo, hi), tol = 1e-10)$root
}
ite_quantile <- function(p, cop, th, m0, m1, s0, s1) find_quantile(function(t) ite_cdf(t, cop, th, m0, m1, s0, s1), p)

## oracle intervals (relative to tau)
widths <- function(m, s, level) { a <- (1 - level) / 2; c(l = -s * m$g(qnorm(a)), u = s * m$g(qnorm(1 - a))) }
cw_oracle <- function(rho, m0, m1, s0, s1, alpha = ALPHA) {
  w0 <- widths(m0, s0, 1 - alpha); w1 <- widths(m1, s1, 1 - alpha)
  c(L = -Drho(w1[["l"]], w0[["u"]], rho), U = Drho(w0[["l"]], w1[["u"]], rho))
}
naive_oracle <- function(m0, m1, s0, s1, alpha = ALPHA) {
  w0 <- widths(m0, s0, 1 - alpha / 2); w1 <- widths(m1, s1, 1 - alpha / 2)
  c(L = -(w1[["l"]] + w0[["u"]]), U = w0[["l"]] + w1[["u"]])
}
analyse <- function(cop, th, rho, m0, m1, s0, s1, alpha = ALPHA) {
  b <- cw_oracle(rho, m0, m1, s0, s1, alpha); F <- unname(ite_cdf(unname(b), cop, th, m0, m1, s0, s1))
  qlo <- ite_quantile(alpha / 2, cop, th, m0, m1, s0, s1); qhi <- ite_quantile(1 - alpha / 2, cop, th, m0, m1, s0, s1)
  nv <- naive_oracle(m0, m1, s0, s1, alpha)
  c(cov = F[2] - F[1], miss_lo = F[1], miss_up = 1 - F[2], width_cw = b[["U"]] - b[["L"]],
    width_true = qhi - qlo, width_naive = nv[["U"]] - nv[["L"]], qlo = qlo, qhi = qhi)
}

## ----------------------------------------------------------------------------------- latent-factor models
## eps0 = lam H + sqrt(1-lam^2) E0,  eps1 = sgn lam H + sqrt(1-lam^2) E1,  cor(eps0, eps1) = sgn lam^2
twopoint <- list(name = "twopoint", label = "two-point", nodes = list(x = c(-1, 1), w = c(0.5, 0.5)),
                 rvs = function(n) sample(c(-1, 1), n, replace = TRUE))
LAT_LAWS <- c(MARG, list(twopoint = twopoint, negexp = negated(m_exp)))
LAT_CONFIGS <- list("H=exp,E=gauss" = c("exp", "gauss"), "H=gauss,E=t3" = c("gauss", "t3"),
                    "H=t3,E=gauss" = c("t3", "gauss"), "H=exp,E=negexp" = c("exp", "negexp"),
                    "H=twopoint,E=gauss" = c("twopoint", "gauss"), "H=gauss,E=exp" = c("gauss", "exp"))
ZQ <- seq(-9, 9, length.out = 361); WQ <- trap_w(ZQ); WQ <- WQ / sum(WQ)
law_nodes <- function(L) if (!is.null(L$nodes)) L$nodes else list(x = L$g(ZQ), w = WQ)
.tab <- new.env()
law_table <- function(E) {
  if (is.null(.tab[[E$name]])) { y <- sinh(seq(-asinh(1e7), asinh(1e7), length.out = 60001))
  assign(E$name, list(y = y, F = E$cdf(y), S = E$sf(y)), envir = .tab) }
  .tab[[E$name]]
}
tcdf <- function(E, x) { T <- law_table(E); approx(T$y, T$F, x, rule = 2)$y }
tsf  <- function(E, x) { T <- law_table(E); approx(T$y, T$S, x, rule = 2)$y }
numeric_marg <- function(a, H, b, E, ny = 20001) {   # full law of a H + b E on a grid (used in E2)
  hn <- law_nodes(H)
  ymax <- abs(a) * max(abs(hn$x)) + abs(b) * max(abs(law_nodes(E)$x)) + 1
  y <- sinh(seq(-asinh(ymax), asinh(ymax), length.out = ny))
  U <- outer(y, -a * hn$x, "+") / b
  P <- as.vector(matrix(tcdf(E, U), ny) %*% hn$w); Q <- as.vector(matrix(tsf(E, U), ny) %*% hn$w)
  zg <- cummax(score_from_cdf_sf(P, Q) + 1e-13 * seq_along(y))
  new_marg("numeric", "numeric", function(z) approx(zg, y, z, rule = 2)$y, function(yy) approx(y, zg, yy, rule = 2)$y)
}
conv_marg <- function(a, H, b, E) {           # light version: quantiles by root-finding (used in E1)
  hn <- law_nodes(H)
  cdf1 <- function(y) sum(hn$w * tcdf(E, (y - a * hn$x) / b))
  list(name = "conv", g = function(z) sapply(z, function(zz) find_quantile(cdf1, pnorm(zz))))
}
latent_marginals <- function(lam, sgn, H, E, full = FALSE) {
  if (lam == 0) return(list(E, E))
  cc <- sqrt(1 - lam^2); mk <- if (full) numeric_marg else conv_marg
  list(mk(lam, H, cc, E), mk(sgn * lam, H, cc, E))
}
latent_ite_cdf <- function(t, lam, sgn, H, E, s0, s1) {
  cc <- sqrt(1 - lam^2); hn <- law_nodes(H); en <- law_nodes(E); k <- (s1 * sgn - s0) * lam
  W <- outer(hn$w, en$w); base <- outer(-k * hn$x, cc * s0 * en$x, "+") / (cc * s1)
  sapply(t, function(tt) sum(W * tcdf(E, tt / (cc * s1) + base)))
}

## ===================================================================================== E1: population level
row_list <- function(...) as.data.frame(list(...), stringsAsFactors = FALSE)
run_E1 <- function(margs = names(MARG_CONFIGS)) {
  out <- list()
  for (mk in margs) {
    cfg <- MARG_CONFIGS[[mk]]; m0 <- MARG[[cfg[1]]]; m1 <- MARG[[cfg[2]]]
    for (fam in names(COPS)) for (rho in RHOS) {
      th <- theta_for_rho(fam, rho, m0, m1)
      for (s1 in SRAT) {
        if (is.na(th)) next
        a <- analyse(true_cop(fam, rho), th, rho, m0, m1, 1, s1)
        out[[length(out) + 1]] <- row_list(model = "copula", family = fam, marg = mk, rho = rho, s1 = s1, theta = th,
                                           cov = a[["cov"]], miss_lo = a[["miss_lo"]], miss_up = a[["miss_up"]],
                                           width_cw = a[["width_cw"]], width_true = a[["width_true"]],
                                           width_naive = a[["width_naive"]], qlo = a[["qlo"]], qhi = a[["qhi"]])
      }
    }
    message("  E1 copula models: ", mk, " done")
  }
  for (lk in names(LAT_CONFIGS)) {
    H <- LAT_LAWS[[LAT_CONFIGS[[lk]][1]]]; E <- LAT_LAWS[[LAT_CONFIGS[[lk]][2]]]
    for (rho in RHOS) {
      lam <- sqrt(abs(rho)); sgn <- if (rho >= 0) 1 else -1; mm <- latent_marginals(lam, sgn, H, E)
      for (s1 in SRAT) {
        b <- cw_oracle(rho, mm[[1]], mm[[2]], 1, s1); Fn <- function(t) latent_ite_cdf(t, lam, sgn, H, E, 1, s1)
        F <- Fn(unname(b)); qlo <- find_quantile(Fn, ALPHA / 2); qhi <- find_quantile(Fn, 1 - ALPHA / 2)
        nv <- naive_oracle(mm[[1]], mm[[2]], 1, s1)
        out[[length(out) + 1]] <- row_list(model = "latent", family = lk, marg = lk, rho = rho, s1 = s1, theta = lam,
                                           cov = F[2] - F[1], miss_lo = F[1], miss_up = 1 - F[2],
                                           width_cw = b[["U"]] - b[["L"]], width_true = qhi - qlo,
                                           width_naive = nv[["U"]] - nv[["L"]], qlo = qlo, qhi = qhi)
      }
    }
    message("  E1 latent models: ", lk, " done")
  }
  do.call(rbind, out)
}

## coverage of CW(rho) at other nominal levels (parametric families)
run_E1_levels <- function(E1) {
  d <- E1[E1$model == "copula" & E1$family %in% REG, ]
  out <- list()
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]; cfg <- MARG_CONFIGS[[r$marg]]; m0 <- MARG[[cfg[1]]]; m1 <- MARG[[cfg[2]]]
    for (al in c(0.20, 0.05, 0.01)) {
      b <- cw_oracle(r$rho, m0, m1, 1, r$s1, al); F <- ite_cdf(unname(b), true_cop(r$family, r$rho), r$theta, m0, m1, 1, r$s1)
      out[[length(out) + 1]] <- row_list(family = r$family, marg = r$marg, rho = r$rho, s1 = r$s1, alpha = al, cov = F[2] - F[1])
    }
  }
  rbind(do.call(rbind, out), row_list(family = d$family, marg = d$marg, rho = d$rho, s1 = d$s1, alpha = 0.10, cov = d$cov))
}

## copula-based intervals of Proposition `asymptotic_theorem_copula` with the true marginals
run_E1_cop <- function(E1) {
  d <- E1[E1$model == "copula", ]; out <- list()
  for (key in unique(paste(d$marg, d$rho, d$s1))) {
    g <- d[paste(d$marg, d$rho, d$s1) == key, ]; q <- g[g$family %in% REG, ]
    if (!("gauss" %in% q$family)) next
    cfg <- MARG_CONFIGS[[g$marg[1]]]; m0 <- MARG[[cfg[1]]]; m1 <- MARG[[cfg[2]]]
    lo_u <- min(q$qlo); hi_u <- max(q$qhi); lo_g <- q$qlo[q$family == "gauss"]; hi_g <- q$qhi[q$family == "gauss"]
    for (i in seq_len(nrow(g))) {
      r <- g[i, ]; F <- ite_cdf(c(lo_u, hi_u, lo_g, hi_g), true_cop(r$family, r$rho), r$theta, m0, m1, 1, r$s1)
      out[[length(out) + 1]] <- row_list(marg = r$marg, rho = r$rho, s1 = r$s1, family = r$family, inside = r$family %in% REG,
                                         cov_cw = r$cov, cov_ga = F[4] - F[3], cov_union = F[2] - F[1],
                                         width_true = r$width_true, width_cw = r$width_cw,
                                         width_ga = hi_g - lo_g, width_union = hi_u - lo_u)
    }
  }
  do.call(rbind, out)
}

## ===================================================================================== E3: worst case (LP)
## Smallest coverage of the oracle CW(rho) over m x m patchwork copulas with Gaussian marginals and
## correlation rho, optionally with cell masses within a factor Gamma of the Gaussian copula.
gauss_cell_masses <- function(r, ze) {
  m <- length(ze) - 1; zq <- seq(-9, 9, length.out = 6001); wq <- dnorm(zq) * (zq[2] - zq[1])
  zc <- pmin(pmax(ze, -38), 38)
  Hc <- pnorm(outer(-r * zq, zc, "+") / sqrt(1 - r^2)) * wq
  cum <- rbind(0, apply(Hc, 2, cumsum))
  Cuv <- cum[findInterval(pmin(pmax(ze, -9), 9), zq) + 1, , drop = FALSE]
  Cuv[1, ] <- 0; Cuv[m + 1, ] <- pnorm(zc); Cuv[, 1] <- 0; Cuv[, m + 1] <- pnorm(zc)
  P <- pmax(t(diff(t(diff(Cuv)))), 0)
  for (it in 1:200) { P <- P * ((1 / m) / rowSums(P)); P <- t(t(P) * ((1 / m) / colSums(P))) }
  P
}
solve_lp <- function(obj, i, j, x, nrow, rhs, lower, upper) {
  if (requireNamespace("highs", quietly = TRUE)) {
    A <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(nrow, length(obj)))
    return(highs::highs_solve(L = obj, lower = lower, upper = upper, A = A, lhs = rhs, rhs = rhs, maximum = TRUE)$objective_value)
  }
  if (!requireNamespace("Rglpk", quietly = TRUE)) stop("Install the 'highs' or the 'Rglpk' package for E3.")
  A <- slam::simple_triplet_matrix(i, j, x, nrow = nrow, ncol = length(obj))
  nv <- length(obj); fin <- is.finite(upper)
  res <- Rglpk::Rglpk_solve_LP(obj, A, rep("==", nrow), rhs, max = TRUE,
                               bounds = list(lower = list(ind = seq_len(nv), val = lower),
                                             upper = list(ind = which(fin), val = upper[fin])))
  res$optimum
}
run_E3 <- function(m = M_LP, k = 8, gammas = c(1, 1.5, 2, Inf)) {
  ze <- qnorm(seq(0, 1, length.out = m + 1)); ze[1] <- -38; ze[m + 1] <- 38
  gbar <- (dnorm(ze[-(m + 1)]) - dnorm(ze[-1])) * m          # exact cell means of a N(0,1) variable
  uu <- as.vector(outer((seq_len(k) - 0.5) / k, 0:(m - 1), "+")) / m
  wsub <- qnorm(uu)                                         # Gaussian marginals: eps = Phi^{-1}(u)
  nv <- m * m; I <- rep(seq_len(m), times = m); J <- rep(seq_len(m), each = m)   # variable (I, J), column-major
  cons_i <- c(I, m + J, rep(2 * m + 1, nv)); cons_j <- rep(seq_len(nv), 3)
  out <- list()
  for (s1 in SRAT) for (rho in RHOS) {
    b <- cw_oracle(rho, m_gauss, m_gauss, 1, s1)
    D <- outer(-wsub, s1 * wsub, "+")                       # s1*eps1 - eps0 on the sub-grid
    a <- apply(array((D < b[["L"]]) | (D > b[["U"]]), dim = c(k, m, k, m)), c(2, 4), mean)
    corr <- function(P) sum(P * outer(gbar, gbar))
    rc <- if (rho == 0) 0 else uniroot(function(r) corr(gauss_cell_masses(r, ze)) - rho,
                                       if (rho > 0) c(rho, min(0.9995, rho + 0.09)) else c(max(-0.9995, rho - 0.09), rho),
                                       tol = 1e-7)$root
    P <- gauss_cell_masses(rc, ze)
    for (G in gammas) {
      val <- if (G == 1) sum(P * a) else {
        lower <- if (is.finite(G)) as.vector(P) / G else rep(0, nv)
        upper <- if (is.finite(G)) pmin(G * as.vector(P), 1 / m) else rep(Inf, nv)
        solve_lp(as.vector(a), cons_i, cons_j, c(rep(1, 2 * nv), as.vector(outer(gbar, gbar))), 2 * m + 1,
                 c(rep(1 / m, 2 * m), rho), lower, upper)
      }
      out[[length(out) + 1]] <- row_list(s1 = s1, rho = rho, Gamma = G, worst_cov = 1 - val)
    }
    message(sprintf("  E3: s1 = %g, rho = %5.2f done", s1, rho))
  }
  do.call(rbind, out)
}

## ===================================================================================== E2: finite samples
E2_COMBOS <- list(
  list(label = "Gaussian / Gaussian",       kind = "copula", fam = "gauss",    marg = "gauss/gauss",     tex = "Gaussian / Gaussian"),
  list(label = "t4 / t3",                   kind = "copula", fam = "t4",       marg = "t3/t3",           tex = "$t_4$ / $t_3$"),
  list(label = "Clayton / Laplace",         kind = "copula", fam = "clayton",  marg = "laplace/laplace", tex = "Clayton / Laplace"),
  list(label = "Gumbel / Uniform",          kind = "copula", fam = "gumbel",   marg = "unif/unif",       tex = "Gumbel / Uniform"),
  list(label = "Frank / Bimodal",           kind = "copula", fam = "frank",    marg = "bimodal/bimodal", tex = "Frank / Bimodal"),
  list(label = "Gauss. mixture / Gaussian", kind = "copula", fam = "gmix",     marg = "gauss/gauss",     tex = "Gauss. mixture / Gaussian"),
  list(label = "Frechet M/W / Gaussian",    kind = "copula", fam = "frechet",  marg = "gauss/gauss",     tex = "Fr\\'echet $M/W$ / Gaussian"),
  list(label = "Gaussian / Log-normal",     kind = "copula", fam = "gauss",    marg = "lognorm/lognorm", tex = "Gaussian / Log-normal"),
  list(label = "Latent: H exp., E Gauss.",  kind = "latent", fam = "H=exp,E=gauss", marg = NA,         tex = "Latent: $H$ exp., $\\tilde\\varepsilon$ Gauss."))
E2_LABELS <- sapply(E2_COMBOS, `[[`, "label")
theta_E1 <- function(E1, fam, mk, rho) {
  h <- E1[E1$model == "copula" & E1$family == fam & E1$marg == mk & abs(E1$rho - rho) < 1e-9, "theta"]
  if (length(h)) h[1] else NA
}
make_tau <- function() {
  a <- rnorm(3); b <- rnorm(3)
  function(x) { s <- (x + 1) * pi / 2; 1 + as.vector(sapply(1:3, function(j) (a[j] * sin(j * s) + b[j] * cos(j * s)) / j) %*% rep(1, 3)) }
}
## ---- CW(rho) / CW^{+CI}(rho): port of `Case study.R` (spline CQR, centred bootstrap CI, c*)
spline_basis <- function(x, n_knots = 6) { q <- unique(quantile(x, seq(0, 1, length.out = n_knots), names = FALSE))
list(lo = q[1], hi = q[length(q)], inner = q[-c(1, length(q))]) }
basis_eval <- function(sb, v) bs(pmin(pmax(v, sb$lo), sb$hi), knots = sb$inner, degree = 3, intercept = TRUE, Boundary.knots = c(sb$lo, sb$hi))
ols_fit <- function(B, y) solve(crossprod(B) + diag(1e-8, ncol(B)), crossprod(B, y))
qr_fit  <- function(B, y, tau) rq.fit(B, y, tau = tau, method = "br")$coefficients
marginal_cqr <- function(x, y, x_new, levels = c(0.90, 0.95), trainprop = 0.8, n_boot = 100) {
  n <- length(y); p <- sample.int(n); ntr <- floor(trainprop * n); tr <- p[1:ntr]; ca <- p[(ntr + 1):n]
  sb <- spline_basis(x[tr]); Btr <- basis_eval(sb, x[tr]); Bca <- basis_eval(sb, x[ca]); Bn <- basis_eval(sb, x_new)
  key <- function(q) sprintf("%.6f", q)
  qs <- sort(unique(round(c((1 - levels) / 2, (1 + levels) / 2), 6)))
  co <- setNames(lapply(qs, function(q) qr_fit(Btr, y[tr], q)), key(qs))
  bc <- ols_fit(Btr, y[tr]); c_c <- as.vector(Bca %*% bc); c_n <- as.vector(Bn %*% bc)
  rs <- abs(y[ca] - c_c); bsig <- ols_fit(Bca, rs)                 # sigma_t: |Y - mu_t| on X, calibration split
  wfun <- function(cc, lo, hi) { fl <- pmax(0.05 * abs(hi - lo), 1e-6); list(l = pmax(cc - lo, fl), u = pmax(hi - cc, fl)) }
  out <- list(m = c_n, sigma = pmax(as.vector(Bn %*% bsig), 0.05 * mean(rs)))
  for (lev in levels) {
    a <- key(round((1 - lev) / 2, 6)); b <- key(round((1 + lev) / 2, 6))
    w <- wfun(c_c, as.vector(Bca %*% co[[a]]), as.vector(Bca %*% co[[b]]))
    sc <- pmax((c_c - y[ca]) / w$l, (y[ca] - c_c) / w$u); kk <- ceiling((length(ca) + 1) * lev)
    gam <- max(sort(sc)[min(kk, length(sc))], 0)
    wn <- wfun(c_n, as.vector(Bn %*% co[[a]]), as.vector(Bn %*% co[[b]]))
    out[[as.character(lev)]] <- list(l = gam * wn$l, u = gam * wn$u)
  }
  ytr <- y[tr]; out$m_star <- t(sapply(seq_len(n_boot), function(bb) { s <- sample.int(ntr, ntr, replace = TRUE)
  as.vector(Bn %*% ols_fit(Btr[s, , drop = FALSE], ytr[s])) }))
  out
}
ci_from_boot <- function(c0, c1, beta = 0.10) {
  q0 <- apply(sweep(c0$m_star, 2, c0$m), 2, quantile, probs = c(beta / 2, 1 - beta / 2))
  q1 <- apply(sweep(c1$m_star, 2, c1$m), 2, quantile, probs = c(beta / 2, 1 - beta / 2))
  list(rl = 2 * pmax(pmax(q1[2, ], -q0[1, ]), 0), ru = 2 * pmax(pmax(-q1[1, ], q0[2, ]), 0))
}
c_star <- function(a, b, r, rho) {
  va <- pmax(a - r / 2, 0); vb <- pmax(b - r / 2, 0)
  cc <- ifelse(r > 0, (Drho(va, vb, rho) + a - va + b - vb - Drho(a, b, rho)) / ifelse(r > 0, r, 1), 0)
  pmin(pmax(cc, 0), 1)
}
## CW^{+CI}_K(rho) of the current paper: D_{rho,K}(a,b) = max{D_rho(a,b), lambda a, lambda b},
## lambda^2 = 1 - 1/K^2 - (2/K)(rho - 1/K)_+ (lambda = 1 if K = Inf), and c = c*(x), so that
## U = tau + r_u + D_{rho,K}((l0 - r_u/2)_+, (u1 - r_u/2)_+), L = tau - r_l - D_{rho,K}((l1 - r_l/2)_+, (u0 - r_l/2)_+).
lambdaK <- function(rho, K) if (is.infinite(K)) 1 else sqrt(max(1 - 1 / K^2 - (2 / K) * max(rho - 1 / K, 0), 0))
DrhoK <- function(a, b, rho, K) pmax(Drho(a, b, rho), lambdaK(rho, K) * pmax(a, b))
## K_hat (Appendix "Choice of K"): 90% quantile over the test points of max{kappa_hat, 1/kappa_hat}, kappa_hat = sigma_1/sigma_0
K_hat <- function(c0, c1) { k <- c1$sigma / c0$sigma; max(1, quantile(pmax(k, 1 / k), 0.9, names = FALSE)) }
cw_interval <- function(tau, w0, w1, rho, ci = NULL, K = Inf) {
  if (is.null(ci)) return(list(L = tau - Drho(w1$l, w0$u, rho), U = tau + Drho(w0$l, w1$u, rho)))
  list(L = tau - ci$rl - DrhoK(pmax(w1$l - ci$rl / 2, 0), pmax(w0$u - ci$rl / 2, 0), rho, K),
       U = tau + ci$ru + DrhoK(pmax(w0$l - ci$ru / 2, 0), pmax(w1$u - ci$ru / 2, 0), rho, K))
}
## ---- copula-based intervals: spline location-scale estimate of F_t, copula draws independent of the data
fit_ls <- function(x, y, x_new) {
  sb <- spline_basis(x); B <- basis_eval(sb, x); beta <- ols_fit(B, y); e <- as.vector(y - B %*% beta)
  bs_ <- ols_fit(B, abs(e)); fl <- 0.05 * mean(abs(e)); s_in <- pmax(as.vector(B %*% bs_), fl); Bn <- basis_eval(sb, x_new)
  list(mu = as.vector(Bn %*% beta), s = pmax(as.vector(Bn %*% bs_), fl), z = sort(e / s_in))
}
qfun <- function(zs) { pp <- (seq_along(zs) - 0.5) / length(zs); function(u) approx(pp, zs, u, rule = 2)$y }
cop_interval <- function(drawsets, tau_n, s0n, s1n, Q0, Q1) {
  r <- s1n / s0n; rg <- if (diff(range(r)) > 1e-9) seq(min(r), max(r), length.out = 9) else min(r)
  lo <- rep(Inf, length(tau_n)); hi <- rep(-Inf, length(tau_n))
  for (ds in drawsets) {
    A <- Q1(ds$U1); B <- Q0(ds$U0)
    qq <- matrix(sapply(rg, function(rr) quantile(rr * A - B, c(ALPHA / 2, 1 - ALPHA / 2), names = FALSE)), nrow = 2)
    ql <- if (length(rg) > 1) approx(rg, qq[1, ], r)$y else rep(qq[1, 1], length(r))
    qh <- if (length(rg) > 1) approx(rg, qq[2, ], r)$y else rep(qq[2, 1], length(r))
    lo <- pmin(lo, tau_n + s0n * ql); hi <- pmax(hi, tau_n + s0n * qh)
  }
  list(L = lo, U = hi)
}
uniforms <- function(cop, th, n) { s <- cop$sample_z(n, th); list(U0 = pnorm(s$zu), U1 = pnorm(s$zv)) }
build_noise <- function(cb, rho, E1) {        # noise sampler + copula draw sets defining Theta
  if (cb$kind == "copula") {
    cfg <- MARG_CONFIGS[[cb$marg]]; m0 <- MARG[[cfg[1]]]; m1 <- MARG[[cfg[2]]]
    th0 <- theta_E1(E1, cb$fam, cb$marg, rho); if (is.na(th0)) return(NULL)
    cop0 <- true_cop(cb$fam, rho)
    sampler <- function(n) { s <- cop0$sample_z(n, th0); list(e0 = m0$g(s$zu), e1 = m1$g(s$zv)) }
    draws <- list(`COP-true` = list(uniforms(cop0, th0, M_DRAWS)),
                  `COP-Ga` = list(uniforms(COPS$gauss, theta_E1(E1, "gauss", cb$marg, rho), M_DRAWS)),
                  `COP-union` = lapply(Filter(function(f) !is.na(theta_E1(E1, f, cb$marg, rho)), REG),
                                       function(f) uniforms(true_cop(f, rho), theta_E1(E1, f, cb$marg, rho), M_DRAWS)))
  } else {
    lc <- LAT_CONFIGS[[cb$fam]]; H <- LAT_LAWS[[lc[1]]]; E <- LAT_LAWS[[lc[2]]]
    lam <- sqrt(abs(rho)); sgn <- if (rho >= 0) 1 else -1; cc <- sqrt(1 - lam^2)
    sampler <- function(n) { h <- H$rvs(n); list(e0 = lam * h + cc * E$rvs(n), e1 = sgn * lam * h + cc * E$rvs(n)) }
    mm <- latent_marginals(lam, sgn, H, E, full = TRUE); s <- sampler(M_DRAWS)
    thg <- if (rho == 0) 0 else theta_for_rho("gauss", rho, mm[[1]], mm[[2]])
    draws <- list(`COP-true` = list(list(U0 = mm[[1]]$cdf(s$e0), U1 = mm[[2]]$cdf(s$e1))),
                  `COP-Ga` = list(uniforms(COPS$gauss, thg, M_DRAWS)))
  }
  list(sampler = sampler, draws = draws)
}
e2_rep <- function(nz) {
  tau <- make_tau(); x <- runif(N_TRAIN + N_TEST, -1, 1); e <- nz$sampler(N_TRAIN + N_TEST)
  y0 <- 5 + 5 * x + e$e0; y1 <- 5 + 5 * x + tau(x) + 2 * e$e1
  tr <- seq_len(N_TRAIN); T <- rbinom(N_TRAIN, 1, (1 + abs(x[tr])) / 4); y <- ifelse(T == 1, y1[tr], y0[tr])
  xa <- x[tr]; xt <- x[-tr]; ite <- (y1 - y0)[-tr]
  c0 <- marginal_cqr(xa[T == 0], y[T == 0], xt); c1 <- marginal_cqr(xa[T == 1], y[T == 1], xt)
  tauh <- c1$m - c0$m; a0 <- c0[["0.9"]]; a1 <- c1[["0.9"]]; n0 <- c0[["0.95"]]; n1 <- c1[["0.95"]]
  Kh <- K_hat(c0, c1)
  ints <- list(CW = cw_interval(tauh, a0, a1, NZRHO), `CW+CI` = cw_interval(tauh, a0, a1, NZRHO, ci_from_boot(c0, c1), K = Kh),
               naive = list(L = tauh - (n1$l + n0$u), U = tauh + (n0$l + n1$u)))
  f0 <- fit_ls(xa[T == 0], y[T == 0], xt); f1 <- fit_ls(xa[T == 1], y[T == 1], xt); Q0 <- qfun(f0$z); Q1 <- qfun(f1$z)
  for (mth in names(nz$draws)) ints[[mth]] <- cop_interval(nz$draws[[mth]], f1$mu - f0$mu, f0$s, f1$s, Q0, Q1)
  do.call(rbind, lapply(names(ints), function(mth) { I <- ints[[mth]]
  row_list(method = mth, coverage = mean(ite >= I$L & ite <= I$U), width = mean(I$U - I$L), Khat = Kh) }))
}
run_E2 <- function(E1) {
  out <- list()
  for (ci in seq_along(E2_COMBOS)) for (ri in seq_along(RHOS)) {
    cb <- E2_COMBOS[[ci]]; rho <- RHOS[ri]
    set.seed(1e5 + 100 * ci + ri); nz <- build_noise(cb, rho, E1); if (is.null(nz)) next
    NZRHO <<- rho
    for (rep in seq_len(R_REPS)) {
      set.seed(1e6 + 1e4 * ci + 100 * ri + rep)
      out[[length(out) + 1]] <- cbind(row_list(combo = cb$label, rho = rho, rep = rep), e2_rep(nz))
    }
    message(sprintf("  E2: %-26s rho = %5.2f done", cb$label, rho))
  }
  do.call(rbind, out)
}

## ===================================================================================== figures and tables
FAM_COL <- c(gauss = "#000000", t4 = "#7f7f7f", frank = "#1f77b4", plackett = "#17becf", clayton = "#2ca02c",
             gumbel = "#ff7f0e", joe = "#d62728", gmix = "#9467bd", frechet = "#8c564b")
FAM_SHP <- c(gauss = 16, t4 = 15, frank = 17, plackett = 6, clayton = 18, gumbel = 3, joe = 4, gmix = 8, frechet = 1)
MARG_LABELS <- sapply(MARG_CONFIGS, `[`, 3)
theme_cw <- function() theme_minimal(base_size = 8.5) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold"))

fig1 <- function(E1, file) {                    # Figure 1: CW(rho), sigma1 = 2 sigma0
  d <- E1[E1$model == "copula" & E1$s1 == 2, ]
  d$panel <- factor(MARG_LABELS[d$marg], levels = MARG_LABELS); d$family <- factor(d$family, levels = names(COPS))
  p <- ggplot(d, aes(rho, cov, colour = family, shape = family)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.88, ymax = 0.92, fill = "grey85", alpha = 0.6) +
    geom_hline(yintercept = 0.9, linewidth = 0.4) + geom_line(linewidth = 0.4) + geom_point(size = 1.3) +
    facet_wrap(~panel, ncol = 3) + coord_cartesian(ylim = c(0.78, 1)) +
    scale_x_continuous(breaks = c(-0.75, 0, 0.75), minor_breaks = RHOS) +
    scale_colour_manual(values = FAM_COL, labels = FAM_LABEL[names(COPS)], name = NULL) +
    scale_shape_manual(values = FAM_SHP, labels = FAM_LABEL[names(COPS)], name = NULL) +
    labs(x = expression("cross-world correlation " * rho), y = expression("coverage of " * CW(rho))) +
    theme_cw() + guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  ggsave(file, p, width = 7.2, height = 6.6)
}

fig2 <- function(E1C, file) {                   # Figure 2: CW(rho) vs copula-based intervals
  ins <- E1C[E1C$inside, ]; out <- list()
  for (mk in names(MARG_CONFIGS)) for (s1 in SRAT) for (v in c("cov_cw", "cov_ga", "cov_union")) {
    q <- ins[ins$marg == mk & ins$s1 == s1, v]
    if (length(q)) out[[length(out) + 1]] <- row_list(marg = mk, s1 = s1, method = v, lo = min(q), hi = max(q))
  }
  d <- do.call(rbind, out)
  d$y <- match(d$marg, rev(names(MARG_CONFIGS))) + c(cov_cw = 0.25, cov_ga = 0, cov_union = -0.25)[d$method]
  d$method <- factor(d$method, levels = c("cov_cw", "cov_ga", "cov_union"))
  d$panel <- factor(paste0("sigma[1]/sigma[0] == ", d$s1))
  p <- ggplot(d) + annotate("rect", xmin = 0.88, xmax = 0.92, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
    geom_vline(xintercept = 0.9, linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = method), linewidth = 2.6) +
    facet_wrap(~panel, labeller = label_parsed) +
    scale_y_continuous(breaks = seq_along(MARG_CONFIGS), labels = rev(MARG_LABELS)) + coord_cartesian(xlim = c(0.8, 0.97)) +
    scale_colour_manual(values = c(cov_cw = "#1f77b4", cov_ga = "#9467bd", cov_union = "#ff7f0e"), name = NULL,
                        labels = c(expression(CW(rho)), expression(C[ITE]^cop * ", Gaussian family"),
                                   expression(C[ITE]^cop * ", union of 7 families"))) +
    labs(x = "oracle coverage", y = NULL) + theme_cw()
  ggsave(file, p, width = 7.2, height = 3.6)
}

pop_ref <- function(E1, cb, rho) {              # population coverage and true width (sigma1 = 2 sigma0)
  r <- if (cb$kind == "copula") E1[E1$model == "copula" & E1$family == cb$fam & E1$marg == cb$marg & abs(E1$rho - rho) < 1e-9 & E1$s1 == 2, ]
  else E1[E1$model == "latent" & E1$family == cb$fam & abs(E1$rho - rho) < 1e-9 & E1$s1 == 2, ]
  if (nrow(r)) c(cov = r$cov[1], width_true = r$width_true[1]) else c(cov = NA, width_true = NA)
}
METH <- c("CW", "CW+CI", "COP-true", "COP-Ga", "COP-union")
fig3 <- function(E2, E1, file) {                # Figure 3: finite samples
  s <- aggregate(coverage ~ combo + rho + method, E2[E2$method %in% METH, ],
                 function(v) c(m = mean(v), lo = unname(quantile(v, 0.1)), hi = unname(quantile(v, 0.9))))
  s <- data.frame(s[, 1:3], s$coverage); s$y <- 6 - match(s$method, METH)
  s$method <- factor(s$method, levels = METH); s$combo <- factor(s$combo, levels = E2_LABELS)
  s$rlab <- factor(paste0("rho == ", s$rho), levels = paste0("rho == ", RHOS))
  grid <- expand.grid(combo = E2_LABELS, rho = RHOS, stringsAsFactors = FALSE)
  grid$pop <- mapply(function(cl, r) pop_ref(E1, E2_COMBOS[[match(cl, E2_LABELS)]], r)[["cov"]], grid$combo, grid$rho)
  grid$has <- paste(grid$combo, grid$rho) %in% paste(s$combo, s$rho)
  grid$combo <- factor(grid$combo, levels = E2_LABELS); grid$rlab <- factor(paste0("rho == ", grid$rho), levels = paste0("rho == ", RHOS))
  grid$x <- 0.88; grid$y <- 3; grid$label <- "Undefined"
  p <- ggplot(s) + geom_vline(xintercept = 0.9, colour = "grey60", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = method), linewidth = 0.5) +
    geom_point(aes(x = m, y = y, colour = method), size = 0.9) +
    geom_segment(data = grid[grid$has & !is.na(grid$pop), ], aes(x = pop, xend = pop, y = 4.6, yend = 5.4), linewidth = 0.5) +
    facet_grid(combo ~ rlab, labeller = labeller(rlab = label_parsed, combo = label_wrap_gen(14))) +
    scale_y_continuous(breaks = NULL, limits = c(0.5, 5.5)) + scale_x_continuous(breaks = c(0.8, 0.9, 1)) +
    coord_cartesian(xlim = c(0.76, 1)) +
    scale_colour_manual(values = c(CW = "#1f77b4", `CW+CI` = "#d62728", `COP-true` = "#2ca02c", `COP-Ga` = "#9467bd", `COP-union` = "#ff7f0e"),
                        labels = c(expression(CW(rho)), expression(CW^{"+CI"} * (rho)), expression(C[ITE]^cop * ", true family"),
                                   expression(C[ITE]^cop * ", Gaussian family"), expression(C[ITE]^cop * ", union of 7")), name = NULL) +
    labs(x = "coverage of Y(1) - Y(0)   (dot: mean over replications; line: 10%-90% range; black tick: population coverage of CW)", y = NULL) +
    theme_cw() + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 6.5), axis.text.x = element_text(size = 6),
                       panel.spacing.x = unit(0.25, "lines"))
  if (any(!grid$has))                           # cells without data (e.g. a family that cannot reach rho)
    p <- p + geom_text(data = grid[!grid$has, ], aes(x = x, y = y, label = label),
                       inherit.aes = FALSE, size = 2.3, fontface = "italic", colour = "grey35")
  ggsave(file, p, width = 7.4, height = 5.6)
}

fmt3 <- function(v) ifelse(is.finite(v), ifelse(v < 0.86, sprintf("\\textbf{%.3f}", v), sprintf("%.3f", v)), "--")
table1 <- function(E1C, file) {                 # Table 1 body (by rho)
  ins <- E1C[E1C$inside, ]; ac <- ins[ins$family != "frechet", ]
  lines <- sapply(RHOS, function(r) {
    q <- ins[ins$rho == r, ]; w <- ac[ac$rho == r, ]; o <- E1C[!E1C$inside & E1C$rho == r, ]
    rng <- function(v) sprintf("%.3f--%.3f", min(v), max(v))
    paste0("$", r, "$ & ", rng(q$cov_cw), " & ", rng(q$cov_ga), " & ", rng(q$cov_union), " & ",
           sprintf("%.3f & %.3f & %.3f", median(w$width_cw / w$width_true), median(w$width_ga / w$width_true),
                   median(w$width_union / w$width_true)), " & ", if (nrow(o)) sprintf("%.3f", min(o$cov_union)) else "--", " \\\\")
  })
  writeLines(lines, file)
}
table2 <- function(E2, E1, file) {              # Table 2 body (per noise model)
  lines <- sapply(seq_along(E2_COMBOS), function(ci) {
    cb <- E2_COMBOS[[ci]]; a <- E2[E2$combo == cb$label, ]
    cov_cells <- sapply(METH, function(m) { q <- tapply(a$coverage[a$method == m], a$rho[a$method == m], mean)
    if (length(q)) sprintf("%s (%s)", fmt3(mean(q)), fmt3(min(q))) else "--" })
    w_cells <- sapply(c("CW+CI", "COP-union", "naive"), function(m) {
      q <- tapply(a$width[a$method == m], a$rho[a$method == m], mean)
      if (!length(q)) return("--")
      wt <- sapply(as.numeric(names(q)), function(r) pop_ref(E1, cb, r)[["width_true"]]); sprintf("%.2f", mean(q / wt)) })
    paste0(cb$tex, " & ", paste(c(cov_cells, w_cells), collapse = " & "), " \\\\")
  })
  writeLines(lines, file)
}

## ===================================================================================== numbers quoted in the text
text_numbers <- function(E1, E1L, E1C, E2, E3, file) {
  L <- character(0); say <- function(...) L <<- c(L, sprintf(...))
  cop <- E1[E1$model == "copula", ]; reg <- cop[cop$family %in% REG, ]; lat <- E1[E1$model == "latent", ]
  rng <- function(v) sprintf("[%.3f, %.3f]", min(v), max(v))
  say("POPULATION LEVEL (%d attainable copula models, %d latent-factor models)", nrow(cop), nrow(lat))
  say("Gaussian copula + Gaussian marginals (ellipticity): coverage %s", rng(cop$cov[cop$family == "gauss" & cop$marg == "gauss/gauss"]))
  say("Gaussian marginals, seven parametric families, all rho: %s", rng(reg$cov[reg$marg == "gauss/gauss"]))
  say("Gaussian marginals, hidden-subgroup mixture and Frechet: %s", rng(cop$cov[cop$marg == "gauss/gauss" & cop$family %in% c("gmix", "frechet")]))
  sp <- aggregate(cov ~ marg + rho + s1, reg, function(v) max(v) - min(v))
  for (mk in names(MARG_CONFIGS)) say("  spread over the seven families, %-16s |rho|<=0.5: %.3f   |rho|=0.75: %.3f", mk,
                                      max(sp$cov[sp$marg == mk & abs(sp$rho) <= 0.5]), max(sp$cov[sp$marg == mk & abs(sp$rho) == 0.75]))
  ind <- cop[cop$family == "frank" & cop$rho == 0 & cop$s1 == 1, ]
  for (mk in c("unif/unif", "bimodal/bimodal", "t3/t3", "laplace/laplace")) say("  independence, sigma1 = sigma0, %-16s: %.3f", mk, ind$cov[ind$marg == mk])
  sk <- reg[reg$marg %in% c("lognorm/lognorm", "exp/exp", "gauss/lognorm"), ]
  say("skewed marginals, |rho| <= 0.5: %s", rng(sk$cov[abs(sk$rho) <= 0.5]))
  for (mk in c("lognorm/lognorm", "exp/exp")) { q <- reg[reg$marg == mk & reg$rho == 0.75 & reg$s1 == 2, ]
  say("  %s, rho = 0.75, sigma1 = 2 sigma0: min coverage %.3f, max lower miscoverage %.3f", mk, min(q$cov), max(q$miss_lo)) }
  say("  skewed (log-normal, exponential), rho = 0.75, sigma1 = sigma0: %s", rng(reg$cov[reg$marg %in% c("lognorm/lognorm", "exp/exp") & reg$rho == 0.75 & reg$s1 == 1]))
  lv <- function(f, s) lat$cov[lat$family == f & lat$rho == 0.75 & lat$s1 == s]
  say("latent, rho = 0.75: exponential H (sigma1 = 2 sigma0) %.3f; t3 H (sigma1 = sigma0) %.3f; two-point H (sigma1 = sigma0) %.3f",
      lv("H=exp,E=gauss", 2), lv("H=t3,E=gauss", 1), lv("H=twopoint,E=gauss", 1))
  ac <- cop[cop$family != "frechet", ]
  say("CW width / true width, median %.3f (IQR %.3f-%.3f)", median(ac$width_cw / ac$width_true),
      quantile(ac$width_cw / ac$width_true, 0.25), quantile(ac$width_cw / ac$width_true, 0.75))
  for (al in c(0.2, 0.1, 0.05, 0.01)) { q <- E1L[E1L$alpha == al, ]
  say("  level %.2f: share with miscoverage <= 1.25 alpha = %.0f%%, median miscoverage/alpha = %.2f", 1 - al, 100 * mean(q$cov >= 1 - 1.25 * al), median((1 - q$cov) / al)) }
  ins <- E1C[E1C$inside, ]; out <- E1C[!E1C$inside, ]; aci <- ins[ins$family != "frechet", ]
  say("copula-based, Gaussian family, |rho| <= 0.5: %s (CW: %s)", rng(ins$cov_ga[abs(ins$rho) <= 0.5]), rng(ins$cov_cw[abs(ins$rho) <= 0.5]))
  lowga <- ins[ins$rho == 0.75, ]; lowga <- lowga[order(lowga$cov_ga), ][1:3, c("marg", "family", "s1", "cov_ga")]
  say("  lowest Gaussian-family coverage at rho = 0.75: %s", paste(sprintf("%s/%s/s1=%g: %.3f", lowga$marg, lowga$family, lowga$s1, lowga$cov_ga), collapse = "; "))
  say("  union of seven families (true copula among them): %s", rng(ins$cov_union))
  say("  median width / true: CW %.3f, Gaussian family %.3f, union %.3f; union / CW %.3f", median(aci$width_cw / aci$width_true),
      median(aci$width_ga / aci$width_true), median(aci$width_union / aci$width_true), median(ins$width_union / ins$width_cw))
  say("  union, true copula outside the family: min %.3f", min(out$cov_union))
  say(""); say("FINITE SAMPLES (baseline design, n = %d, %d replications)", N_TRAIN, R_REPS)
  cell <- aggregate(cbind(coverage, width) ~ combo + rho + method, E2, mean)
  cell$pop <- NA; cell$wt <- NA
  for (i in seq_len(nrow(cell))) { pr <- pop_ref(E1, E2_COMBOS[[match(cell$combo[i], E2_LABELS)]], cell$rho[i]); cell$pop[i] <- pr[["cov"]]; cell$wt[i] <- pr[["width_true"]] }
  cw <- cell[cell$method == "CW", ]
  say("CW minus population coverage: mean %+.4f (range %+.3f to %+.3f)", mean(cw$coverage - cw$pop), min(cw$coverage - cw$pop), max(cw$coverage - cw$pop))
  for (m in c("CW+CI", "COP-true", "COP-Ga", "COP-union", "CW")) { q <- cell[cell$method == m, ]
  say("  %-9s cell means: min %.3f, max %.3f, mean %.3f (%d configurations)", m, min(q$coverage), max(q$coverage), mean(q$coverage), nrow(q)) }
  for (m in c("CW+CI", "COP-union", "naive", "CW", "COP-true")) { q <- cell[cell$method == m, ]; r <- tapply(q$width / q$wt, q$combo, mean)
  say("  width / true, %-9s: per-model averages %.2f-%.2f, median over cells %.3f", m, min(r), max(r), median(q$width / q$wt)) }
  b <- cell[cell$combo == "Frank / Bimodal" & cell$rho == 0, ]
  say("bimodal marginals, rho = 0: CW %.3f, true family %.3f, Gaussian family %.3f", b$coverage[b$method == "CW"],
      b$coverage[b$method == "COP-true"], b$coverage[b$method == "COP-Ga"])
  say(""); say("WORST CASE (patchwork copulas, m = %d)", M_LP)
  say("unrestricted: %s (lowest at rho = %s)", rng(E3$worst_cov[is.infinite(E3$Gamma)]), E3$rho[is.infinite(E3$Gamma)][which.min(E3$worst_cov[is.infinite(E3$Gamma)])])
  for (G in c(1.5, 2)) say("Gamma = %.1f: min %.3f", G, min(E3$worst_cov[E3$Gamma == G]))
  writeLines(L, file); cat(L, sep = "\n")
}

## ===================================================================================== main
if (sys.nframe() == 0L || isTRUE(getOption("cw.run.all", FALSE))) {
  message("E1: population-level sweep"); E1 <- cached("E1", run_E1())
  message("E1: other nominal levels");   E1L <- cached("E1_levels", run_E1_levels(E1))
  message("E1: copula-based intervals"); E1C <- cached("E1_cop", run_E1_cop(E1))
  message("E3: worst case");             E3 <- cached("E3", run_E3())
  message("E2: finite samples");         E2 <- cached("E2", run_E2(E1))
  fig1(E1, file.path(OUT_DIR, "figures", "robust_E1_copulas_sym.pdf"))
  fig2(E1C, file.path(OUT_DIR, "figures", "robust_E1_cop_sym.pdf"))
  fig3(E2, E1, file.path(OUT_DIR, "figures", "robust_E2_summary_sym.pdf"))
  table1(E1C, file.path(OUT_DIR, "tables", "table_E1cop_sym_body.tex"))
  table2(E2, E1, file.path(OUT_DIR, "tables", "table_E2_sym_body.tex"))
  text_numbers(E1, E1L, E1C, E2, E3, file.path(OUT_DIR, "text_numbers.txt"))
  message("Done. Outputs in ", normalizePath(OUT_DIR))
}