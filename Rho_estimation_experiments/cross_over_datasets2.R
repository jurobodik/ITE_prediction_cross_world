###############################################################################
# Additional 19 datasets for the empirical analysis of the
# cross-world correlation  rho = cor(Y(1), Y(0)).
#
# Companion to cross-over_dataset.R (there are 30 replicateBE reference datasets).
# Here we collect further publicly available studies in which every unit is
# exposed to BOTH a control/reference condition and an active condition:
#   * split-body trials (two parts of the same unit treated differently),
#   * clinical crossover trials (mostly placebo- or no-treatment-controlled),
#   * bioavailability crossover studies not contained in replicateBE,
#   * one animal switchback trial,
#   * within-subject behavioural / educational / exercise experiments.

# Reproduces
#   * Table "tab_cross_over_appendix"    : 30 replicateBE bioequivalence datasets
#   * Table "tab_additional_within_unit" : 19 additional studies (24 contrasts)
#   * Table "tab_conditional_rho"        : conditional correlations rho(x)
# and all numbers quoted in the text about these studies (Section 2 and the
# appendix "Cross-over, split-body, and within-subject studies").

# Requirements: R (>= 4.1) with the packages dplyr, tidyr and readr, and an
# internet connection. All data are downloaded from GitHub, pinned to fixed
# commits so that the results cannot change:
#   * replicateBE reference datasets: github.com/Helmut01/replicateBE
#   * all other datasets: read-only GitHub mirrors of CRAN packages
#     (github.com/cran/<package>).
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})
if (!exists("save_tables")) save_tables <- FALSE

# =============================================================================
# 0. Helpers
# =============================================================================

replicateBE_commit <- "08d523145bb6a4c44b7e3f0890cf5f6fd4cd2567"

mirror_commit <- c(
  HistData   = "92bbf20c9cd4cf88ff212f6ae0cffcd0e2e06855",
  ISwR       = "b841213021ddcd8bb0ca07bc511da5a511dae69b",
  Sleuth3    = "d51acd84e9d592762f3c0bbd1ad1ae4cfcba0613",
  Stat2Data  = "776655fd20ff0f65402a557432efe96e881acf10",
  asbio      = "e52cad0a4b12b881a000dc2fb93a93de053aa20d",
  CrossCarry = "2b34cb5d037df535f8a0c4566fdc0482b0bef283",
  daewr      = "8fb3ef99f63c905ee2350ae125ddf6574bebda8a",
  survival   = "40b954baeb1b8ca9c67faaca9f9fdbcd242ee903",
  MASS       = "09cbdbe741b14d88191a14972e7044ed31285746",
  agridat    = "6807cd2f74d7c053077203474c2f2992d1904414",
  Lock5Data  = "af76f4aae2178e6da24e722e8de7cd1a8b0d0228",
  openintro  = "7593731d25102625362da30a655b9f36e7880a21",
  faraway    = "5c64099c3de32aedb00c27c9686ab6fda84805cf",
  afex       = "a05ec093fa9055711481390ac62d8a4f23b9531b",
  PairedData = "7076f81a19dde00a613727ea2a6f7679b782f9eb"
)

# Download data/<file> of a CRAN package from its pinned GitHub mirror and load it.
load_cran_data <- function(pkg, file, object = sub("\\.(rda|RData)$", "", file)) {
  url <- sprintf("https://raw.githubusercontent.com/cran/%s/%s/data/%s",
                 pkg, mirror_commit[[pkg]], file)
  tmp <- file.path(tempdir(), paste0(pkg, "_", file))
  if (!file.exists(tmp)) download.file(url, tmp, mode = "wb", quiet = TRUE)
  if (grepl("\\.(rda|RData)$", file)) {
    env <- new.env()
    load(tmp, envir = env)
    as.data.frame(get(object, envir = env))
  } else {
    read.table(tmp, header = TRUE)  # also reads gzip-compressed .txt.gz files
  }
}

# Pearson correlation with Fisher-z 95% CI, plus Spearman rank correlation.
paired_rho <- function(y0, y1) {
  ok <- is.finite(y0) & is.finite(y1)
  y0 <- y0[ok]; y1 <- y1[ok]
  n  <- length(y0)
  r  <- cor(y1, y0)
  se <- 1 / sqrt(n - 3)
  tibble(n = n, rho = r,
         ci_lower = tanh(atanh(r) - qnorm(0.975) * se),
         ci_upper = tanh(atanh(r) + qnorm(0.975) * se),
         spearman = cor(y1, y0, method = "spearman"))
}

# Remove period main effects estimated within unit (model y ~ unit + trt + period).
remove_period_effects <- function(d) {
  fit <- lm(y ~ factor(unit) + factor(trt) + factor(period), data = d)
  b <- coef(fit)
  b <- b[grepl("^factor\\(period\\)", names(b))]
  b[is.na(b)] <- 0
  key <- paste0("factor(period)", d$period)
  d$y - ifelse(key %in% names(b), b[key], 0)
}

# Average to one value per unit and condition; keep complete pairs.
unit_pairs <- function(d, control, active) {
  d %>%
    filter(trt %in% c(control, active), is.finite(y)) %>%
    group_by(unit, trt) %>%
    summarise(y = mean(y), .groups = "drop") %>%
    pivot_wider(names_from = trt, values_from = y) %>%
    filter(!is.na(.data[[control]]), !is.na(.data[[active]]))
}

# Long data (unit, trt, y [, period]) -> paired estimates (+ period-adjusted rho).
rho_from_long <- function(d, control, active) {
  p <- unit_pairs(d, control, active)
  res <- paired_rho(p[[control]], p[[active]])
  res$rho_period_adj <- NA_real_
  if ("period" %in% names(d) && length(unique(d$period)) > 1) {
    d_adj <- d %>% filter(is.finite(y))
    d_adj$y <- remove_period_effects(d_adj)
    pa <- unit_pairs(d_adj, control, active)
    res$rho_period_adj <- cor(pa[[active]], pa[[control]])
  }
  res
}

wide_to_long <- function(unit, y0, y1, control, active) {
  n_units <- length(unit)
  tibble(unit = rep(unit, 2), trt = rep(c(control, active), each = n_units), y = c(y0, y1))
}

r2 <- function(x) round(x, 2)
r4 <- function(x) round(x, 4)

# =============================================================================
# 1. Thirty replicateBE bioequivalence datasets  (Table tab_cross_over_appendix)
#    Subjects receive the reference (R) and test (T) formulation, possibly
#    repeatedly; Y(0), Y(1) = log PK response under R and T, averaged within
#    subject over replicate periods.
# =============================================================================

estimate_rho_ds <- function(ds_id) {
  url <- sprintf("https://raw.githubusercontent.com/Helmut01/replicateBE/%s/inst/extdata/DS%02d.csv",
                 replicateBE_commit, ds_id)
  tmp <- file.path(tempdir(), sprintf("replicateBE_DS%02d.csv", ds_id))
  if (!file.exists(tmp)) download.file(url, tmp, mode = "wb", quiet = TRUE)
  # missing values are coded as "", "." or "NA" in the reference datasets
  raw <- read_csv(tmp, comment = "#", na = c("", ".", "NA"), show_col_types = FALSE)
  pairs <- raw %>%
    mutate(subject = as.character(subject), treatment = as.character(treatment),
           Y = if ("logPK" %in% names(.)) as.numeric(logPK) else log(as.numeric(PK))) %>%
    filter(treatment %in% c("R", "T"), !is.na(Y)) %>%
    group_by(subject, treatment) %>%
    summarise(Y = mean(Y), .groups = "drop") %>%
    pivot_wider(names_from = treatment, values_from = Y) %>%
    filter(!is.na(R), !is.na(T))
  bind_cols(tibble(dataset = sprintf("%02d", ds_id)), paired_rho(pairs$R, pairs$T))
}

# Origin of each dataset as documented in the replicateBE package
# ("Measured": public-domain data of non-simulated origin; "Simulated": edited
# or derived from simulated datasets).
simulated_ids <- c(7, 8, 9, 12, 13, 14, 15, 18, 19, 20, 22, 25, 27, 29, 30)

rho_tab_BE <- bind_rows(lapply(1:30, estimate_rho_ds)) %>%
  mutate(type = ifelse(as.integer(dataset) %in% simulated_ids, "Simulated", "Measured")) %>%
  arrange(desc(rho)) %>%
  select(dataset, type, n, rho, ci_lower, ci_upper)

# =============================================================================
# 2. Nineteen additional studies  (Table tab_additional_within_unit)
#    Inclusion criteria (fixed before any correlation was computed):
#    (i) each unit receives both a control/reference and an active condition,
#    sequentially or simultaneously on two parts of the same unit; (ii) the
#    condition is assigned by the investigator; (iii) real (non-simulated)
#    individual-level data from a citable source, not part of replicateBE;
#    (iv) at least 10 units with both outcomes. Screened but excluded
#    candidates are listed at the end of this file.
#    Replicates of a condition are averaged within unit; repeated post-dose
#    measurements are averaged (pre-dose values discarded); PK exposure
#    endpoints are log-transformed; multi-arm studies: each active arm vs
#    control.
# =============================================================================

results <- list()
add_result <- function(study, group, contrast, res) {
  results[[length(results) + 1]] <<- bind_cols(
    tibble(study = study, group = group, contrast = contrast), res)
}

## Clinical crossover trials -------------------------------------------------

# Cushny & Peebles (1905): hours of sleep without (control) and with hypnotics.
cp <- load_cran_data("HistData", "CushnyPeebles.RData")
cp$unit <- seq_len(nrow(cp))
cp_long <- cp %>% pivot_longer(-unit, names_to = "trt", values_to = "y")
for (drug in c("L_hyoscyamine", "L_hyoscine", "DL_hyoscine"))
  add_result("Cushny & Peebles (1905)", "Clinical", paste(drug, "vs none"),
             rho_from_long(cp_long, "Control", drug))

# Ashina et al. (1999): L-NMMA vs placebo, headache pain (VAS); grp 1 = placebo first.
ash <- load_cran_data("ISwR", "ashina.txt.gz")
ash$unit <- seq_len(nrow(ash))
ash_long <- bind_rows(
  tibble(unit = ash$unit, trt = "placebo", y = ash$vas.plac,   period = ifelse(ash$grp == 1, 1, 2)),
  tibble(unit = ash$unit, trt = "L-NMMA",  y = ash$vas.active, period = ifelse(ash$grp == 1, 2, 1)))
add_result("Ashina et al. (1999)", "Clinical", "L-NMMA vs placebo",
           rho_from_long(ash_long, "placebo", "L-NMMA"))

# Swain et al. (1990): oat bran vs low-fiber wheat, serum cholesterol.
ob <- load_cran_data("Sleuth3", "case1602.rda")
ob$unit <- seq_len(nrow(ob))
ob_long <- bind_rows(
  tibble(unit = ob$unit, trt = "low fiber",  y = ob$LoFiber, period = ifelse(ob$Order == "LH", 1, 2)),
  tibble(unit = ob$unit, trt = "high fiber", y = ob$HiFiber, period = ifelse(ob$Order == "HL", 1, 2)))
add_result("Swain et al. (1990)", "Clinical", "Oat bran vs low fiber",
           rho_from_long(ob_long, "low fiber", "high fiber"))

# Chang et al. (1979): THC vs placebo, vomiting episodes.
mj <- load_cran_data("Sleuth3", "ex0432.rda")
add_result("Chang et al. (1979)", "Clinical", "THC vs placebo",
           rho_from_long(wide_to_long(mj$Subject, mj$Placebo, mj$Marijuana, "placebo", "THC"),
                         "placebo", "THC"))

# Bradstreet & Panebianco (2004): drug D vs placebo, log AUC of ethinyl estradiol.
oc <- load_cran_data("Stat2Data", "Contraceptives.rda")
oc_long <- tibble(unit = oc$ID, trt = as.character(oc$Treatment), y = log(oc$EE), period = oc$StudyPeriod)
add_result("Bradstreet & Panebianco (2004)", "Clinical", "Drug D vs placebo",
           rho_from_long(oc_long, "Placebo", "Drug"))

# Littell et al. (2002): asthma drugs a, c vs placebo p; FEV1 averaged over 8 hours.
fev <- load_cran_data("asbio", "asthma.rda")
fev_long <- tibble(unit = fev$PATIENT, trt = as.character(fev$DRUG),
                   y = rowMeans(fev[, paste0("FEV1", 1:8, "H")]))
for (drug in c("a", "c"))
  add_result("Littell et al. (2002)", "Clinical", paste("Drug", drug, "vs placebo"),
             rho_from_long(fev_long, "p", drug))

# Jones & Kenward (2015): 20 mg (B), 40 mg (C) vs placebo (A); post-dose diastolic BP.
art <- load_cran_data("CrossCarry", "Arterial.rda")
art_long <- art %>% filter(Time > 0) %>%
  group_by(unit = Subject, trt = Treatment, period = Period) %>%
  summarise(y = mean(Pressure), .groups = "drop")
for (dose in c("B", "C"))
  add_result("Jones & Kenward (2015)", "Clinical", c(B = "20 mg vs placebo", C = "40 mg vs placebo")[[dose]],
             rho_from_long(art_long, "A", dose))

## Bioavailability crossover studies -----------------------------------------

# Lawson (2015): two-period crossover of an antifungal agent (response 'pl').
af <- load_cran_data("daewr", "antifungal.rda")
af_long <- tibble(unit = af$Subject, trt = as.character(af$Treat), y = af$pl, period = af$Period)
add_result("Lawson (2015), antifungal", "Bioavailability", "Antifungal B vs A",
           rho_from_long(af_long, "A", "B"))

# Lawson (2015): extra-period (ABB/BAA) bioequivalence study, log scale.
be <- load_cran_data("daewr", "bioequiv.rda")
be_long <- tibble(unit = be$Subject, trt = as.character(be$Treat), y = log(be$y), period = be$Period)
add_result("Lawson (2015), bioequivalence", "Bioavailability", "Formulation B vs A",
           rho_from_long(be_long, "A", "B"))

## Split-body trials ---------------------------------------------------------

# Diabetic Retinopathy Study: one eye lasered, fellow eye untreated; blindness
# within 24 months (status known if event before 24 months or followed beyond).
drs <- load_cran_data("survival", "diabetic.rda")
drs_long <- drs %>%
  mutate(y = case_when(status == 1 & time <= 24 ~ 1, time > 24 ~ 0, TRUE ~ NA_real_),
         trt = ifelse(trt == 1, "laser", "untreated")) %>%
  select(unit = id, trt, y)
add_result("Diabetic Retinopathy Study", "Split-body", "Laser vs untreated",
           rho_from_long(drs_long, "untreated", "laser"))

# Box, Hunter & Hunter (1978): sole materials A, B on the two feet of 10 boys.
sh <- load_cran_data("MASS", "shoes.rda")
add_result("Box, Hunter & Hunter (1978)", "Split-body", "Material B vs A",
           rho_from_long(wide_to_long(1:10, sh$A, sh$B, "A", "B"), "A", "B"))

## Animal switchback trial ---------------------------------------------------

# Brandt (1938): two grain rations, T1-T2-T1 / T2-T1-T2 switchback; milk yield.
bs <- load_cran_data("agridat", "brandt.switchback.txt")
bs_long <- tibble(unit = bs$cow, trt = bs$trt, y = bs$yield, period = bs$period)
add_result("Brandt (1938)", "Animal", "Ration T2 vs T1", rho_from_long(bs_long, "T1", "T2"))

## Within-subject experiments ------------------------------------------------

# Chard et al. (2019): drinking water vs usual school day; letter-cancellation
# score LCS (= correct - incorrect, range -14 to 38).
wa <- load_cran_data("CrossCarry", "Water.rda")
wa_long <- tibble(unit = wa$ID, trt = ifelse(wa$Treatment == 1, "water", "control"),
                  y = wa$LCS, period = wa$Period)
add_result("Chard et al. (2019)", "Within-subject", "Water vs usual day",
           rho_from_long(wa_long, "control", "water"))

# Varao-Sousa & Kingstone (2015): live vs video lecture; memory test score.
ll <- load_cran_data("openintro", "lecture_learning.rda")
ll_long <- tibble(unit = ll$student, trt = ll$method, y = ll$memory)
add_result("Varao-Sousa & Kingstone (2015)", "Within-subject", "Video vs live lecture",
           rho_from_long(ll_long, "Live", "Video"))

# Lin et al. (2020): high- vs low-demand task, mean Stroop reaction time over
# all non-missing trials (aggregation as in the afex vignettes).
st <- load_cran_data("afex", "stroop.rda")
st_long <- st %>% filter(!is.na(rt)) %>%
  transmute(unit = as.character(pno), trt = as.character(condition), y = rt)
add_result("Lin et al. (2020)", "Within-subject", "High vs low demand",
           rho_from_long(st_long, "control", "deplete"))

# Steele (1998): foam vs normal surface; CTSIB balance score averaged over
# vision conditions and replicates.
ct <- load_cran_data("faraway", "ctsib.rda")
ct_long <- tibble(unit = ct$Subject, trt = as.character(ct$Surface), y = ct$CTSIB)
add_result("Steele (1998)", "Within-subject", "Foam vs normal surface",
           rho_from_long(ct_long, "norm", "foam"))

# de Lucas et al. (2000): wetsuit vs swimsuit; 1500 m swim velocity.
ws <- load_cran_data("Lock5Data", "Wetsuits.rda")
add_result("de Lucas et al. (2000)", "Within-subject", "Wetsuit vs swimsuit",
           rho_from_long(wide_to_long(seq_len(nrow(ws)), ws$NoWetsuit, ws$Wetsuit,
                                      "no wetsuit", "wetsuit"), "no wetsuit", "wetsuit"))

# Chipman, in Lawson (2015): Williams design, three sprint conditions (1 = reference).
cs <- load_cran_data("daewr", "chipman.rda")
cs_long <- tibble(unit = cs$Subject, trt = as.character(cs$Treat), y = cs$Time, period = cs$Period)
for (k in c("2", "3"))
  add_result("Lawson (2015), sprint", "Within-subject", paste("Condition", k, "vs 1"),
             rho_from_long(cs_long, "1", k))

# Wetherill (1972): blink rate, oscillating vs straight track.
bl <- load_cran_data("PairedData", "Blink.rda")
add_result("Wetherill (1972)", "Within-subject", "Oscillating vs straight",
           rho_from_long(wide_to_long(bl$Subject, bl$Straight, bl$Oscillating,
                                      "Straight", "Oscillating"), "Straight", "Oscillating"))

rho_tab_additional <- bind_rows(results)

# =============================================================================
# 3. Conditional correlations rho(x)  (Table tab_conditional_rho)
#    Studies of Section 2 that record pre-treatment covariates X (the replicateBE
#    datasets record none; post-treatment variables are not used).
#    rho(x) is estimated by the partial correlation of Y(1) and Y(0) given X:
#    the correlation of the residuals of arm-specific linear regressions on X,
#    with Fisher-z CI on n - 3 - k degrees of freedom (k covariate columns).
#    Strata: correlations within the levels of one covariate (baseline values
#    split at the median).
# =============================================================================

partial_rho <- function(y0, y1, X) {
  X  <- as.data.frame(X)
  ok <- is.finite(y0) & is.finite(y1) & complete.cases(X)
  y0 <- y0[ok]; y1 <- y1[ok]; X <- X[ok, , drop = FALSE]
  e0 <- resid(lm(y0 ~ ., data = cbind(y0 = y0, X)))
  e1 <- resid(lm(y1 ~ ., data = cbind(y1 = y1, X)))
  n  <- length(y0); k <- ncol(model.matrix(~ ., data = X)) - 1
  r  <- cor(e1, e0); se <- 1 / sqrt(n - 3 - k)
  tibble(n = n, k = k, rho_marginal = cor(y1, y0), rho_x = r,
         ci_lower = tanh(atanh(r) - qnorm(0.975) * se),
         ci_upper = tanh(atanh(r) + qnorm(0.975) * se))
}
strata_rho <- function(y0, y1, g) {
  bind_rows(lapply(sort(unique(g)), function(l) {
    i <- which(g == l)
    bind_cols(tibble(stratum = as.character(l)), paired_rho(y0[i], y1[i]))
  }))
}
median_split <- function(x) ifelse(x <= median(x), "below", "above")

adj <- list(); strat <- list()
add_cond <- function(study, contrast, covariates, y0, y1, X, g = NULL) {
  adj[[length(adj) + 1]] <<- bind_cols(tibble(study = study, contrast = contrast,
                                              covariates = covariates), partial_rho(y0, y1, X))
  if (!is.null(g))
    strat[[length(strat) + 1]] <<- bind_cols(tibble(study = study, contrast = contrast),
                                             strata_rho(y0, y1, g))
}

# Swain et al. (1990): baseline cholesterol.
add_cond("Swain et al. (1990)", "Oat bran vs low fiber", "Baseline cholesterol",
         ob$LoFiber, ob$HiFiber, data.frame(baseline = ob$Baseline), median_split(ob$Baseline))

# Littell et al. (2002): baseline FEV1 (patient mean of the pre-dose values).
fev_base <- fev %>% group_by(unit = PATIENT) %>% summarise(base = mean(BASEFEV1))
for (drug in c("a", "c")) {
  p <- unit_pairs(fev_long, "p", drug) %>% left_join(fev_base, by = "unit")
  add_cond("Littell et al. (2002)", paste("Drug", drug, "vs placebo"), "Baseline FEV1",
           p$p, p[[drug]], data.frame(base = p$base), median_split(p$base))
}

# Jones & Kenward (2015): baseline diastolic pressure (mean of pre-dose values).
art_base <- art %>% filter(Time <= 0) %>% group_by(unit = Subject) %>% summarise(base = mean(Pressure))
for (dose in c("B", "C")) {
  p <- unit_pairs(art_long, "A", dose) %>% left_join(art_base, by = "unit")
  add_cond("Jones & Kenward (2015)", c(B = "20 mg vs placebo", C = "40 mg vs placebo")[[dose]],
           "Baseline pressure", p$A, p[[dose]], data.frame(base = p$base))
}

# Diabetic Retinopathy Study: type of diabetes, laser, treated side, age, risk scores.
ret <- load_cran_data("survival", "retinopathy.rda")
ret$Y <- ifelse(ret$status == 1 & ret$futime <= 24, 1, ifelse(ret$futime > 24, 0, NA))
dd <- merge(ret[ret$trt == 1, c("id", "type", "laser", "eye", "age", "risk", "Y")],
            ret[ret$trt == 0, c("id", "risk", "Y")], by = "id", suffixes = c("_t", "_c"))
dd <- dd[!is.na(dd$Y_t) & !is.na(dd$Y_c), ]
add_cond("Diabetic Retinopathy Study", "Laser vs untreated", "Type, laser, side, age, risk scores",
         dd$Y_c, dd$Y_t, dd[, c("type", "laser", "eye", "age", "risk_t", "risk_c")],
         as.character(dd$type))

# Chard et al. (2019): school. The treatment sequence is determined by school
# (A: control first, B: water first), so conditioning on school also removes
# the period (practice) effect. Age and sex are missing for many pupils.
wa_cov <- wa %>% group_by(unit = ID) %>% summarise(school = first(school))
p <- unit_pairs(wa_long, "control", "water") %>% left_join(wa_cov, by = "unit")
add_cond("Chard et al. (2019)", "Water vs usual day", "School",
         p$control, p$water, data.frame(school = p$school), as.character(p$school))

# Varao-Sousa & Kingstone (2015): gender (interest and motivation are post-treatment).
ll_cov <- ll %>% group_by(unit = student) %>% summarise(gender = first(gender))
p <- unit_pairs(ll_long, "Live", "Video") %>% left_join(ll_cov, by = "unit")
add_cond("Varao-Sousa & Kingstone (2015)", "Video vs live lecture", "Gender",
         p$Live, p$Video, data.frame(gender = p$gender), as.character(p$gender))

# Lin et al. (2020): constituent study (1-4).
st_cov <- st %>% distinct(unit = as.character(pno), study)
p <- unit_pairs(st_long, "control", "deplete") %>% left_join(st_cov, by = "unit")
add_cond("Lin et al. (2020)", "High vs low demand", "Study",
         p$control, p$deplete, data.frame(study = p$study), as.character(p$study))

# Steele (1998): sex, age, height, weight.
ct_cov <- ct %>% group_by(unit = Subject) %>%
  summarise(sex = first(Sex), age = first(Age), height = first(Height), weight = first(Weight))
p <- unit_pairs(ct_long, "norm", "foam") %>% left_join(ct_cov, by = "unit")
add_cond("Steele (1998)", "Foam vs normal surface", "Sex, age, height, weight",
         p$norm, p$foam, p[, c("sex", "age", "height", "weight")], as.character(p$sex))

# de Lucas et al. (2000): gender, swimmer type.
add_cond("de Lucas et al. (2000)", "Wetsuit vs swimsuit", "Gender, swimmer type",
         ws$NoWetsuit, ws$Wetsuit, ws[, c("Gender", "Type")])

rho_tab_conditional <- bind_rows(adj)
rho_tab_strata      <- bind_rows(strat)

# =============================================================================
# 4. Output: tables and numbers quoted in the text
# =============================================================================

options(width = 200, pillar.sigfig = 4)

cat("\n==== Table tab_cross_over_appendix: 30 replicateBE datasets ====\n")
print(as.data.frame(rho_tab_BE %>% mutate(across(c(rho, ci_lower, ci_upper), r4))), row.names = FALSE)
cat("Average rho:", r4(mean(rho_tab_BE$rho)), "\n")

cat("\n==== Table tab_additional_within_unit: 19 studies, 24 contrasts ====\n")
print(as.data.frame(rho_tab_additional %>%
        transmute(study, contrast, n, rho = r2(rho), ci_lower = r2(ci_lower),
                  ci_upper = r2(ci_upper), spearman = r2(spearman))), row.names = FALSE)

cat("\n==== Table tab_conditional_rho: conditional correlations rho(x) ====\n")
strata_summary <- rho_tab_strata %>%
  group_by(study, contrast) %>%
  summarise(strata = if (n() > 2) paste0(sprintf("%.2f", min(rho)), "-", sprintf("%.2f", max(rho)))
                     else paste(sprintf("%.2f", rho), collapse = ", "), .groups = "drop")
print(as.data.frame(rho_tab_conditional %>%
        left_join(strata_summary, by = c("study", "contrast")) %>%
        transmute(study, contrast, covariates, n, rho = r2(rho_marginal), rho_x = r2(rho_x),
                  ci_lower = r2(ci_lower), ci_upper = r2(ci_upper),
                  strata = ifelse(is.na(strata), "--", strata))), row.names = FALSE)
cat("Strata order: baseline (above, below the median), type of diabetes (adult, juvenile),",
    "school (A, B), gender or sex (female, male), study (range over 1-4).\n")

all_rho <- c(rho_tab_BE$rho, rho_tab_additional$rho)
all_lo  <- c(rho_tab_BE$ci_lower, rho_tab_additional$ci_lower)
lowest  <- rho_tab_additional[which.min(rho_tab_additional$rho), ]

cat("\n==== Numbers quoted in the text ====\n")
cat("Studies:", nrow(rho_tab_BE) + n_distinct(rho_tab_additional$study),
    "(30 replicateBE + 19 additional); contrasts:", length(all_rho), "\n")
cat("All 54 estimates positive:", all(all_rho > 0),
    "| median:", r2(median(all_rho)),
    "| > 0.5:", sum(all_rho > 0.5), "| > 0.33:", sum(all_rho > 0.33),
    "| CI excludes 0:", sum(all_lo > 0), "\n")
cat("Lowest:", r2(lowest$rho), "+/-", r2((lowest$ci_upper - lowest$ci_lower) / 2),
    "(", lowest$study, ")\n")
cat("replicateBE: > 0.5:", sum(rho_tab_BE$rho > 0.5), "of 30 | > 0.33:", sum(rho_tab_BE$rho > 0.33),
    "of 30 | values below 0.5 all simulated:", all(rho_tab_BE$type[rho_tab_BE$rho < 0.5] == "Simulated"),
    "| minimum over measured datasets:", r4(min(rho_tab_BE$rho[rho_tab_BE$type == "Measured"])), "\n")
cat("Additional studies: median", r2(median(rho_tab_additional$rho)),
    "| > 0.5:", sum(rho_tab_additional$rho > 0.5), "of 24 | CI excludes 0:",
    sum(rho_tab_additional$ci_lower > 0), "of 24\n")

# Attainable correlation of the binary retinopathy outcomes, given their marginals.
drs_pairs <- unit_pairs(drs_long, "untreated", "laser")
p0 <- mean(drs_pairs$untreated); p1 <- mean(drs_pairs$laser)
lo <- min(p0, p1); hi <- max(p0, p1)
cat(sprintf("Retinopathy: blindness %.1f%% (untreated) vs %.1f%% (laser); maximal attainable correlation %.2f\n",
            100 * p0, 100 * p1, sqrt(lo * (1 - hi) / (hi * (1 - lo)))))

pa <- rho_tab_additional %>% filter(!is.na(rho_period_adj))
cat("Period adjustment: contrasts with period labels:", nrow(pa),
    "| rho increased in all:", all(pa$rho_period_adj > pa$rho), "\n")
print(as.data.frame(pa %>% transmute(study, contrast, rho = r2(rho), rho_period_adjusted = r2(rho_period_adj))),
      row.names = FALSE)

cat("Conditional: adjusted rho(x) positive:", sum(rho_tab_conditional$rho_x > 0), "of", nrow(rho_tab_conditional),
    "| CI excludes 0:", sum(rho_tab_conditional$ci_lower > 0),
    "| stratum-specific positive:", sum(rho_tab_strata$rho > 0), "of", nrow(rho_tab_strata), "\n")

if (save_tables) {
  write.csv(rho_tab_BE, "table_replicateBE.csv", row.names = FALSE)
  write.csv(rho_tab_additional, "table_additional_studies.csv", row.names = FALSE)
  write.csv(rho_tab_conditional, "table_conditional_rho.csv", row.names = FALSE)
  write.csv(rho_tab_strata, "table_conditional_rho_strata.csv", row.names = FALSE)
}

# =============================================================================
# Screened but excluded candidates (reason in brackets)
# =============================================================================
# Fewer than 10 complete pairs: agridat::lucas.switchback and
#   agridat::patterson.switchback [each cow receives only 2 of the 3/4
#   treatments]; nlme::Tetracycline1/2 [5 subjects]; daewr::bioeqv [3];
#   PairedData::Tobacco and ::Sewage [8]; BSDA::Coffee [9]; nlme::ergoStool [9];
#   nlme::Machines [6]; nlme::Meat [<= 3 per contrast].
# Not a within-unit contrast: agridat::depalluel.sheep [one diet per animal];
#   exactRankTests::ears [treatment assigned per child]; heplots::Headache,
#   Lock5Data::CaffeineTaps [between-subject].
# Simulated / hypothetical: BE::NCAResult4BE, Lock5Data::TrafficFlow,
#   heplots::ReactTime and afex::md_12.1 [hypothetical data].
# Duplicate of a replicateBE dataset: BElikelihood::dat [same Cmax values as DS26].
# Only within-subject differences available: Stat2Data::Eyes, datasets::sleep.
# Matched pairs of distinct units: HistData::ZeaMays, PairedData::Corn,
#   DAAG::mignonette, PairedData::ChickWeight, PairedData::BloodLead,
#   PairedData::Barley, ::Grain, ::Grain2, Sleuth3::case0202 [twins].
# Method comparison / no intervention: PairedData::Iron, ::Meat, PASWR2::TENNIS,
#   Lock5Data::QuizPulse10, PairedData::GrapeFruit, PairedData::Anorexia
#   [before-after], product ratings (Lock5Data::WaterTaste,
#   crossdes::biertest.dat, Lock5Data::StorySpoilers).
