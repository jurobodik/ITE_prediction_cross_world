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
#
# Data are downloaded from the read-only GitHub mirrors of the CRAN packages
# (https://github.com/cran/<package>), pinned to fixed commits.

# Inclusion criteria (fixed before any correlation was computed):
#   (i)   each unit receives both a control/reference and an active condition,
#         either sequentially (crossover / switchback) or simultaneously on two
#         parts of the same unit (split-body);
#   (ii)  the condition is an intervention assigned by the investigator
#         (not a measurement method, a naturally occurring exposure, or time);
#   (iii) real (non-simulated, non-hypothetical) individual-level data from a
#         citable source that do not duplicate a replicateBE dataset;
#   (iv)  at least 10 units with both outcomes observed.
# Screened but excluded candidates are listed at the bottom of this file.
#
# Analysis (as in cross-over_dataset.R):
#   * construct unit-level pairs (Y0, Y1); replicate measurements of the same
#     condition are averaged within unit, repeated post-dose measurements over
#     time are averaged (pre-dose / baseline measurements are discarded);
#   * pharmacokinetic exposure endpoints (AUC) on log scale, all others raw;
#   * Pearson correlation with Fisher-z 95% confidence interval;
#   * multi-arm studies: every active arm versus the control arm;
#   * binary outcome (diabetic retinopathy): Pearson (phi) correlation of
#     blindness within a pre-specified horizon of 24 months;
#   * robustness: Spearman correlation, and Pearson after removing period
#     effects (estimated within subject) whenever period labels are available.
###############################################################################

library(dplyr)
library(tidyr)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

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

# Download data/<file> of a CRAN package from its GitHub mirror and load it.
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

# Pearson correlation with Fisher-z CI (identical to cross-over_dataset.R),
# plus the Spearman rank correlation as a robustness check.
paired_rho <- function(y0, y1) {
  ok <- is.finite(y0) & is.finite(y1)
  y0 <- y0[ok]; y1 <- y1[ok]
  n  <- length(y0)
  r  <- cor(y1, y0)
  se <- 1 / sqrt(n - 3)
  tibble(
    n        = n,
    rho      = r,
    ci_lower = tanh(atanh(r) - qnorm(0.975) * se),
    ci_upper = tanh(atanh(r) + qnorm(0.975) * se),
    spearman = cor(y1, y0, method = "spearman")
  )
}

# Remove period main effects, estimated within unit from
# y ~ unit + treatment + period (the standard crossover model).
remove_period_effects <- function(d) {
  fit <- lm(y ~ factor(unit) + factor(trt) + factor(period), data = d)
  b <- coef(fit)
  b <- b[grepl("^factor\\(period\\)", names(b))]
  b[is.na(b)] <- 0
  key <- paste0("factor(period)", d$period)
  d$y - ifelse(key %in% names(b), b[key], 0)
}

# Average to one value per unit and condition, keep complete pairs.
unit_pairs <- function(d, control, active) {
  d %>%
    filter(trt %in% c(control, active), is.finite(y)) %>%
    group_by(unit, trt) %>%
    summarise(y = mean(y), .groups = "drop") %>%
    pivot_wider(names_from = trt, values_from = y) %>%
    filter(!is.na(.data[[control]]), !is.na(.data[[active]]))
}

# Long data (columns unit, trt, y and optionally period) -> one result row.
rho_from_long <- function(d, control, active) {
  p <- unit_pairs(d, control, active)
  res <- paired_rho(p[[control]], p[[active]])
  res$rho_period_adj <- NA_real_
  if ("period" %in% names(d) && length(unique(d$period)) > 1) {
    d_adj <- d %>% filter(is.finite(y))
    d_adj$y <- remove_period_effects(d_adj)
    p <- unit_pairs(d_adj, control, active)
    res$rho_period_adj <- cor(p[[active]], p[[control]])
  }
  res
}

wide_to_long <- function(unit, y0, y1, control = "control", active = "active") {
  n_units <- length(unit)
  tibble(unit = rep(unit, 2),
         trt  = rep(c(control, active), each = n_units),
         y    = c(y0, y1))
}

results <- list()
add_result <- function(study, domain, design, outcome, control, active, res) {
  results[[length(results) + 1]] <<- bind_cols(
    tibble(study = study, domain = domain, design = design, outcome = outcome,
           contrast = paste(active, "vs", control)), res)
}

# ---------------------------------------------------------------------------
# (A) Clinical crossover trials with placebo / no-treatment control
# ---------------------------------------------------------------------------

## A1. Cushny & Peebles (1905): hours of sleep without hypnotic (control) and
##     with three hyoscine/hyoscyamine hydrobromides; averages over nights.
##     Patient 11 has no control observation.
cp <- load_cran_data("HistData", "CushnyPeebles.RData")
cp$unit <- seq_len(nrow(cp))
cp_long <- cp %>% pivot_longer(-unit, names_to = "trt", values_to = "y")
for (drug in c("L_hyoscyamine", "L_hyoscine", "DL_hyoscine"))
  add_result("Cushny & Peebles (1905)", "Clinical", "within-patient, alternate nights",
             "hours of sleep", "Control", drug, rho_from_long(cp_long, "Control", drug))

## A2. Ashina et al. (1999): NO-synthase inhibitor L-NMMA vs placebo in chronic
##     tension-type headache; summed change of pain VAS from baseline.
##     grp = 1: placebo first, grp = 2: active first.
ash <- load_cran_data("ISwR", "ashina.txt.gz")
ash$unit <- seq_len(nrow(ash))
ash_long <- bind_rows(
  tibble(unit = ash$unit, trt = "placebo", y = ash$vas.plac,   period = ifelse(ash$grp == 1, 1, 2)),
  tibble(unit = ash$unit, trt = "L-NMMA",  y = ash$vas.active, period = ifelse(ash$grp == 1, 2, 1)))
add_result("Ashina et al. (1999)", "Clinical", "2-period crossover",
           "pain score (VAS change)", "placebo", "L-NMMA", rho_from_long(ash_long, "placebo", "L-NMMA"))

## A3. Swain et al. (1990): high-fiber oat bran vs low-fiber wheat supplement;
##     total serum cholesterol after each 6-week period.
ob <- load_cran_data("Sleuth3", "case1602.rda")
ob$unit <- seq_len(nrow(ob))
ob_long <- bind_rows(
  tibble(unit = ob$unit, trt = "low fiber",  y = ob$LoFiber, period = ifelse(ob$Order == "LH", 1, 2)),
  tibble(unit = ob$unit, trt = "high fiber", y = ob$HiFiber, period = ifelse(ob$Order == "HL", 1, 2)))
add_result("Swain et al. (1990)", "Clinical", "2-period crossover",
           "serum cholesterol", "low fiber", "high fiber", rho_from_long(ob_long, "low fiber", "high fiber"))

## A4. Chang et al. (1979): THC (marijuana) vs placebo during chemotherapy;
##     number of vomiting and retching episodes (period labels not available).
mj <- load_cran_data("Sleuth3", "ex0432.rda")
add_result("Chang et al. (1979)", "Clinical", "2-period crossover",
           "vomiting episodes", "placebo", "THC",
           rho_from_long(wide_to_long(mj$Subject, mj$Placebo, mj$Marijuana, "placebo", "THC"), "placebo", "THC"))

## A5. Bradstreet & Panebianco (2004): drug D vs placebo co-administered with an
##     oral contraceptive; ethinyl estradiol AUC (log scale).
oc <- load_cran_data("Stat2Data", "Contraceptives.rda")
oc_long <- tibble(unit = oc$ID, trt = as.character(oc$Treatment), y = log(oc$EE), period = oc$StudyPeriod)
add_result("Bradstreet & Panebianco (2004)", "Clinical", "2-period crossover",
           "log AUC of ethinyl estradiol", "Placebo", "Drug", rho_from_long(oc_long, "Placebo", "Drug"))

## A6. Littell et al. (2002): two asthma drugs (a, c) and placebo (p);
##     FEV1 averaged over the 8 hourly post-dose measurements
##     (the 11-h pre-dose baseline is discarded). Period labels not available.
##     Within-patient pairing is supported by the data: baseline FEV1 of the
##     same patient correlates 0.85-0.89 across the three treatment periods.
fev <- load_cran_data("asbio", "asthma.rda")
fev_long <- tibble(unit = fev$PATIENT, trt = as.character(fev$DRUG),
                   y = rowMeans(fev[, paste0("FEV1", 1:8, "H")]))
for (drug in c("a", "c"))
  add_result("Littell et al. (2002)", "Clinical", "3-period crossover",
             "FEV1 (mean over 8 h)", "p", drug, rho_from_long(fev_long, "p", drug))

## A7. Jones & Kenward (2015): placebo (A) vs 20 mg (B) and 40 mg (C) of a test
##     drug; diastolic blood pressure averaged over post-dose times (Time > 0).
art <- load_cran_data("CrossCarry", "Arterial.rda")
art_long <- art %>% filter(Time > 0) %>%
  group_by(unit = Subject, trt = Treatment, period = Period) %>%
  summarise(y = mean(Pressure), .groups = "drop")
for (dose in c("B", "C"))
  add_result("Jones & Kenward (2015)", "Clinical", "3-period crossover",
             "diastolic blood pressure", "A", dose, rho_from_long(art_long, "A", dose))

# ---------------------------------------------------------------------------
# (B) Bioavailability crossover studies (not part of replicateBE)
# ---------------------------------------------------------------------------

## B1. Lawson (2015, Ch. 9): two-period crossover of an antifungal agent
##     (treatments A, B); response 'pl' on its original scale.
af <- load_cran_data("daewr", "antifungal.rda")
af_long <- tibble(unit = af$Subject, trt = as.character(af$Treat), y = af$pl, period = af$Period)
add_result("Lawson (2015), antifungal", "Bioavailability", "2-period crossover",
           "drug response", "A", "B", rho_from_long(af_long, "A", "B"))

## B2. Lawson (2015, Ch. 9): extra-period (ABB/BAA) bioequivalence study;
##     log scale, replicates of the same formulation averaged within subject.
be <- load_cran_data("daewr", "bioequiv.rda")
be_long <- tibble(unit = be$Subject, trt = as.character(be$Treat), y = log(be$y), period = be$Period)
add_result("Lawson (2015), bioequivalence", "Bioavailability", "3-period crossover (ABB/BAA)",
           "log PK response", "A", "B", rho_from_long(be_long, "A", "B"))

# ---------------------------------------------------------------------------
# (C) Split-body trials
# ---------------------------------------------------------------------------

## C1. Diabetic Retinopathy Study (Huster et al., 1989): in each patient one eye
##     was randomized to laser photocoagulation, the fellow eye left untreated.
##     Binary outcome Y = 1{blindness within tau months}; an eye contributes if
##     its status at tau is known (event before tau, or followed beyond tau).
##     Censoring is at the patient level (both eyes share censoring times).
drs <- load_cran_data("survival", "diabetic.rda")
drs_binary <- function(tau) {
  drs %>%
    mutate(y = case_when(status == 1 & time <= tau ~ 1,
                         time > tau ~ 0,
                         TRUE ~ NA_real_),
           trt = ifelse(trt == 1, "laser", "untreated")) %>%
    select(unit = id, trt, y)
}
add_result("Diabetic Retinopathy Study", "Split-body", "paired eyes",
           "blindness within 24 months", "untreated", "laser",
           rho_from_long(drs_binary(24), "untreated", "laser"))
drs_sensitivity <- bind_rows(lapply(c(12, 24, 36, 48), function(tau)
  bind_cols(tibble(horizon_months = tau),
            rho_from_long(drs_binary(tau), "untreated", "laser"))))

## C2. Box, Hunter & Hunter (1978): shoe-sole materials A and B randomly
##     allocated to the left/right foot of 10 boys; amount of wear.
sh <- load_cran_data("MASS", "shoes.rda")
add_result("Box, Hunter & Hunter (1978)", "Split-body", "paired feet",
           "sole wear", "A", "B",
           rho_from_long(wide_to_long(1:10, sh$A, sh$B, "A", "B"), "A", "B"))

# ---------------------------------------------------------------------------
# (D) Animal switchback trial
# ---------------------------------------------------------------------------

## D1. Brandt (1938): two grain mixtures fed to 10 dairy cows in a
##     T1-T2-T1 / T2-T1-T2 switchback design; milk yield (replicates averaged).
bs <- load_cran_data("agridat", "brandt.switchback.txt")
bs_long <- tibble(unit = bs$cow, trt = bs$trt, y = bs$yield, period = bs$period)
add_result("Brandt (1938)", "Animal", "3-period switchback",
           "milk yield", "T1", "T2", rho_from_long(bs_long, "T1", "T2"))


# ---------------------------------------------------------------------------
# (E) Within-subject behavioural / educational / exercise experiments
# ---------------------------------------------------------------------------

## E1. Chard et al. (2019): supplementary drinking water vs normal school day
##     in rural Mali; letter-cancellation score (visual attention).
wa <- load_cran_data("CrossCarry", "Water.rda")
wa_long <- tibble(unit = wa$ID, trt = ifelse(wa$Treatment == 1, "water", "control"),
                  y = wa$LCS, period = wa$Period)
add_result("Chard et al. (2019)", "Behavioural", "2-period crossover",
           "letter-cancellation score", "control", "water", rho_from_long(wa_long, "control", "water"))

## E2. Varao-Sousa & Kingstone (2015): live vs pre-recorded video lecture;
##     proportion correct on the memory test (format and session coincide).
ll <- load_cran_data("openintro", "lecture_learning.rda")
ll_long <- tibble(unit = ll$student, trt = ll$method, y = ll$memory)
add_result("Varao-Sousa & Kingstone (2015)", "Behavioural", "2-session within-subject",
           "memory test score", "Live", "Video", rho_from_long(ll_long, "Live", "Video"))

## E3. Lin et al. (2020): high-demand ('deplete') vs low-demand ('control')
##     initial task, followed by a Stroop task; mean reaction time over all
##     non-missing Stroop trials (the aggregation used in the afex vignettes).
st <- load_cran_data("afex", "stroop.rda")
st_long <- st %>% filter(!is.na(rt)) %>%
  transmute(unit = as.character(pno), trt = as.character(condition), y = rt)
add_result("Lin et al. (2020)", "Behavioural", "within-subject (4 studies)",
           "Stroop reaction time", "control", "deplete", rho_from_long(st_long, "control", "deplete"))
stroop_by_study <- bind_rows(lapply(levels(st$study), function(s)
  bind_cols(tibble(study = s),
            rho_from_long(st_long %>% filter(unit %in% as.character(st$pno[st$study == s])),
                          "control", "deplete"))))

## E4. Steele (1998): standing on foam vs normal surface; CTSIB balance score
##     (1-4) averaged over the three vision conditions and two replicates.
ct <- load_cran_data("faraway", "ctsib.rda")
ct_long <- tibble(unit = ct$Subject, trt = as.character(ct$Surface), y = ct$CTSIB)
add_result("Steele (1998)", "Behavioural", "within-subject",
           "balance score (CTSIB)", "norm", "foam", rho_from_long(ct_long, "norm", "foam"))

## E5. de Lucas et al. (2000): 1500 m swim with vs without wetsuit, random order;
##     maximal velocity (m/s).
ws <- load_cran_data("Lock5Data", "Wetsuits.rda")
add_result("de Lucas et al. (2000)", "Behavioural", "within-subject, random order",
           "swim velocity", "no wetsuit", "wetsuit",
           rho_from_long(wide_to_long(seq_len(nrow(ws)), ws$NoWetsuit, ws$Wetsuit,
                                      "no wetsuit", "wetsuit"), "no wetsuit", "wetsuit"))

## E6. Chipman, reported in Lawson (2015, Ch. 9): Williams crossover design with
##     three sprint conditions (1, 2, 3); sprint time. No designated control,
##     so condition 1 serves as the reference.
cs <- load_cran_data("daewr", "chipman.rda")
cs_long <- tibble(unit = cs$Subject, trt = as.character(cs$Treat), y = cs$Time, period = cs$Period)
for (k in c("2", "3"))
  add_result("Lawson (2015), sprint", "Behavioural", "3-period Williams crossover",
             "sprint time", "1", k, rho_from_long(cs_long, "1", k))

## E7. Wetherill (1972), via Preece (1982): blink rate while steering a pencil
##     along a straight vs an oscillating track.
bl <- load_cran_data("PairedData", "Blink.rda")
add_result("Wetherill (1972)", "Behavioural", "within-subject",
           "blink rate", "Straight", "Oscillating",
           rho_from_long(wide_to_long(bl$Subject, bl$Straight, bl$Oscillating,
                                      "Straight", "Oscillating"), "Straight", "Oscillating"))

# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

rho_tab_additional <- bind_rows(results)

options(width = 200, pillar.sigfig = 3)
print(rho_tab_additional %>% select(study, contrast, n, rho, ci_lower, ci_upper, spearman, rho_period_adj),
      n = Inf)

cat("\nNumber of dataset entries:", n_distinct(rho_tab_additional$study),
    " | number of treatment contrasts:", nrow(rho_tab_additional), "\n")
cat("Contrasts with rho > 0   :", sum(rho_tab_additional$rho > 0), "\n")
cat("Contrasts with rho > 0.33:", sum(rho_tab_additional$rho > 0.33), "\n")
cat("Contrasts with rho > 0.5 :", sum(rho_tab_additional$rho > 0.5), "\n")
cat("Median rho:", round(median(rho_tab_additional$rho), 3),
    " | mean rho:", round(mean(rho_tab_additional$rho), 3), "\n")
cat("Contrasts whose approximate Fisher-z 95% CI excludes 0:",
    sum(rho_tab_additional$ci_lower > 0 | rho_tab_additional$ci_upper < 0), "\n")

cat("\nDiabetic Retinopathy Study, sensitivity to the horizon:\n")
print(drs_sensitivity)
cat("\nLin et al. (2020), separately by study:\n")
print(stroop_by_study)

# write.csv(rho_tab_additional, "rho_additional_datasets.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Screened but excluded candidates (reason in brackets)
# ---------------------------------------------------------------------------
# Fewer than 10 complete pairs: agridat::lucas.switchback and
#   agridat::patterson.switchback [each cow receives only 2 of the 3/4
#   treatments, <= 4 pairs per contrast]; nlme::Tetracycline1/2 [5 subjects];
#   daewr::bioeqv [3]; PairedData::Tobacco and ::Sewage [8]; BSDA::Coffee [9];
#   nlme::ergoStool [9]; nlme::Machines [6]; nlme::Meat [<= 3 per contrast].
# Not a within-unit contrast: agridat::depalluel.sheep [one diet per animal];
#   exactRankTests::ears [treatment assigned per child, both ears same arm];
#   heplots::Headache, Lock5Data::CaffeineTaps [between-subject].
# Simulated / hypothetical: BE::NCAResult4BE, Lock5Data::TrafficFlow,
#   heplots::ReactTime and afex::md_12.1 [Maxwell & Delaney hypothetical data].
# Duplicate of a replicateBE dataset: BElikelihood::dat [same Cmax values as DS26].
# Only within-subject differences available: Stat2Data::Eyes, datasets::sleep.
# Matched pairs of distinct units (not the same unit): Darwin's maize
#   (HistData::ZeaMays, PairedData::Corn, DAAG::mignonette), sibling chicks
#   (PairedData::ChickWeight), PairedData::BloodLead, paired plots
#   (PairedData::Barley, ::Grain, ::Grain2), Sleuth3::case0202 [twins].
# Method comparison / no intervention: PairedData::Iron, ::Meat, PASWR2::TENNIS,
#   Lock5Data::QuizPulse10, PairedData::GrapeFruit, PairedData::Anorexia
#   [before-after], product-rating designs (Lock5Data::WaterTaste,
#   crossdes::biertest.dat, Lock5Data::StorySpoilers).
