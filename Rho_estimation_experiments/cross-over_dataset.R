# Cross-over study where each subject receives both treatments (R and T) in different periods. 
# We simply compute the direct paired empirical correlation cor(Y1, Y0)

library(dplyr); 
library(tidyr); 
library(readr)

url <- "https://raw.githubusercontent.com/Helmut01/replicateBE/master/inst/extdata/DS01.csv"

read_csv(url, skip = 2, show_col_types = FALSE) %>%
  mutate(subject = as.character(subject), treatment = as.character(treatment), logPK = as.numeric(logPK)) %>%
  filter(treatment %in% c("R", "T"), !is.na(logPK)) %>%
  group_by(subject, treatment) %>%
  summarise(Y = mean(logPK), .groups = "drop") %>%
  pivot_wider(names_from = treatment, values_from = Y) %>%
  rename(Y0 = R, Y1 = T) %>%
  filter(!is.na(Y0), !is.na(Y1)) %>%
  summarise(
    rho = cor(Y1, Y0),
    ci_lower = tanh(atanh(rho) - qnorm(0.975) / sqrt(n() - 3)),
    ci_upper = tanh(atanh(rho) + qnorm(0.975) / sqrt(n() - 3))
  )




# replicate for all datasets DS01-DS30.
library(dplyr)
library(tidyr)
library(readr)

estimate_rho_ds <- function(ds_id) {
  url <- sprintf(
    "https://raw.githubusercontent.com/Helmut01/replicateBE/master/inst/extdata/DS%02d.csv",
    ds_id
  )
  
  tryCatch({
    raw_dat <- read_csv(
      url,
      comment = "#",           # ignore leading comment lines automatically
      na = c("", "."),
      show_col_types = FALSE
    )
    
    if (!all(c("subject", "treatment") %in% names(raw_dat))) {
      stop("required columns not found")
    }
    
    paired_dat <- raw_dat %>%
      mutate(
        subject   = as.character(subject),
        treatment = as.character(treatment),
        Y = if ("logPK" %in% names(.)) {
          as.numeric(logPK)
        } else if ("PK" %in% names(.)) {
          log(as.numeric(PK))
        } else {
          stop("neither logPK nor PK found")
        }
      ) %>%
      filter(treatment %in% c("R", "T"), !is.na(Y)) %>%
      group_by(subject, treatment) %>%
      summarise(Y = mean(Y), .groups = "drop") %>%
      pivot_wider(names_from = treatment, values_from = Y) %>%
      rename(Y0 = R, Y1 = T) %>%
      filter(!is.na(Y0), !is.na(Y1))
    
    n_pairs <- nrow(paired_dat)
    rho <- if (n_pairs >= 2) cor(paired_dat$Y1, paired_dat$Y0) else NA_real_
    
    if (is.finite(rho) && n_pairs > 3 && abs(rho) < 1) {
      se <- 1 / sqrt(n_pairs - 3)
      z  <- atanh(rho)
      ci_lower <- tanh(z - qnorm(0.975) * se)
      ci_upper <- tanh(z + qnorm(0.975) * se)
    } else {
      ci_lower <- NA_real_
      ci_upper <- NA_real_
    }
    
    tibble(
      dataset  = sprintf("DS%02d", ds_id),
      n_pairs  = n_pairs,
      rho      = rho,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      error    = NA_character_
    )
    
  }, error = function(e) {
    tibble(
      dataset  = sprintf("DS%02d", ds_id),
      n_pairs  = NA_integer_,
      rho      = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      error    = conditionMessage(e)
    )
  })
}

rho_tab <- bind_rows(lapply(1:30, estimate_rho_ds)) %>%
  arrange(desc(rho))

rho_tab
