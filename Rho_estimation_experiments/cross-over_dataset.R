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
