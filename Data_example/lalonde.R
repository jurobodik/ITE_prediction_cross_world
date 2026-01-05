# library(causaldata)
# data(nsw_mixtape)
# 
# 
# library(MatchIt)
# data("lalonde")


# Lalonde Sample (no re74)

colnames_nsw <- c(
  "treat",     # treatment indicator
  "age",
  "educ",
  "is_black",
  "is_hispanic",
  "married",
  "nodegree",
  "re75",      # earnings 1975 (pre-treatment)
  "re78"       # earnings 1978 (outcome)
)

nsw_control <- read.table(
  "/Users/yaxuan/PycharmProjects/ITE_prediction_cross_world/Data_example/nsw_control.txt",
  header = FALSE
)

nsw_treated <- read.table(
  "/Users/yaxuan/PycharmProjects/ITE_prediction_cross_world/Data_example/nsw_treated.txt",
  header = FALSE
)

colnames(nsw_control) <- colnames_nsw
colnames(nsw_treated) <- ccolnames_nsw
nsw_full <- rbind(nsw_control, nsw_treated)


# Lalonde Sample (no re74)

colnames_nswre74 <- c(
  "treat",     # treatment indicator
  "age",
  "educ",
  "is_black",
  "is_hispanic",
  "married",
  "nodegree",
  "re74",
  "re75",      # earnings 1975 (pre-treatment)
  "re78"       # earnings 1978 (outcome)
)

nswre74_control <- read.table(
  "/Users/yaxuan/PycharmProjects/ITE_prediction_cross_world/Data_example/nswre74_control.txt",
  header = FALSE
)

nswre74_treated <- read.table(
  "/Users/yaxuan/PycharmProjects/ITE_prediction_cross_world/Data_example/nswre74_treated.txt",
  header = FALSE
)

colnames(nswre74_control) <- colnames_nswre74
colnames(nswre74_treated) <- colnames_nswre74
nswre74_full <- rbind(nswre74_control, nswre74_treated)





