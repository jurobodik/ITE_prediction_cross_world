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
colnames(nsw_treated) <- colnames_nsw
nsw_full <- rbind(nsw_control, nsw_treated)


# Lalonde Sample (with re74)

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


nswre74_full$treat <- as.integer(nswre74_full$treat)
xvars <- c("age","educ","is_black","is_hispanic","married","nodegree","re74","re75")


# model fitting 

fml_y <- as.formula(paste("re78 ~", paste(xvars, collapse = " + ")))

lm1 <- lm(fml_y, data = nswre74_treated)
lm0 <- lm(fml_y, data = nswre74_control)

nswre74_full$mu1_hat <- predict(lm1, newdata = nswre74_full)
nswre74_full$mu0_hat <- predict(lm0, newdata = nswre74_full)


# matching

library(MatchIt)

fml_x <- as.formula(paste("treat ~", paste(xvars, collapse = " + ")))
m.out <- matchit(
  formula  = fml_x,
  data     = nswre74_full,
  method   = "nearest",
  distance = "mahalanobis", # covariate matching
  replace  = TRUE,
  ratio    = 1
)

print(summary(m.out))
mm <- m.out$match.matrix
treated_ids <- rownames(mm)
control_ids <- as.character(mm[, 1])

# library(dplyr)

res1 <- nswre74_full[treated_ids, "re78"] - nswre74_full[treated_ids, "mu1_hat"]
res0 <- nswre74_full[control_ids, "re78"] - nswre74_full[control_ids, "mu0_hat"]

rho_proxy <- cor(res1, res0, use = "complete.obs")
cat("Number of matched pairs:", length(res1), "\n")
cat("Proxy residual correlation (Mahalanobis matched pairs):", round(rho_proxy, 3), "\n")



