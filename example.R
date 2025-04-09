"
This script serves as a demo - example to try out and possibly change for your dataset
First, we generate synthetic data using the `data_synthetic()` function. It simulates realistic treatment-outcome dependencies using various copulas and marginal noise distributions.
Then, we apply the `D_rho_intervals` function to estimate prediction intervals for Individual Treatment Effects (ITE).
The `D_rho_intervals` function can be found in the `D_rho_intervals_function.R` file

Function source: source(D_rho_intervals_function.R)

This function is applied on the synthetic dataset and plotted the results below
"

library(MASS)
library(akima)
library(rgl)
library(raster)
library(ambient)
library(dplyr)      
library(plotly)
library(ggplot2)
library(gridExtra)
library(forcats)
library(patchwork)
library(purrr)
library(grid)
library(copula)

####################################################################################################
################# Generating benchmark datasets ####################################################
####################################################################################################

data_synthetic <- function(n = 1000, 
                           d = 2, 
                           rho, 
                           sigma_1 = 1, 
                           sigma_2 = 4, 
                           constant_propensity = FALSE, 
                           copula_type = "gaussian", #choices: "gaussian", "clayton", "t", "gumbel"
                           marginal = "gaussian"){ #choices: "gaussian", "laplace", "t", "laplace", "chisq"
  # generate random 1D and 2D functions
  random_function_1d <- function(freq = 0.1) {
    s = seq(-10, 10, length.out = 1001)
    f <- long_grid(s)
    f$noise <- rep(0, length(s))
    while (all(f$noise[300:600] == 0)) {
      f$noise <- gen_perlin(f$x, frequency = freq, fractal = 'rigid-multi')
    }
    amplitude = 10 / (max(f$noise) + 1)
    return(f$noise * amplitude)
  }
  
  evaluation_of_f_1d <- function(f, X1) {
    minim = min(X1)
    maxim = max(X1)
    sapply(X1, function(x) f[round((x - minim) * 1000 / (maxim - minim) + 1)])
  }
  
  random_function_2d <- function(freq = 0.001) {
    f = noise_perlin(c(1001, 1001), frequency = freq, fractal = 'fbm', octaves = 2, lacunarity = 2, gain = 0.4)
    f = f^2
    amplitude = 10 / max(f)
    return(f * amplitude)
  }
  
  evaluation_of_f_2d <- function(f, X1, X2) {
    minim_x = min(X1)
    maxim_x = max(X1)
    minim_y = min(X2)
    maxim_y = max(X2)
    
    mapply(function(x, y) {
      xx = round((x - minim_x) * 1000 / (maxim_x - minim_x) + 1)
      yy = round((y - minim_y) * 1000 / (maxim_y - minim_y) + 1)
      f[xx, yy]
    }, X1, X2)
  }
  
  # Covariates and CATE
  if (d == 1) {
    X = runif(n, -1, 1)
    CATE = evaluation_of_f_1d(random_function_1d(), X)
    mu = 5 + 5 * X
  } else {
    Sigma <- matrix(0.25, nrow = d, ncol = d)
    diag(Sigma) <- rep(1, d)
    tilde_X <- mvrnorm(n, mu = rep(0, d), Sigma = Sigma)
    X <- pnorm(tilde_X)
    tau = random_function_2d()
    CATE = evaluation_of_f_2d(tau, X[,1], X[,2])
    beta = rnorm(d)
    mu = X %*% beta
  }
  
  # Error generation using copulas
  generate_errors <- function(n, rho, sigma_1 , sigma_2, copula_type = "gaussian", marginal = "gaussian") {
    qlaplace <- function(p) {ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))}
    
    tau = (2 / pi) * asin(rho); #relation between correlation and parameter in clayton/gumbel copula
    if(tau <= -1) tau = -0.99; if(tau >= 1) tau = 0.99 #soft fix for edge cases
    cop <- switch(copula_type,
                  "gaussian" = normalCopula(param = rho, dim = 2),  # directly use rho
                  "t"        = tCopula(param = rho, dim = 2, df = 3),  # also takes rho directly
                  "clayton"  = claytonCopula(param = iTau(claytonCopula(), tau), dim = 2), #parameter and rho have this relation: (2 / pi) * asin(rho), althoguh rho is kendels tau not pearson
                  "gumbel"   = gumbelCopula(param = iTau(gumbelCopula(), tau), dim = 2),   #parameter and rho have this relation: (2 / pi) * asin(rho)
                  stop("Unsupported copula type.")
    )
    
    u <- rCopula(n, cop)
    
    transform_marginal <- function(u_vec, type) {
      switch(type,
             "gaussian" = qnorm(u_vec),
             "t" = qt(u_vec, df = 3),
             "laplace" = qlaplace(u_vec),
             "chisq" = qchisq(u_vec, df = 3),
             stop("Unsupported marginal type.")
      )
    }
    
    eps1 <- transform_marginal(u[, 1], marginal) * sqrt(sigma_1)
    eps2 <- transform_marginal(u[, 2], marginal) * sqrt(sigma_2)
    
    data.frame(eps1, eps2)
  }
  
  eps <- generate_errors(n, rho, sigma_1, sigma_2, copula_type, marginal)
  epsilon1 <- eps[,1]
  epsilon2 <- eps[,2]
  
  Y0 <- mu + epsilon1
  Y1 <- mu + CATE + epsilon2
  
  if (!constant_propensity) {
    propensity_score = if (d == 1) (1 + abs(X)) / 4 else (1 + abs(X[,1])) / 4
    treatment = rbinom(n, 1, 1 - propensity_score)
  } else {
    treatment = sample(c(0, 1), n, replace = TRUE)
  }
  
  Y_obs = ifelse(treatment == 1, Y1, Y0)
  
  if (d == 1) {
    return(data.frame(X, Y0, Y1, Y_obs, treatment))
  } else {
    return(data.frame(X, Y0 = Y0, Y1 = Y1, Y_obs = Y_obs, treatment = treatment))
  }
}





rho_true = 0.1
rho_used = rho_true

d=1
n=1000
data_all = data_synthetic(n = 2000, d = d, rho = rho_used,
                          copula_type = 'gaussian',
                          marginal = 'gaussian',
                          constant_propensity = FALSE)
data = data_all[1:n,]
data_test = data_all[(n+1):nrow(data_all),]
new_points = data.frame(data_test[,1:d]); names(new_points) = paste0('X', 1:d)

X = data[, 1:d]
Y0 = data$Y0
Y1 = data$Y1
Y = data$Y_obs
treatment = data$treatment




D_rho = D_rho_intervals(X, Y, treatment, new_points, rho = rho_used, 
                                weighted_conformal = TRUE, #should we estimate propensity scores and use weighted conformal prediction?
                                conformal = 'CQR', #choices: 'CQR', 'PCS', 'CLEAR', 'naive' 
                                add_confidence_intervals = TRUE, #If TRUE, then we compute confidence intervals for CATE and add them to the prediction intervals
                                number_of_bootstraps = 20)




true_ITE <- data_test$Y1 - data_test$Y0
coverage_ITE <- mean(true_ITE >= D_rho$lower & true_ITE <= D_rho$upper)
length_ITE <- mean(D_rho$upper - D_rho$lower)

cat(" Coverage for ITE: ", coverage_ITE, "\n", '  Average Length: ',length_ITE, "\n")


############# Plotting test dataset and our prediction intervals #############
####Its long because we also show CQR prediction intervals for the treated and untreated units separatelly

if(d==1){
  
  library(ggplot2)
  
  # Add label column for legend mapping
  df_plot$point_type <- "Observed"
  df_true <- df_plot
  df_true$point_type <- "True ITE"
  
  # Combine observed and true ITE data for unified plotting
  df_all <- rbind(df_plot, df_true)
  
  ggplot(df_all, aes(x = X)) +
    # True ITE points (green)
    geom_point(data = subset(df_all, point_type == "True ITE"),
               aes(y = true_ITE, color = point_type), alpha = 0.5, size = 1.2) +
    # Observed outcomes colored by treatment
    geom_point(data = subset(df_all, point_type == "Observed"),
               aes(y = Y_obs, color = treatment), alpha = 0.5, size = 1.5) +
    # Estimated CATE line
    geom_line(data = df_plot, aes(y = CATE), color = "black", size = 1) +
    # Interval ribbon
    geom_ribbon(data = df_plot, aes(ymin = lower, ymax = upper), fill = "green", alpha = 0.3) +
    # Reference line at 0
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    # Labels and theme
    labs(
      x = "X",
      y = "Y_obs / CATE / true ITE",
      color = "Legend",
      title = bquote("ITE prediction intervals using " * rho == .(rho_used) *
                       " | Coverage: " * .(round(coverage_ITE, 3)) *
                       " | Avg. Length: " * .(round(length_ITE, 3)))
    ) +
    scale_color_manual(
      values = c("0" = "blue", "1" = "red", "True ITE" = "green"),
      labels = c("Control (t=0)", "Treated (t=1)", "True ITE")
    ) +
    theme_minimal()


  
  
  
  
  
  
  
  
  
  # Prepare treated and untreated training sets
  X_treated <- X[treatment == 1]
  Y_treated <- Y[treatment == 1]
  
  X_untreated <- X[treatment == 0]
  Y_untreated <- Y[treatment == 0]
  
  # Run CQR separately
  CQR_treated <- CQR(X_treated, Y_treated, new_points,
                     desired_coverage = 0.9,
                     CQR_qr = 'qGAM',  # or 'auto'
                     CQR_type = 'multiplicative')
  
  CQR_untreated <- CQR(X_untreated, Y_untreated, new_points,
                       desired_coverage = 0.9,
                       CQR_qr = 'qGAM',  # or 'auto'
                       CQR_type = 'multiplicative')
  
  
  
  # Prepare data
  df_combined <- data.frame(
    X = D_rho$new_points[,1],
    CATE = as.numeric(D_rho$CATE),
    lower_CATE = as.numeric(D_rho$lower),
    upper_CATE = as.numeric(D_rho$upper),
    Y_obs = data_test$Y_obs,
    treatment = factor(data_test$treatment),
    true_ITE = data_test$Y1 - data_test$Y0,
    pred_Treated = CQR_treated$hat_f,
    lower_Treated = CQR_treated$lower,
    upper_Treated = CQR_treated$upper,
    pred_Untreated = CQR_untreated$hat_f,
    lower_Untreated = CQR_untreated$lower,
    upper_Untreated = CQR_untreated$upper
  )
  
  # Sort by X
  df_combined <- df_combined[order(df_combined$X), ]
  
  # Plot
  ggplot(df_combined, aes(x = X)) +
    # True ITE points
    geom_point(aes(y = true_ITE, color = "True ITE"), alpha = 0.5, size = 1.2) +
    
    # Observed outcomes
    geom_point(aes(y = Y_obs, color = treatment), alpha = 0.5, size = 1.5) +
    
    # CATE line and interval
    geom_line(aes(y = CATE), color = "black", size = 1) +
    geom_ribbon(aes(ymin = lower_CATE, ymax = upper_CATE, fill = "CATE Interval"), alpha = 0.3) +
    
    # Treated and untreated intervals
    geom_ribbon(aes(ymin = lower_Treated, ymax = upper_Treated, fill = "Treated Interval"), alpha = 0.2) +
    geom_line(aes(y = pred_Treated), color = "red", linetype = "solid", size = 0.8) +
    
    geom_ribbon(aes(ymin = lower_Untreated, ymax = upper_Untreated, fill = "Control Interval"), alpha = 0.2) +
    geom_line(aes(y = pred_Untreated), color = "blue", linetype = "solid", size = 0.8) +
    
    # Labels and scales
    labs(
      x = "X",
      y = "Outcome / CATE / ITE",
      color = "Points",
      fill = "Prediction Intervals",
      title = bquote("ITE and Potential Outcome Prediction Intervals (ρ = " * .(rho_used) * ")"),
      subtitle = bquote("Coverage on test set = " * .(round(coverage_ITE, 3)) *
                          " | Avg. Length = " * .(round(length_ITE, 3)))
    ) +
    scale_color_manual(
      values = c("0" = "blue", "1" = "red", "True ITE" = "green"),
      labels = c("Control (t=0)", "Treated (t=1)", "True ITE")
    ) +
    scale_fill_manual(
      values = c("CATE Interval" = "green", 
                 "Treated Interval" = "red", 
                 "Control Interval" = "blue"),
      guide = guide_legend(override.aes = list(alpha = 0.3))
    ) +
    theme_minimal()
  

}
