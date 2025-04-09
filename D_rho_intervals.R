################################################################################
# D_rho_intervals.R
#
# Title: ITE Prediction Intervals under Cross-World Assumptions 
# Author: jurajbodik.com
# Description:
#   This script provides an implementation of prediction intervals 
#   for Individual Treatment Effects (ITE) using D_rho_intervals under
#   the Cross-World Assumption (correlation parameter `rho`). 
#   It supports several conformal prediction methods (CQR, PCS, CLEAR, Naive) 
#   and two CATE estimation strategies (T-learner and Causal Forest).
#
#   The core function is `D_rho_intervals`, which:
#     - Estimates ITE prediction intervals at new covariate values
#     - Combines uncertainty from both treated and control outcomes
#     - Optionally adds bootstrap-based confidence intervals (recommended, but a bit slower)
#
#   Dependencies: randomForest, quantregForest, qgam, mgcv, locfit, gbm, matrixStats, grf, Hmisc, MASS
#   Input:
#     - X: covariates (data.frame)
#     - Y: observed outcomes (numeric vector)
#     - treatment: treatment indicator (0/1)
#     - new_points: covariates where prediction is needed (data.frame)
#     - rho: assumed correlation between potential outcomes (number)
#
#   Output: Prediction intervals for ITE (CATE ± uncertainty), optionally with added CI (recommended, but a bit slower)

################################################################################
# Example Usage:
#
# ################# Simulate synthetic data
# set.seed(0)
#  n <- 1000
#  d <- 5
#  rho=0.5
#  
#  X <- as.data.frame(matrix(rnorm(n * d), ncol = d))
#  treatment <- rbinom(n, 1, 0.5)
#  epsilons <- mvrnorm(n = n, mu = c(0,0), Sigma = matrix(c(1, rho, rho, 1), nrow = 2)) 
#  Y0 <- X[,1]              + epsilons[,1]
#  Y1 <- X[,1] + 1+ X[,2]^2 + epsilons[,2]
#  Y <- ifelse(treatment == 1, Y1, Y0)
# 
# ################# Define new points to predict ITEs
#  new_points <- as.data.frame(matrix(rnorm(10 * d), ncol = d))
# #################  Compute D_rho_intervals
#  result <- D_rho_intervals(X = X,
#                            Y = Y,
#                            treatment = treatment,
#                            new_points = new_points,
#                            rho = rho,
#                            conformal = "CQR", 
#                            weighted_conformal = FALSE,
#                            add_confidence_intervals = FALSE)
#
# print(result$CATE)       # Predicted CATEs
# print(result$lower)      # Lower prediction bounds
# print(result$upper)      # Upper prediction bounds
###print(result$CATE - 1-new_points[,2]^2) #print error: hat{CATE} - CATE 
###############################################################################


library(randomForest)
library(quantregForest)
library(qgam)
library(mgcv)
library(locfit)
library(gbm)
library(matrixStats)
library(grf)
library(Hmisc)
library(MASS)


D_rho_intervals  <- function(X, Y, treatment, new_points,
                        rho,
                        desired_coverage = 0.9, 
                        train_calib_split = 0.8, 
                        weighted_conformal = FALSE, #should we estimate propensity scores and use weighted conformal prediction?
                        conformal = 'CQR', #choices: 'CQR', 'PCS', 'CLEAR', 'naive' 
                        CQR_qr = 'auto', #choices: 'RF', 'qGAM', 'naive', 'auto' where 'auto = RF' if d>5 and 'auto = qGAM' if d<=5
                        CQR_type = 'multiplicative', #choices: 'multiplicative', 'additive', 'naive'
                        CATE_estimate = 'T-learner', #choices: 'T-learner', 'Causal_forest'
                        add_confidence_intervals = TRUE, #If TRUE, then we compute confidence intervals for CATE via bootstrap and add them to the prediction intervals
                        number_of_bootstraps = 50){ #Number of bootstraps for confidence intervals when add_confidence_intervals==TRUE
  
  
  
  ####################################################################################
  ####################### Some helpful functions #####################################
  ####################################################################################
  rho_distance = function(x, y, rho)return(sqrt(x^2 + y^2 - 2*rho*x*y))
  CI_addition = function(x, rho) return( x*((1+rho)/2)^2 )                      
  
  
  
  
  ####################################################################################
  ####################### Some data handling and splitting ###########################
  ####################################################################################
  data_original = data.frame(X=X, Y=Y, treatment = treatment)
  X=data.frame(X); Y=data.frame(Y);
  
  d = ncol(X); if (is.null(nrow(X))) {d=1}
  names(data_original) <- c(paste0("X", 1:d), "Y", 'treatment')
  names(X) <- c(paste0("X", 1:d))
  new_points = data.frame(new_points); names(new_points) <- c(paste0("X", 1:d))
  my_sequence = rbind(X, new_points)
  
  data = data_original

  ####################################################################################
  ############# Lets clean the hyper-parameters into a nice format ###################
  ####################################################################################
  if(conformal=='Naive')conformal = 'naive'
  if(CQR_qr == 'Naive' )CQR_qr = 'naive'
  if(CQR_qr == 'auto'){if(ncol(data.frame(X))>5){CQR_qr = 'RF'}else{CQR_qr = 'qGAM'}}
  if(CQR_qr != 'RF' & CQR_qr != 'qGAM' & CQR_qr != 'naive')stop("Quantile regression method not recognized. Choose 'RF' or 'qGAM' or 'naive'")

  
  
  ####################################################################################
  ############# Now we upload CQR and PCS and CLEAR uncertainty estimates ############
  ####################################################################################
  
  
  
  
  
  
  CQR_weighted <- function(X, Y, new_points, treatment, X_all,
                           our_weights = '1/e',  # or '1/(1-e)'
                           desired_coverage = 0.9,
                           train_calib_split = 0.8,
                           CQR_qr = 'auto',
                           CQR_type = 'multiplicative', 
                           propensity_clipping = 0.2) {
    
    # ------------------ Setup and checks ------------------
    if (CQR_qr == 'Naive') CQR_qr <- 'naive'
    if (CQR_qr == 'auto') CQR_qr <- if (ncol(data.frame(X)) > 5) 'RF' else 'qGAM'
    if (!CQR_qr %in% c('RF', 'qGAM', 'naive')) stop("Unsupported CQR_qr")
    
    treatment <- as.numeric(treatment)
    d <- ncol(data.frame(X))
    X <- data.frame(X)
    Y <- as.numeric(unlist(Y))
    new_points <- data.frame(new_points)
    names(X) <- names(new_points) <- paste0("X", 1:d)
    
    # ------------------ Split train/calibration ------------------
    n <- nrow(X)
    train_idx <- 1:(n * train_calib_split)
    calib_idx <- (n * train_calib_split + 1):n
    
    X_train <- data.frame(X[train_idx, ]); names(X_train) = names(X)
    Y_train <- Y[train_idx]
    X_calib <- data.frame(X[calib_idx, ]); names(X_calib) = names(X)
    Y_calib <- Y[calib_idx]
    treatment_train <- treatment[train_idx]
    treatment_calib <- treatment[calib_idx]
    X_all <- data.frame(X_all); names(X_all) = names(X)
    
    all_X <- rbind(X, new_points)
    
    # ------------------ Fit quantile regression ------------------
    create_formula <- function(X, smoothing = TRUE) {
      terms <- sapply(1:ncol(X), function(j) {
        if (length(unique(X[, j])) > 9 && smoothing) paste0("s(X", j, ")")
        else paste0("as.factor(X", j, ")")
      })
      reformulate(terms, response = "Y")
    }
    
    if (CQR_qr == 'qGAM') {
      formula <- create_formula(X_train)
      invisible(capture.output({
        fit_95 <- qgam(formula, data = data.frame(X_train, Y = Y_train), qu = (1 + desired_coverage) / 2)
        fit_50 <- qgam(formula, data = data.frame(X_train, Y = Y_train), qu = 0.5)
        fit_05 <- qgam(formula, data = data.frame(X_train, Y = Y_train), qu = (1 - desired_coverage) / 2)
      }))
      pred_95 <- predict(fit_95, newdata = all_X)
      pred_50 <- predict(fit_50, newdata = all_X)
      pred_05 <- predict(fit_05, newdata = all_X)
    } else {
      qrf <- quantregForest(x = X_train, y = Y_train, ntree = 1000, nodesize = 5)
      pred_95 <- predict(qrf, newdata = all_X, what = (1 + desired_coverage) / 2)
      pred_50 <- predict(qrf, newdata = all_X, what = 0.5)
      pred_05 <- predict(qrf, newdata = all_X, what = (1 - desired_coverage) / 2)
    }
    
    # ------------------ Enforce order of quantiles ------------------
    pred_95 <- pmax(pred_95, pred_50)
    pred_05 <- pmin(pred_05, pred_50)
    
    # ------------------ Propensity scores ------------------
    data_ps <- data.frame(treatment_train = treatment_train, X_all[train_idx, ]); names(data_ps) <- c("treatment_train", names(X_all))
    
    ps_model <- gbm(treatment_train ~ ., data = data_ps,
                    distribution = "bernoulli", n.trees = 300, interaction.depth = 3,
                    shrinkage = 0.05, verbose = FALSE)
    newdata_ps <- data.frame(X_all[calib_idx, ]);    colnames(newdata_ps) <- names(data_ps)[-1]  # Just to be safe
    e_hat <- predict(ps_model, newdata = newdata_ps, n.trees = 300, type = "response")
    e_hat <- pmin(pmax(e_hat, propensity_clipping), 1-propensity_clipping)
    
    # ------------------ Calibration ------------------
    estimate_gamma <- function(type) {
      if (type == "additive") {
        scores <- mapply(function(i) {
          idx <- calib_idx[i]
          max(pred_05[idx] - Y_calib[i], Y_calib[i] - pred_95[idx])
        }, seq_along(Y_calib))
      } else if (type == "multiplicative") {
        l <- pred_50 - pred_05
        u <- pred_95 - pred_50
        scores <- mapply(function(i) {
          idx <- calib_idx[i]
          max((pred_50[idx] - Y_calib[i]) / l[idx], (Y_calib[i] - pred_50[idx]) / u[idx])
        }, seq_along(Y_calib))
      } else stop("Unknown CQR_type")
      
      weights <- if (our_weights == '1/e') 1 / e_hat else 1 / (1 - e_hat)
      wtd.quantile(scores, weights = weights, probs = desired_coverage)
    }
    
    gamma <- estimate_gamma(CQR_type)
    
    # ------------------ Final intervals ------------------
    if (CQR_type == "multiplicative") {
      lower <- pred_50 - gamma * (pred_50 - pred_05)
      upper <- pred_50 + gamma * (pred_95 - pred_50)
    } else {
      lower <- pred_05 - gamma
      upper <- pred_95 + gamma
    }
    
    # ------------------ Extract new points ------------------
    m <- nrow(new_points)
    lower <- tail(lower, m)
    upper <- tail(upper, m)
    prediction_f <- tail(pred_50, m)
    
    list(
      hat_f = prediction_f,
      lower = lower,
      upper = upper,
      new_points = new_points
    )
  }
  
  
  
  
  
  Naive_conformal = function(X, Y, new_points,
                             desired_coverage = 0.9, 
                             train_calib_split = 0.8, 
                             CQR_qr = 'auto'){
    
    
    
    if(CQR_qr == 'auto'){if(ncol(data.frame(X))>5){CQR_qr = 'RF'}else{CQR_qr = 'qGAM'}}
    if(CQR_qr == 'Naive' )CQR_qr = 'RF'
    if(CQR_qr != 'RF' & CQR_qr != 'qGAM' & CQR_qr != 'naive')stop("Quantile regression method not recognized. Choose 'RF' or 'qGAM' or 'naive'")
    ###########################################################
    #################Data handling and splitting###############
    ###########################################################
    
    data_original = data.frame(X=X, Y=Y)
    X=data.frame(X); Y=data.frame(Y);
    
    d = ncol(X); if (is.null(nrow(X))) {d=1}
    names(data_original) <- c(paste0("X", 1:d), "Y")
    names(X) <- c(paste0("X", 1:d))
    new_points = data.frame(new_points); names(new_points) <- c(paste0("X", 1:d))
    my_sequence = rbind(X, new_points)
    n=nrow(Y);
    train = 1:(n*train_calib_split); calib = (n*train_calib_split+1):n
    
    X_train = X[train,]; Y_train = Y[train,];
    X_calib = X[calib,]; Y_calib = Y[calib,];
    
    n_train = length(Y_train); n_calib = length(Y_calib);
    index_sequences = (n_train+1):(n_train+n_calib) #This is used to extract the new points from the combined sequence
    
    
    ###########################################################
    ################# hat_f ###################################
    ###########################################################
    
    create_formula = function(X, smoothing=TRUE){
      form="Y ~ "
      d= ncol(X); if (is.null(nrow(X))) {d=1}
      if(d==1){if(smoothing==FALSE)return(as.formula('Y~X1'))else return(as.formula('Y~s(X1)'))}
      if(smoothing==TRUE){
        for (i in 1:d) {#Discrete variable is if it contains <=9 different values
          if (  length(unique(X[,i]))>9  ) {form=paste(form, paste0("s(X", i,  ")+"))}
          if (  length(unique(X[,i]))<=9) {form=paste(form, paste0("as.factor(X", i,  ")+"))  }
        }
      }
      if(smoothing==FALSE){
        for (i in 1:d) {form=paste(form, paste0("X", i, "+"))}
      }
      
      form=substr(form,1,nchar(form)-1)
      return(as.formula(form))
    }
    
    X = X_train; Y = Y_train
    n=length(Y_train);
    data = data.frame(X=X_train, Y=Y_train);   names(data) <- c(paste0("X", 1:d), "Y")
    if(d==1){newdata = data.frame(X1 = c(my_sequence))}else{newdata = my_sequence}
    
    if(CQR_qr == 'qGAM'){
      fit_f <- gam(create_formula(X), data = data)
      hat_f = predict(fit_f, newdata = newdata)
    }
    if(CQR_qr == 'RF' || CQR_qr == 'naive'){
      rf_model <- randomForest(x = data.frame(x = X), y = Y, ntree = 1000, nodesize = 5)
      names(newdata) <- names(rf_model$forest$xlevels)
      hat_f <- predict(rf_model, newdata = data.frame(newdata))
    }
    
    
    ###########################################################
    ################# Calibration #############################
    ###########################################################
    X = X_calib; Y = Y_calib
    n=length(Y);
    
    score = c()
    for(i in 1:n){
      index = index_sequences[i]
      s = abs(hat_f[index] - Y[i])
      score = c(score, s)
    }
    gamma = quantile(score, desired_coverage)
    
    lower = hat_f - gamma
    upper = hat_f + gamma
    
    
    
    ###########################################################
    ################# Output ##################################
    ###########################################################
    n=length(data_original$Y)
    if(!is.null(new_points)){
      m=length(new_points$X1)
      hat_f = hat_f[(n+1):(n+m)]
      lower = lower[(n+1):(n+m)]
      upper = upper[(n+1):(n+m)]
    }
    
    
    return(list(hat_f = hat_f,
                lower = lower, 
                upper = upper,
                new_points = new_points))
  }
  
  
  CQR <- function(X, Y, new_points,
                  desired_coverage = 0.9,
                  train_calib_split = 0.8,
                  CQR_qr = 'auto',
                  CQR_type = 'multiplicative') {
    
    # ------------------ Setup ------------------
    if (CQR_qr == 'Naive') CQR_qr <- 'naive'
    if (CQR_qr == 'auto') CQR_qr <- if (ncol(data.frame(X)) > 5) 'RF' else 'qGAM'
    if (!CQR_qr %in% c('RF', 'qGAM', 'naive')) stop("Unsupported CQR_qr method")
    
    d <- ncol(data.frame(X))
    X <- data.frame(X); Y <- as.numeric(Y)
    new_points <- data.frame(new_points)
    names(X) <- names(new_points) <- paste0("X", 1:d)
    
    # ------------------ Split ------------------
    n <- nrow(X)
    train_idx <- 1:(n * train_calib_split)
    calib_idx <- (n * train_calib_split + 1):n
    
    X_train <- data.frame(X[train_idx, ]); names(X_train) = names(X)
    Y_train <- Y[train_idx]
    X_calib <- data.frame(X[calib_idx, ]); names(X_calib) = names(X)
    Y_calib <- Y[calib_idx]
    
    all_X <- rbind(X, new_points)
    
    # ------------------ Formula ------------------
    create_formula <- function(X, smoothing = TRUE) {
      terms <- sapply(1:ncol(X), function(j) {
        if (length(unique(X[, j])) > 9 && smoothing) paste0("s(X", j, ")")
        else paste0("as.factor(X", j, ")")
      })
      reformulate(terms, response = "Y")
    }
    
    # ------------------ Fit model ------------------
    if (CQR_qr == 'qGAM') {
      formula <- create_formula(X_train)
      invisible(capture.output({
        fit_95 <- qgam(formula, data = data.frame(X_train, Y = Y_train), qu = (1 + desired_coverage) / 2)
        fit_50 <- qgam(formula, data = data.frame(X_train, Y = Y_train), qu = 0.5)
        fit_05 <- qgam(formula, data = data.frame(X_train, Y = Y_train), qu = (1 - desired_coverage) / 2)
      }))
      pred_95 <- predict(fit_95, newdata = all_X)
      pred_50 <- predict(fit_50, newdata = all_X)
      pred_05 <- predict(fit_05, newdata = all_X)
    } else {
      qrf <- quantregForest(x = X_train, y = Y_train, ntree = 1000, nodesize = 5)
      pred_95 <- predict(qrf, newdata = all_X, what = (1 + desired_coverage) / 2)
      pred_50 <- predict(qrf, newdata = all_X, what = 0.5)
      pred_05 <- predict(qrf, newdata = all_X, what = (1 - desired_coverage) / 2)
    }
    
    # ------------------ Enforce quantile order ------------------
    pred_95 <- pmax(pred_95, pred_50)
    pred_05 <- pmin(pred_05, pred_50)
    
    # ------------------ Calibration ------------------
    n_calib <- length(Y_calib)
    calib_idx_all <- calib_idx
    
    if (CQR_type == "multiplicative") {
      l <- pred_50 - pred_05
      u <- pred_95 - pred_50
      scores <- mapply(function(i) {
        idx <- calib_idx_all[i]
        max((pred_50[idx] - Y_calib[i]) / l[idx], (Y_calib[i] - pred_50[idx]) / u[idx])
      }, seq_len(n_calib))
      gamma <- quantile(scores, probs = desired_coverage)
      lower <- pred_50 - gamma * (pred_50 - pred_05)
      upper <- pred_50 + gamma * (pred_95 - pred_50)
      
    } else if (CQR_type == "additive") {
      scores <- mapply(function(i) {
        idx <- calib_idx_all[i]
        max(pred_05[idx] - Y_calib[i], Y_calib[i] - pred_95[idx])
      }, seq_len(n_calib))
      gamma <- quantile(scores, probs = desired_coverage)
      lower <- pred_05 - gamma
      upper <- pred_95 + gamma
      
    } else if (CQR_type == "naive") {
      scores <- abs(pred_50[calib_idx_all] - Y_calib)
      gamma <- quantile(scores, probs = desired_coverage)
      lower <- pred_50 - gamma
      upper <- pred_50 + gamma
      
    } else {
      stop("Unsupported CQR_type. Use 'additive', 'multiplicative', or 'naive'")
    }
    
    # ------------------ Output only new points ------------------
    m <- nrow(new_points)
    prediction_f <- tail(pred_50, m)
    lower <- tail(lower, m)
    upper <- tail(upper, m)
    
    list(
      hat_f = prediction_f,
      lower = lower,
      upper = upper,
      new_points = new_points
    )
  }
  
  
  
  PCS = function(X, Y, new_points, 
                 desired_coverage = 0.9, train_calib_split = 0.8, 
                 number_of_bootstraps = 30, symmetric_noise = TRUE) {
    
    # Data formatting and splitting
    data_original = data.frame(X = X, Y = Y)
    X = data.frame(X); Y = data.frame(Y)
    d = ncol(X); if (is.null(nrow(X))) d = 1
    names(data_original) = c(paste0("X", 1:d), "Y")
    new_points = data.frame(new_points); names(new_points) = c(paste0("X", 1:d))
    my_sequence = rbind(X, new_points)
    
    n = nrow(Y)
    train = 1:(n * train_calib_split); calib = (n * train_calib_split + 1):n
    X_train = X[train,]; Y_train = Y[train,]
    X_calib = X[calib,]; Y_calib = Y[calib,]
    n_train = length(Y_train); n_calib = length(Y_calib)
    index_sequences = (n_train + 1):(n_train + n_calib)
    
    create_formula = function(X, smoothing = TRUE) {
      form = "Y ~ "
      d = ncol(X); if (is.null(nrow(X))) d = 1
      for (i in 1:d) {
        if (smoothing && length(unique(X[, i])) > 9) {
          form = paste(form, paste0("s(X", i, ")+"))
        } else {
          form = paste(form, paste0("as.factor(X", i, ")+"))
        }
      }
      form = substr(form, 1, nchar(form) - 1)
      return(as.formula(form))
    }
    
    # Epistemic estimation
    X = X_train; Y = Y_train
    n = length(Y)
    data = data.frame(X = X_train, Y = Y_train)
    names(data) = c(paste0("X", 1:d), "Y")
    
    PCS_estimates_gam = list(); PCS_estimates_locfit = list()
    for (i in 1:number_of_bootstraps) {
      data_b = data[sample(1:n, n, replace = TRUE),]
      names(data_b) = c(paste0("X", 1:d), "Y")
      fit = gam(create_formula(X), sp = 0, data = data_b)
      preds = predict(fit, newdata = my_sequence)
      PCS_estimates_gam[[i]] = preds
      
      if (d < 3) {
        fit = locfit(create_formula(X, FALSE), data = data_b)
        preds = predict(fit, newdata = my_sequence)
        PCS_estimates_locfit[[i]] = preds
      }
    }
    
    PCS_mat = do.call(cbind, c(PCS_estimates_gam, PCS_estimates_locfit))
    prediction_f = apply(PCS_mat, 1, median)
    PCS_lower = apply(PCS_mat, 1, quantile, probs = 0.05)
    PCS_upper = apply(PCS_mat, 1, quantile, probs = 0.95)
    
    estimate_gamma = function(f = prediction_f) {
      score = sapply(1:n_calib, function(i) {
        index = index_sequences[i]
        max((f[index] - Y_calib[i]) / (f[index] - PCS_lower[index]),
            (Y_calib[i] - f[index]) / (PCS_upper[index] - f[index]))
      })
      return(quantile(score, desired_coverage))
    }
    
    gamma = estimate_gamma()
    lower = prediction_f - gamma * (prediction_f - PCS_lower)
    upper = prediction_f + gamma * (PCS_upper - prediction_f)
    
    # Restrict to prediction points
    m = nrow(new_points)
    prediction_f = prediction_f[(n+1):(n+m)]
    lower = lower[(n+1):(n+m)]
    upper = upper[(n+1):(n+m)]
    my_sequence = my_sequence[(n+1):(n+m),]
    
    return(list(hat_f = prediction_f,
                lower = lower, 
                upper = upper,
                new_points = new_points))
  }
  
  
  CLEAR = function(X, Y, new_points, 
                   minimum_lambda = 1, maximum_lambda = 10,  
                   number_of_bootstraps = 30, desired_coverage = 0.9, 
                   train_calib_split = 0.8, symmetric_noise = TRUE) {
    
    # Data formatting and splitting
    data_original = data.frame(X = X, Y = Y)
    X = data.frame(X); Y = data.frame(Y)
    d = ncol(X); if (is.null(nrow(X))) d = 1
    names(data_original) = c(paste0("X", 1:d), "Y")
    new_points = data.frame(new_points); names(new_points) = c(paste0("X", 1:d))
    my_sequence = rbind(X, new_points)
    
    n = nrow(Y)
    train = 1:(n * train_calib_split); calib = (n * train_calib_split + 1):n
    X_train = X[train,]; Y_train = Y[train,]
    X_calib = X[calib,]; Y_calib = Y[calib,]
    n_train = length(Y_train); n_calib = length(Y_calib)
    index_sequences = (n_train + 1):(n_train + n_calib)
    
    # Formula generator
    create_formula = function(X, smoothing = TRUE) {
      form = "Y ~ "
      d = ncol(X); if (is.null(nrow(X))) d = 1
      for (i in 1:d) {
        if (smoothing && length(unique(X[, i])) > 9) {
          form = paste(form, paste0("s(X", i, ")+"))
        } else {
          form = paste(form, paste0("as.factor(X", i, ")+"))
        }
      }
      form = substr(form, 1, nchar(form) - 1)
      return(as.formula(form))
    }
    
    # Step 1: Estimate aleatoric quantiles using bootstrapped qGAM
    X = X_train; Y = Y_train
    data = data.frame(X = X_train, Y = Y_train)
    names(data) = c(paste0("X", 1:d), "Y")
    
    q05_mat = q50_mat = q95_mat = list()
    invisible(capture.output({ #This should stop from printing to the console (verboise=FALSE doesnt work)
      for (i in 1:number_of_bootstraps) {
        data_b = data[sample(1:n, n, replace = TRUE),]
        newdata = if (d == 1) data.frame(X1 = c(my_sequence)) else my_sequence
        
        q95_mat[[i]] = predict(qgam(create_formula(X), data = data_b, qu = (1 + desired_coverage) / 2), newdata)
        q50_mat[[i]] = predict(qgam(create_formula(X), data = data_b, qu = 0.5), newdata)
        q05_mat[[i]] = predict(qgam(create_formula(X), data = data_b, qu = (1 - desired_coverage) / 2), newdata)
      }
    }))
    
    q95_mat = do.call(cbind, q95_mat)
    q50_mat = do.call(cbind, q50_mat)
    q05_mat = do.call(cbind, q05_mat)
    
    for (i in 1:nrow(q95_mat)) {
      for (j in 1:ncol(q95_mat)) {
        q95_mat[i, j] = max(q95_mat[i, j], q50_mat[i, j])
        q05_mat[i, j] = min(q50_mat[i, j], q05_mat[i, j])
      }
    }
    
    prediction_quant0.95 = rowMedians(q95_mat)
    prediction_quant0.5 = rowMedians(q50_mat)
    prediction_quant0.05 = rowMedians(q05_mat)
    
    # Step 2: Estimate epistemic uncertainty using gam and locfit
    PCS_mat = list()
    for (i in 1:number_of_bootstraps) {
      data_b = data[sample(1:n, n, replace = TRUE),]
      fit = gam(create_formula(X), sp = 0, data = data_b)
      PCS_mat[[length(PCS_mat) + 1]] = predict(fit, newdata = my_sequence)
      
      if (d < 3) {
        fit = locfit(create_formula(X, FALSE), data = data_b)
        PCS_mat[[length(PCS_mat) + 1]] = predict(fit, newdata = my_sequence)
      }
    }
    
    PCS_mat = do.call(cbind, PCS_mat)
    prediction_f_PCS = apply(PCS_mat, 1, median)
    PCS_lower = apply(PCS_mat, 1, quantile, 0.05)
    PCS_upper = apply(PCS_mat, 1, quantile, 0.95)
    
    # Step 3: Calibration for PCS_CQR
    estimate_gamma = function(lambda, f = prediction_f_PCS) {
      l = (prediction_quant0.5 - prediction_quant0.05) + lambda * (prediction_f_PCS - PCS_lower)
      u = (prediction_quant0.95 - prediction_quant0.5) + lambda * (PCS_upper - prediction_f_PCS)
      
      score = sapply(1:n_calib, function(i) {
        idx = index_sequences[i]
        max((f[idx] - Y_calib[i]) / l[idx], (Y_calib[i] - f[idx]) / u[idx])
      })
      quantile(score, desired_coverage)
    }
    
    estimate_optimal_lambda = function(min_lambda, max_lambda = 10) {
      losses = sapply(seq(min_lambda, max_lambda, by = 0.1), function(lambda) {
        gamma = estimate_gamma(lambda)
        l = prediction_f_PCS - gamma * ((prediction_quant0.5 - prediction_quant0.05) + lambda * (prediction_f_PCS - PCS_lower))
        u = prediction_f_PCS + gamma * ((prediction_quant0.95 - prediction_quant0.5) + lambda * (PCS_upper - prediction_f_PCS))
        
        mean(sapply(1:n_calib, function(i) {
          idx = index_sequences[i]
          0.05 * max(u[idx] - Y_calib[i], 0) +
            0.95 * max(Y_calib[i] - u[idx], 0) +
            0.95 * max(l[idx] - Y_calib[i], 0) +
            0.05 * max(Y_calib[i] - l[idx], 0)
        }))
      })
      seq(min_lambda, max_lambda, by = 0.1)[which.min(losses)]
    }
    
    lambda = estimate_optimal_lambda(minimum_lambda, maximum_lambda)
    gamma = estimate_gamma(lambda)
    
    PCS_CQR_lower = prediction_f_PCS - gamma * ((prediction_quant0.5 - prediction_quant0.05) + lambda * (prediction_f_PCS - PCS_lower))
    PCS_CQR_upper = prediction_f_PCS + gamma * ((prediction_quant0.95 - prediction_quant0.5) + lambda * (PCS_upper - prediction_f_PCS))
    
    # Keep only prediction points
    m = nrow(new_points)
    prediction_f_PCS = prediction_f_PCS[(n + 1):(n + m)]
    PCS_CQR_lower = PCS_CQR_lower[(n + 1):(n + m)]
    PCS_CQR_upper = PCS_CQR_upper[(n + 1):(n + m)]
    my_sequence = my_sequence[(n + 1):(n + m), ]
    
    return(list(hat_f = prediction_f_PCS,
                lower = PCS_CQR_lower,
                upper = PCS_CQR_upper,
                new_points = new_points))
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  ####################################################################################
  ############# Now we estimate prediction intervals for T=1 and T=0 and CATE ########
  ####################################################################################
  data0 = data[data$treatment==0,]
  X_untreated = data0[, 1:d]; Y_untreated = data0$Y
  data1 = data[data$treatment==1,]
  X_treated = data1[, 1:d]; Y_treated = data1$Y
  
  if(conformal=='CQR' & weighted_conformal==TRUE){
    estimate0 = CQR_weighted(X_untreated, Y_untreated, 
                             treatment = data$treatment,
                             X_all = X,  new_points,
                    desired_coverage=desired_coverage, 
                    our_weights = '1/e',
                    train_calib_split = train_calib_split, 
                    CQR_qr = CQR_qr, 
                    CQR_type = CQR_type)
    estimate1 = CQR_weighted(X_treated, Y_treated,
                             treatment = data$treatment,
                             X_all = X,    new_points,
                    desired_coverage=desired_coverage, 
                    our_weights = '1/(1-e)',
                    train_calib_split = train_calib_split, 
                    CQR_qr = CQR_qr, 
                    CQR_type = CQR_type)
  }
  if(conformal=='CQR' & weighted_conformal==FALSE){
    estimate0 = CQR(X_untreated, Y_untreated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    CQR_qr = CQR_qr, 
                    CQR_type = CQR_type)
    estimate1 = CQR(X_treated, Y_treated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    CQR_qr = CQR_qr, 
                    CQR_type = CQR_type)
  }
  if(conformal=='naive'){
    estimate0 = Naive_conformal(X_untreated, Y_untreated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    CQR_qr = CQR_qr)
    estimate1 = Naive_conformal(X_treated, Y_treated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    CQR_qr = CQR_qr)
  }
  if(conformal=='PCS'){
    estimate0 = PCS(X_untreated, Y_untreated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    number_of_bootstraps = number_of_bootstraps)
    estimate1 = PCS(X_treated, Y_treated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    number_of_bootstraps = number_of_bootstraps)
  }
  if(conformal=='CLEAR'){
    estimate0 = CLEAR(X_untreated, Y_untreated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    number_of_bootstraps = number_of_bootstraps)
    estimate1 = CLEAR(X_treated, Y_treated, new_points,
                    desired_coverage=desired_coverage, 
                    train_calib_split = train_calib_split, 
                    number_of_bootstraps = number_of_bootstraps)
  }
  

  
  hat_f0 = estimate0$hat_f
  lower0 = estimate0$lower
  upper0 = estimate0$upper
  lower_u0 = hat_f0 - estimate0$lower 
  upper_u0 = estimate0$upper - hat_f0
  my_sequence0 = estimate0$new_points
  
  # Extract second estimate details
  hat_f1 = estimate1$hat_f
  lower1 = estimate1$lower
  upper1 = estimate1$upper
  lower_u1 = hat_f1 - estimate1$lower 
  upper_u1 = estimate1$upper - hat_f1
  my_sequence1 = estimate1$new_points
  
  #combine sequences: if univariate then just sort and unique
  if(is.null(ncol(my_sequence0)))my_sequence = sort(unique(c(my_sequence0, my_sequence1)))
  if(!is.null(ncol(my_sequence0)))my_sequence = unique(rbind(my_sequence0, my_sequence1)) 
  
  
  
  if(CATE_estimate=='T-learner'){CATE = hat_f1 - hat_f0}
  
  if(CATE_estimate=='Causal_forest'){
    cf = causal_forest(X = as.matrix(X), 
                       Y = data_original$Y, 
                       W = as.numeric(treatment))
    
    CATE =  predict(cf,newdata = my_sequence)$predictions 
  }

  
  ########################################################################################
  ############## Finally, we compute the CATE confidence intervals #######################
  ########################################################################################
  if(add_confidence_intervals==TRUE){
    
    cate_hat_boot=list()
    for(k in 1:number_of_bootstraps){
      data_boot = data[sample(1:nrow(data_original), nrow(data_original), replace = TRUE),]
      
      if(CATE_estimate=='Causal_forest'){
      cate_hat_boot[[k]] = predict(causal_forest(X=as.matrix(data_boot[, 1:d]), 
                                                 Y = data_boot$Y, 
                                                 W = data_boot$treatment),
                                   newdata = my_sequence)$predictions
      }
      if(CATE_estimate=='T-learner'){
        data0 = data_boot[data_boot$treatment==0,]
        X_untreated = data0[, 1:d]; Y_untreated = data0$Y
        data1 = data_boot[data_boot$treatment==1,]
        X_treated = data1[, 1:d]; Y_treated = data1$Y
        
        estimate0 = Naive_conformal(X_untreated, Y_untreated, my_sequence,
                        desired_coverage=desired_coverage, 
                        train_calib_split = train_calib_split, 
                        CQR_qr = CQR_qr)
        estimate1 = Naive_conformal(X_treated, Y_treated, my_sequence,
                        desired_coverage=desired_coverage, 
                        train_calib_split = train_calib_split, 
                        CQR_qr = CQR_qr)
        cate_hat_boot[[k]] = estimate1$hat_f - estimate0$hat_f
      }
    } 
    #Produce 5% and 95% quantiles from all bootstraps
    cate_hat_boot = do.call(cbind, cate_hat_boot)
    cate_hat_boot = apply(cate_hat_boot, 1, quantile, probs = c(0.05, 0.5, 0.95))
    
    CATE=cate_hat_boot[2,]
    CATE_lower_CI = cate_hat_boot[2,] - cate_hat_boot[1,]
    CATE_upper_CI = cate_hat_boot[3,] - cate_hat_boot[2,]
    
    }


  
  
  if(add_confidence_intervals==FALSE){
    lower = CATE - rho_distance(lower_u0, upper_u1, rho)
    upper = CATE + rho_distance(upper_u0, lower_u1, rho)
  }else{
    lower = CATE - rho_distance(lower_u0, upper_u1, rho) - CI_addition(CATE_lower_CI, rho)
    upper = CATE + rho_distance(upper_u0, lower_u1, rho) + CI_addition(CATE_upper_CI, rho)
  }
  names(upper) = names(lower) = names(CATE) = paste0("new_point_", 1:nrow(data.frame(new_points)))
  return(list(new_points = new_points, CATE = CATE, lower = lower, upper = upper))
}







#Do you want to compute D_rho intervals for different Uncertainty estimation, or different CATE estimate? Its easy!
#input estimate0, estimate1 and CATE in the following form:
#   estimate0 is a data-frame with estimate0$f_hat;  estimate0$lower;  estimate0$upper; estimate0$new_points
#   estimate1 is a data-frame with estimate1$f_hat;  estimate1$lower;  estimate1$upper; estimate1$new_points
#   here, new_points are the points where we want to estimate ITEs (i.e. covariates in the test dataset)
#   CATE is a data-frame with CATE$CATE; (and CATE$lower; CATE$upper if we want to add confidence intervals)
#Output is a data-frame with lower and upper prediction intervals for ITE, together with CATE and new-points



D_rho_intervals_compute_yourself = function(estimate0, estimate1, CATE = NULL, 
                           rho = 0, CATE_CI = NULL) {
  rho_distance = function(x, y, rho)return(sqrt(x^2 + y^2 - 2*rho*x*y))
  CI_addition = function(x, rho) return( x*((1+rho)/2)^2 )       
  
  prediction_f0 = estimate0$f_hat
  lower0 = estimate0$lower
  upper0 = estimate0$upper
  lower_u0 = prediction_f0 - estimate0$lower 
  upper_u0 = estimate0$upper - prediction_f0
  my_sequence0 = estimate0$new_points
  
  # Extract second estimate details
  prediction_f1 = estimate1$f_hat
  lower1 = estimate1$lower
  upper1 = estimate1$upper
  lower_u1 = prediction_f1 - estimate1$lower 
  upper_u1 = estimate1$upper - prediction_f1
  my_sequence1 = estimate1$new_points
  
  #combine sequences
  if(is.null(ncol(my_sequence0)))new_points = sort(unique(c(my_sequence0, my_sequence1)))
  if(!is.null(ncol(my_sequence0)))new_points = unique(rbind(my_sequence0, my_sequence1)) 
  
  
  
  if(is.null(CATE)){CATE = prediction_f1 - prediction_f0}#We use simple t-estimator unless specified otherwise
  CATE = data.frame(CATE)
  
  
  if(is.null(CATE$lower)){
    lower = CATE - rho_distance(lower_u0, upper_u1, rho)
    upper = CATE + rho_distance(upper_u0, lower_u1, rho)
  }else{
    CI_lower = CATE - CATE$lower
    CI_upper = CATE$upper - CATE
    lower = CATE - rho_distance(lower_u0, upper_u1, rho) - CI_addition(CI_lower, rho)
    upper = CATE + rho_distance(upper_u0, lower_u1, rho) + CI_addition(CI_upper, rho)
  }
  
  return(list(new_points = new_points, CATE = CATE$CATE, lower = lower, upper = upper))
}



