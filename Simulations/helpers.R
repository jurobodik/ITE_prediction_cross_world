#Some helpful functions needed for simulations, such as:
#generating benchmark datasets
#uploading IHDP_data and generating new counterfactual with different rho
#wrapper for Lei et al method
#wrapper for Jonkers et al method and Alaa et al method from python using Reticulate library

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

#You need to close 'https://github.com/predict-idlab/cct-cmc.git' and use virtual environment to use python in R
Jonkers_and_Alaa_wrappers_from_python <- function(X_train, Y, T, new_points, confidence = 0.9) {
  library(reticulate)
  #system("git clone https://github.com/predict-idlab/cct-cmc.git")
  #virtualenv_create("cctcmc_env")
  #use_python("C:/Users/jurob/OneDrive/Dokumenty/.virtualenvs/cctcmc_env/Scripts/python.exe", required = TRUE)
  use_virtualenv("cctcmc_env", required = TRUE)
  #virtualenv_install("cctcmc_env", packages = c("numpy", 'scipy', "pandas", "scikit-learn", "crepes", "crepes-weighted", "seaborn", "ipykernel", "jupyterlab", "tqdm"))
  
  # Import Python modules
  np <- import("numpy", delay_load = TRUE)
  sklearn <- import("sklearn.ensemble", delay_load = TRUE)
  rf <- sklearn$RandomForestRegressor
  
  
  
  #############################################################
  ######################CMC Jonkers############################
  #############################################################
  # Import the CMC_T_Learner class
  repo_path <- "cct-cmc"
  cmc_path <- file.path(repo_path, "src/cmc_metalearners")
  CMC <- import_from_path("cmc_metalearners", path = cmc_path)
  CMC_T_Learner <- CMC$CMC_T_Learner
  
  # Convert data to Python types
  X_py <- r_to_py(as.matrix(X_train))
  Y_py <- r_to_py(np$array(as.numeric(as.vector(Y))), convert = FALSE)
  T_py <- r_to_py(np$array(as.integer(as.vector(T))), convert = FALSE)
  
  # Fit the learner
  learner <- CMC_T_Learner(
    rf(n_estimators = as.integer(100)),
    rf(n_estimators = as.integer(100)),
    pseudo_MC = TRUE,
    MC_samples = as.integer(100)
  )
  
  py_run_string("import warnings; warnings.filterwarnings('ignore')")
  learner$fit(X_py, Y_py, T_py)
  # Predict intervals
  intervals <- learner$predict_int( r_to_py(as.matrix(new_points)), confidence = confidence)
  intervals_CMC = py_to_r(intervals)
  
  
  #############################################################
  ######################Alaa DR################################
  #############################################################
  
  py_run_string(paste0("import sys; sys.path.append('", file.path("cct-cmc", 'src'), "')"))
  DR <- import("conformal_metalearners.drlearner")
  
  
  X_py <- r_to_py(as.matrix(X_train))
  Y_py <- r_to_py(np$array(as.numeric(Y)), convert = FALSE)
  T_py <- r_to_py(np$array(as.integer(T)), convert = FALSE)
  ps_py <- r_to_py(np$array(as.numeric(rep(mean(T), length(T)))), convert = FALSE)
  
  # Instantiate
  dr <- DR$conformalMetalearner(
    n_folds = as.integer(5),
    alpha = 1- confidence,
    base_learner = "RF",
    quantile_regression = FALSE,
    metalearner = "DR"
  )
  
  # Fit
  dr$fit(X_py, T_py, Y_py, ps_py)
  dr$conformalize(
    alpha = 1- confidence,
    X_calib = X_py,
    W_calib = T_py,
    Y_calib = Y_py,
    pscores_calib = ps_py
  )
  
  
  pred_r <- py_to_r(dr$predict(r_to_py(as.matrix(new_points))))  # returns list of [point_est, lower, upper]
  
  alaa = list(upper =  pred_r[[3]], lower = pred_r[[2]], new_points = as.numeric(new_points[[1]]))
  jonkers = list(upper = intervals_CMC[,2], lower = intervals_CMC[,1], new_points = as.numeric(new_points[[1]]))
  
  return(list(jonkers = jonkers, alaa = alaa))
}



estimate_coverage = function(X_new, Y_new, l, u){ # l = lower interval, u=upper interval
  coverage = 0
   for(i in 1:nrow(X_new)){
     if(l[i] < Y_new[i] & Y_new[i] < u[i]){
       coverage = coverage + 1
     }}

  return(coverage/nrow(X_new))
}


estimate_average_length = function(l, u){ return(median(u-l ))}



Lei_ITE = function(X, Y, T, new_points, exact = TRUE, outfun = "quantRF"){
  
  X = matrix(unlist(X), nrow = nrow(data.frame(X)))
  
  new_points <- matrix(unlist(new_points), ncol=ncol(data.frame(new_points)))
  
  CIfun <- conformalIte(X, Y, T, alpha = 0.1, 
                        algo = "nest", exact = exact, type = "CQR",
                        quantiles = c(0.05, 0.95), 
                        outfun = "quantRF", useCV = FALSE)
  result = CIfun(new_points)
  return(list(upper = result$upper, lower = result$lower, new_points = new_points))
}



#################Generating benchmark datasets##############

data_synthetic <- function(n = 1000, 
                           d = 2, 
                           rho, 
                           sigma_1 = 1, 
                           sigma_2 = 4, 
                           constant_propensity = FALSE, 
                           copula_type = "gaussian", #choices: "gaussian", "t", "frank"
                           marginal = "gaussian"){ #choices: "gaussian", "laplace", "t", "laplace", "chisq"
  # generate random 1D and 2D functions
  random_function_1d <- function(freq = 0.15) {
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
  
  
  generate_errors <- function(n, rho, sigma_1, sigma_2, copula_type = "gaussian", marginal = "gaussian") {
    library(copula)
    
    qlaplace <- function(p) {
      ifelse(p < 0.5, log(2 * p), -log(2 * (1 - p)))
    }
    
    transform_marginal <- function(u_vec, type) {
      switch(type,
             "gaussian" = qnorm(u_vec),
             "t" = qt(u_vec, df = 3),
             "laplace" = qlaplace(u_vec),
             "chisq" = qchisq(u_vec, df = 3),
             stop("Unsupported marginal type.")
      )
    }
    
    find_tau_for_rho <- function(target_rho, copula_type = "frank", marginal = "gaussian",
                                 sigma_1 = 1, sigma_2 = 1, n = 1e4) {
      
      obj_fn <- function(tau) {
        # Cap tau to avoid extreme values
        tau <- max(min(tau, 0.99), -0.99)
        
        theta <- switch(copula_type,
                        "frank" = iTau(frankCopula(), tau),
                        stop("Unsupported copula type.")
        )
        
        cop <- switch(copula_type,
                      "frank" = frankCopula(param = theta, dim = 2)
        )
        
        u <- rCopula(n, cop)
        x <- transform_marginal(u[, 1], marginal) * sqrt(sigma_1)
        y <- transform_marginal(u[, 2], marginal) * sqrt(sigma_2)
        
        cor(x, y) - target_rho
      }
      
      # Check feasibility
      tryCatch({
        uniroot(obj_fn, lower = -0.95, upper = 0.95)$root
      }, error = function(e) {
        stop("Could not match the desired Pearson correlation with the specified copula and marginal.")
      })
    }
    
    # Determine copula object based on copula_type
    cop <- switch(copula_type,
                  "gaussian" = normalCopula(param = rho, dim = 2),
                  "t"        = tCopula(param = rho, dim = 2, df = 3),
                  "frank"    = {
                    tau <- find_tau_for_rho(rho, copula_type = "frank", marginal = marginal,
                                            sigma_1 = sigma_1, sigma_2 = sigma_2)
                    frankCopula(param = iTau(frankCopula(), tau), dim = 2)
                  },
                  stop("Unsupported copula type.")
    )
    
    u <- rCopula(n, cop)
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




IHDP_with_rho <- function(rho, setup = 'B', load_csv_file=FALSE) {
  if(load_csv_file){ihdp  =  read.csv("ihdp_data.csv")}
  data = ihdp
  data$treatment = as.integer(as.logical(data$treatment) )
  X = as.matrix(data[, -c(1,2,3,4,5)])
  
  p <- ncol(X)
  n <- nrow(X)
  
  if (setup == 'A') {
    # ---------------------------
    # Response Surface A (Linear)
    # ---------------------------
    
    beta_A <- sample(0:4, p, replace = TRUE, prob = c(0.5, 0.2, 0.15, 0.1, 0.05))
    
    mu0 <- X %*% beta_A
    mu1 <- mu0 + 4  # constant treatment effect
    
  } else if (setup == 'B') {
    # ------------------------------
    # Response Surface B (Nonlinear)
    # ------------------------------
    
    beta_B <- sample(c(0, 0.1, 0.2, 0.3, 0.4), p, replace = TRUE, prob = c(0.6, 0.1, 0.1, 0.1, 0.1))
    
    W <- matrix(0.5, nrow = n, ncol = p)
    
    mu0_raw <- exp((X + W) %*% beta_B)
    mu1_raw <- X %*% beta_B
    
    # Center treatment effect to have ATE = 4
    omega_B <- mean(mu1_raw - mu0_raw) - 4
    mu0 <- mu0_raw
    mu1 <- mu1_raw - omega_B
    
  } else {
    stop("setup must be 'A' or 'B'")
  }
  
  # Correlated noise
  epsilons <- MASS::mvrnorm(n, mu = c(0, 0), 
                            Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  
  Y0 <- as.numeric(mu0 + epsilons[,1])
  Y1 <- as.numeric(mu1 + epsilons[,2])
  Y_obs = ifelse(data$treatment == 1, Y1, Y0)
  X1 = data$x1; X2 = data$x2; X3 = data$x3; X4 = data$x4; X5 = data$x5; X6 = data$x6; X7 = data$x7; X8 = data$x8; X9 = data$x9; X10 = data$x10
  X11= data$x11; X12 = data$x12; X13 = data$x13; X14 = data$x14; X15 = data$x15; X16 = data$x16; X17 = data$x17; X18 = data$x18; X19 = data$x19; X20 = data$x20
  X21= data$x21; X22 = data$x22; X23 = data$x23; X24 = data$x24; X25 = data$x25
  return(data.frame(X1=X1,X2=X2,X3=X3,X4=X4,X5=X5,X6=X6,X7=X7,X8=X8,X9=X9,X10=X10,
                    X11=X11,X12=X12,X13=X13,X14=X14,X15=X15,X16=X16,X17=X17,X18=X18,X19=X19,X20=X20,
                    X21=X21,X22=X22,X23=X23,X24=X24,X25=X25, treatment = data$treatment, Y_obs = Y_obs, Y0 = Y0, Y1 = Y1))
}


upload_twins_dataset = function(){
  
  # Read the datasets
  x <- read_csv("https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/master/datasets/TWINS/twin_pairs_X_3years_samesex.csv")
  y <- read_csv("https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/master/datasets/TWINS/twin_pairs_Y_3years_samesex.csv")
  t <- read_csv("https://raw.githubusercontent.com/AMLab-Amsterdam/CEVAE/master/datasets/TWINS/twin_pairs_T_3years_samesex.csv")
  
  # Define column names
  lighter_columns <- c('pldel', 'birattnd', 'brstate', 'stoccfipb', 'mager8',
                       'ormoth', 'mrace', 'meduc6', 'dmar', 'mplbir', 'mpre5', 'adequacy',
                       'orfath', 'frace', 'birmon', 'gestat10', 'csex', 'anemia', 'cardiac',
                       'lung', 'diabetes', 'herpes', 'hydra', 'hemo', 'chyper', 'phyper',
                       'eclamp', 'incervix', 'pre4000', 'preterm', 'renal', 'rh', 'uterine',
                       'othermr', 'tobacco', 'alcohol', 'cigar6', 'drink5', 'crace',
                       'data_year', 'nprevistq', 'dfageq', 'feduc6', 'infant_id_0',
                       'dlivord_min', 'dtotord_min', 'bord_0',
                       'brstate_reg', 'stoccfipb_reg', 'mplbir_reg')
  
  heavier_columns <- c('pldel', 'birattnd', 'brstate', 'stoccfipb', 'mager8',
                       'ormoth', 'mrace', 'meduc6', 'dmar', 'mplbir', 'mpre5', 'adequacy',
                       'orfath', 'frace', 'birmon', 'gestat10', 'csex', 'anemia', 'cardiac',
                       'lung', 'diabetes', 'herpes', 'hydra', 'hemo', 'chyper', 'phyper',
                       'eclamp', 'incervix', 'pre4000', 'preterm', 'renal', 'rh', 'uterine',
                       'othermr', 'tobacco', 'alcohol', 'cigar6', 'drink5', 'crace',
                       'data_year', 'nprevistq', 'dfageq', 'feduc6', 'infant_id_1',
                       'dlivord_min', 'dtotord_min', 'bord_1',
                       'brstate_reg', 'stoccfipb_reg', 'mplbir_reg')
  
  # Process data
  data <- list()
  for (i in seq_len(nrow(t))) {
    weights <- unlist(t[i, 2:3])
    if (any(weights >= 2000)) next
    
    lighter <- as.numeric(x[i, lighter_columns])
    heavier <- as.numeric(x[i, heavier_columns])
    
    lighter_row <- c(lighter, weights[1], 0, y[i, 2])
    heavier_row <- c(heavier, weights[2], 1, y[i, 3])
    
    data[[length(data) + 1]] <- lighter_row
    data[[length(data) + 1]] <- heavier_row
  }
  
  # Final column names
  cols <- c('pldel', 'birattnd', 'brstate', 'stoccfipb', 'mager8',
            'ormoth', 'mrace', 'meduc6', 'dmar', 'mplbir', 'mpre5', 'adequacy',
            'orfath', 'frace', 'birmon', 'gestat10', 'csex', 'anemia', 'cardiac',
            'lung', 'diabetes', 'herpes', 'hydra', 'hemo', 'chyper', 'phyper',
            'eclamp', 'incervix', 'pre4000', 'preterm', 'renal', 'rh', 'uterine',
            'othermr', 'tobacco', 'alcohol', 'cigar6', 'drink5', 'crace',
            'data_year', 'nprevistq', 'dfageq', 'feduc6',
            'infant_id', 'dlivord_min', 'dtotord_min', 'bord',
            'brstate_reg', 'stoccfipb_reg', 'mplbir_reg', 'wt', 'treatment', 'outcome')
  
  df <- as.data.frame(do.call(rbind, data))
  colnames(df) <- cols
  
  # Convert types
  df$treatment <- as.numeric(df$treatment)
  df$outcome <- as.numeric(df$outcome)
  
  # Fill missing values
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      col[is.na(col)] <- mean(col, na.rm = TRUE)
    } else {
      col[is.na(col)] <- as.character(stats::na.omit(col)[1])
    }
    return(col)
  })
  
  
  d = ncol(df) - 2
  df[, 1:d] <- lapply(df[, 1:d], function(x) as.numeric(as.character(x)))
  
  Y_obs <- df$outcome
  T <- df$treatment
  
  Y_1 <- df$outcome[df$treatment == 1]
  Y_0 <- df$outcome[df$treatment == 0]
  X = df[ seq(1, nrow(df) , by = 2), 1:d]; names(X) = c(paste0("X", 1:d))
  return(list(X = X, Y_1 = Y_1, Y_0 = Y_0))
}


