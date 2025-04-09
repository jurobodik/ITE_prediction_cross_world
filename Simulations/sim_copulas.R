#Simulations about the coverage of the D_rho intervals for different copula and marginal distributions of epsilon

copula_marginal_combos <- list(
  list(copula = "gaussian", marginal = "gaussian"),
  list(copula = "t",        marginal = "t"),
  list(copula = "clayton",  marginal = "chisq"),
  list(copula = "gumbel",   marginal = "laplace")
)

rhos <- c(-0.75, -0.5, 0.25, 0, 0.25, 0.5, 0.75)
dimension_of_X <- c(10)
reps <- 50





results <- list()
results_with_CI <- list()

counter <- 1
for (combo in copula_marginal_combos) {
  cop <- combo$copula
  marg <- combo$marginal
  
  for (r in rhos) {
    for (d in dimension_of_X) {
      for (rep in 1:reps) {
        cat(sprintf("Running: %d / %d\r", counter, length(copula_marginal_combos) * length(rhos) * length(dimension_of_X) * reps)); flush.console()
        counter <- counter + 1
        
        df_all <- data_synthetic(n = 10000, d = d, rho = r,
                                 copula_type = cop,
                                 marginal = marg,
                                 constant_propensity = FALSE)
        df <- df_all[1:1000, ]
        df_test <- df_all[1001:10000, ]
        df$ITE <- df$Y1 - df$Y0
        df_test$ITE <- df_test$Y1 - df_test$Y0
        
        for (add_CI in c(FALSE, TRUE)) {
          res <- tryCatch({
            D_rho_intervals(
              X = data.frame(df[,1:d]),
              Y = df$Y_obs,
              treatment = df$treatment,
              new_points = if (d == 1) data.frame(X = df_test$X) else df_test[,1:d],
              rho = r,
              add_confidence_intervals = add_CI,
              conformal = "CQR",
              weighted_conformal = TRUE,
              CQR_qr = "auto",
              CATE_estimate = "T-learner"
            )
          }, error = function(e) return(NULL))
          
          if (!is.null(res)) {
            df_test$lower <- res$lower
            df_test$upper <- res$upper
            df_test$covered <- (df_test$ITE >= df_test$lower) & (df_test$ITE <= df_test$upper)
            coverage <- mean(df_test$covered)
            
            out_row <- data.frame(
              copula = cop,
              marginal = marg,
              rho = r,
              d = d,
              rep = rep,
              coverage = coverage,
              ci = add_CI
            )
            
            if (add_CI) {
              results_with_CI[[length(results_with_CI) + 1]] <- out_row
            } else {
              results[[length(results) + 1]] <- out_row
            }
          }
        }
      }
    }
  }
}

results_all <- bind_rows(results, results_with_CI)





