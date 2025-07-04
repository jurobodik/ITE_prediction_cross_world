#Simulations about the coverage of the CW_rho intervals for different copula and marginal distributions of epsilon

set.seed(123)
copula_marginal_combos <- list(
  list(copula = "gaussian", marginal = "gaussian"),
  list(copula = "t",        marginal = "t"),
  list(copula = "t",        marginal = "gaussian"),
  list(copula = "frank",  marginal = "laplace")
)

rhos <- c(-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75)
dimension_of_X <- c(1)
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
                                 constant_propensity = TRUE)
        df <- df_all[1:1000, ]
        df_test <- df_all[1001:10000, ]
        df$ITE <- df$Y1 - df$Y0
        df_test$ITE <- df_test$Y1 - df_test$Y0
        
        for (add_CI in c(FALSE, TRUE)) {
          res <- tryCatch({
            CW_rho_intervals(
              X = data.frame(df[,1:d]),
              Y = df$Y_obs,
              treatment = df$treatment,
              new_points = if (d == 1) data.frame(X = df_test$X) else df_test[,1:d],
              rho = r,
              add_confidence_intervals = add_CI,
              conformal = "CQR",
              weighted_conformal = FALSE,
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



###################### Plotting ####################
generate_final_plot <- function(results_all) {
  library(ggplot2)
  library(dplyr)
  
  # Set Method
  results_all <- results_all %>%
    mutate(
      Method = ifelse(ci, "CW(rho)+CI", "CW(rho)"),
      Method = factor(Method, levels = c("CW(rho)", "CW(rho)+CI"))
    )
  
  # Full copula and marginal names
  copula_names <- c(gaussian = "Gaussian", t = "Student-t", clayton = "Student-t", gumbel = "Frank", frank = "Frank")
  marginal_names <- c(gaussian = "Gaussian", t = "Student-t", chisq = "Gaussian", laplace = "Laplace")
  
  results_all <- results_all %>%
    mutate(
      copula_full = copula_names[copula],
      marginal_full = marginal_names[marginal],
      col_label = paste0("Copula: ", copula_full, "\nMarginal: ", marginal_full)
    )
  
  # Automatically set desired order from copula_marginal_combos
  desired_order <- sapply(copula_marginal_combos, function(combo) {
    paste0(
      "Copula: ", copula_names[[combo$copula]],
      "\nMarginal: ", marginal_names[[combo$marginal]]
    )
  })
  results_all$col_label <- factor(results_all$col_label, levels = desired_order)
  
  # rho facet (parsed math labels)
  results_all <- results_all %>%
    mutate(rho_lab = factor(paste0("rho==", rho), levels = paste0("rho==", sort(unique(rho)))))
  
  # Labels and colors
  method_labels <- c(
    "CW(rho)" = expression(CW(rho)),
    "CW(rho)+CI" = expression(CW(rho)*"+CI")
  )
  custom_colors <- c(
    "CW(rho)" = "#1f77b4",
    "CW(rho)+CI" = "#d62728"
  )
  
  # Plot
  p <- ggplot(results_all, aes(x = coverage, y = Method, fill = Method)) +
    geom_boxplot(
      position = position_dodge(width = 0.6),
      width = 0.5,
      outlier.size = 0.8,
      color = "black",
      linewidth = 0.4
    ) +
    facet_grid(rows = vars(rho_lab), cols = vars(col_label),
               labeller = labeller(rho_lab = label_parsed, col_label = label_value)) +
    geom_vline(xintercept = 0.9, color = "black", size = 1) +
    scale_x_continuous(limits = c(0.8, 1.0), breaks = c(0.8, 0.9)) +
    scale_fill_manual(values = custom_colors) +
    scale_y_discrete(labels = method_labels) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      strip.text.x = element_text(size = 10, face = "bold"),
      strip.text.y = element_text(size = 10),
      legend.position = "none",
      axis.title = element_blank()
    ) +
    labs(title = expression(''))
  
  print(p)
  return(p)
}





final_plot <- generate_final_plot(results_all)
#write.csv(results_all, "results_copula.csv", row.names = FALSE)
ggsave("sim_copula_figure.pdf", plot = final_plot, width = 7.5, height = 9.5, dpi = 500)











