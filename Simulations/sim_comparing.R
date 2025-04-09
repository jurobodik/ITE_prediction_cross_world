#Simulations about the coverage of the D_rho intervals in comparison with other approaches such as Lei et al.

rho_vals <- c(1, 0.5, 0, -0.5, -1)
d = 1
n = 1000
n_iter = 50
data_generation = 'Synthetic Gaussian'



result_D_rho <- result_D_rho_with_CI <-result_D_rho_0.25 <- list()
result_lei_inexact <- result_lei_exact <- result_jonkers <- result_alaa <- list()
df_final = list()
for(rho_true in rho_vals){  
  
  
  for (k in 1:n_iter) {
    if(data_generation == 'Synthetic Gaussian'){
      n=n
      data_all <- data_synthetic(n = 10000, d = d, rho = rho_true, constant_propensity = FALSE)
      data <- data_all[1:n, ]
      data_test <- data_all[(n+1):nrow(data_all), ]
      new_points <- data.frame(data_test[, 1:d]); names(new_points) <- paste0('X', 1:d)
      if(d==1) X = data.frame(X1 = data[,1]) else X = data[, 1:d]
    }
    if(data_generation == 'IHDP'){
      d=25;n=600
      data_all <- IHDP_with_rho(rho_true)
      train_sample = sample(1:nrow(data_all), n)
      data <- data_all[train_sample, ]
      data_test <- data_all[-train_sample, ]
      new_points <- data_test[, 1:d]; names(new_points) <- paste0('X', 1:d)
    }
    if(data_generation!='Synthetic Gaussian' & data_generation!='IHDP') stop('data_generation must be either synthetic or IHDP')
    
    
    
    
    D_rho <- D_rho_intervals(X = data[, 1:d], Y = data$Y_obs, treatment = data$treatment, 
                                     new_points = new_points, rho = rho_true, 
                                     weighted_conformal = TRUE, add_confidence_intervals = FALSE)
    if(rho_true == -1) {D_rho_with_CI <-D_rho_0.25 <- D_rho}else{
      D_rho_with_CI <- D_rho_intervals(X = data[, 1:d], Y = data$Y_obs, treatment = data$treatment, 
                                               new_points = new_points, rho = rho_true, 
                                               weighted_conformal = TRUE, add_confidence_intervals = TRUE)
      
      D_rho_0.25 <- D_rho_intervals(X = data[, 1:d], Y = data$Y_obs, treatment = data$treatment, 
                                            new_points = new_points, rho = rho_true-0.25, 
                                            weighted_conformal = TRUE, add_confidence_intervals = FALSE)
    }
    Lei_inexact <- Lei_ITE(X = data[, 1:d], Y = data$Y_obs, 
                           T = data$treatment, new_points = new_points, 
                           exact = FALSE)
    Lei_exact <- Lei_ITE(X = data[, 1:d], Y = data$Y_obs, 
                         T = data$treatment, new_points = new_points, 
                         exact = TRUE)
    jonkers_and_alaa <- Jonkers_and_Alaa_wrappers_from_python(X_train = data[, 1:d], Y = data$Y_obs, 
                                                              T = data$treatment, new_points = new_points)
    jonkers <- jonkers_and_alaa$jonkers
    alaa <- jonkers_and_alaa$alaa
    
    
    if(d==1){X_new = data.frame(X1 = data_test[,1])}else{
      X_new <- data_test[, 1:d]}
    Y_new <- data_test$Y1 - data_test$Y0
    
    
    coverage_D_rho <- estimate_coverage(X_new, Y_new, 
                                        l = D_rho$lower, u = D_rho$upper)
    coverage_D_rho_with_CI <- estimate_coverage(X_new,  Y_new, 
                                                l = D_rho_with_CI$lower, 
                                                u = D_rho_with_CI$upper)
    coverage_D_rho_0.25 <- estimate_coverage(X_new,  Y_new, 
                                             l = D_rho_0.25$lower, 
                                             u = D_rho_0.25$upper)
    coverage_ITE_Lei_inexact <- estimate_coverage(X_new, Y_new, 
                                                  l= Lei_inexact$lower,
                                                  u= Lei_inexact$upper)
    coverage_ITE_Lei_exact <- estimate_coverage(X_new, Y_new, 
                                                l= Lei_exact$lower,
                                                u= Lei_exact$upper)
    coverage_ITE_jonkers <- estimate_coverage(X_new, Y_new, 
                                              l= jonkers$lower,
                                              u= jonkers$upper)
    coverage_ITE_alaa <- estimate_coverage(X_new, Y_new, 
                                           l= alaa$lower,
                                           u= alaa$upper)
    
    length_D_rho <- estimate_average_length(D_rho$lower, D_rho$upper)
    length_D_rho_with_CI <- estimate_average_length(D_rho_with_CI$lower, D_rho_with_CI$upper)
    length_D_rho_0.25 <- estimate_average_length(D_rho_0.25$lower, D_rho_0.25$upper)
    length_Lei_inexact <- estimate_average_length(Lei_inexact$lower, Lei_inexact$upper)
    length_Lei_exact <- estimate_average_length(Lei_exact$lower, Lei_exact$upper)
    length_jonkers <- estimate_average_length(jonkers$lower, jonkers$upper)
    length_alaa <- estimate_average_length(alaa$lower, alaa$upper)
    
    result_D_rho[[k]] <- c(coverage_D_rho, length_D_rho)
    result_D_rho_with_CI[[k]] <- c(coverage_D_rho_with_CI, length_D_rho_with_CI)
    result_D_rho_0.25[[k]] <- c(coverage_D_rho_0.25, length_D_rho_0.25)
    result_lei_inexact[[k]] <- c(coverage_ITE_Lei_inexact, length_Lei_inexact)
    result_lei_exact[[k]] <- c(coverage_ITE_Lei_exact, length_Lei_exact)
    result_jonkers[[k]] <- c(coverage_ITE_jonkers, length_jonkers)
    result_alaa[[k]] <- c(coverage_ITE_alaa, length_alaa)
    
    cat('Iteration:', k, ' and rho = ', rho_true ,  '\n', 'D_rho:', result_D_rho[[k]], '\n', 'D_rho_with_CI:', result_D_rho_with_CI[[k]], '\n', 'D_rho_0.25:', result_D_rho_0.25[[k]], '\n',
        'Lei_inexact:', result_lei_inexact[[k]], '\n', 'Lei_exact:', result_lei_exact[[k]], '\n', 'Jonkers:', result_jonkers[[k]], '\n', 'Alaa:', result_alaa[[k]], '\n')
  }
  
  # Combine results
  methods <- list(
    D_rho = result_D_rho,
    D_rho_with_CI = result_D_rho_with_CI,
    D_rho_0.25 = result_D_rho_0.25,
    Lei_inexact = result_lei_inexact,
    Lei_exact = result_lei_exact,
    CMC = result_jonkers,
    DR = result_alaa
  )
  
  df <- map2_dfr(methods, names(methods), ~ {
    mat <- do.call(rbind, .x)
    data.frame(Coverage = mat[,1], Width = mat[,2], Method = .y)
  })
  
  df$Rho <- paste0("rho = ", rho_true)
  
  #add df to df_final
  df_final = rbind(df_final, df)
}

all_results <- df_final










################################## PLOTTING ##################################
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)

rho_vals <- rho_vals[5:1]#from back to front
# Set factor levels for methods
all_results$Method <- factor(all_results$Method,
                             levels = c('D_rho_with_CI','D_rho_0.25', "D_rho", "CMC", "DR", "Lei_inexact","Lei_exact")
)
method_labels <- c(
  "D_rho"        = expression(D[rho]),
  "D_rho_with_CI" = expression(D[rho]*" + CI"),
  "D_rho_0.25"   = expression(D[rho-0.25]),
  "CMC"          = "CMC",
  "DR"           = "DR",
  "Lei_inexact"  = "Lei (inexact)",
  "Lei_exact"    = "Lei (naive)"
)

# Define rho values and expression labels
rho_labels <- list(
  `-1` = expression(rho == -1),
  `-0.5` = expression(rho == -0.5),
  `0` = expression(rho == 0),
  `0.5` = expression(rho == 0.5),
  `1` = expression(rho == 1)
)

highlight_methods <- c("D_rho", "D_rho_with_CI", "D_rho_0.25")
custom_colors <- c(
  "D_rho" = "#1f77b4",         # blue
  "D_rho_with_CI" = "#d62728", # red
  "D_rho_0.25" = "#2ca02c",    # green
  "CMC" = "gray40",
  "DR" = "gray60",
  "Lei_inexact" = "gray80",
  "Lei_exact" = "gray99"
)

generate_plots <- function(rho_val, hide_x_axis = TRUE) {
  df_rho <- filter(all_results, Rho == paste0("rho = ", rho_val))
  
  # Compute y-positions for the highlight
  method_levels <- levels(df_rho$Method)
  highlight_indices <- which(method_levels %in% highlight_methods)
  ymin_highlight <- min(highlight_indices) - 0.5
  ymax_highlight <- max(highlight_indices) + 0.5
  
  highlight_rect <- annotate("rect",
                             xmin = -Inf, xmax = Inf,
                             ymin = ymin_highlight, ymax = ymax_highlight,
                             alpha = 0.1, fill = "lightblue")
  
  coverage_plot <- ggplot(df_rho, aes(x = Coverage, y = Method, fill = Method)) +
    highlight_rect +
    geom_boxplot(outlier.size = 1.2) +
    geom_vline(xintercept = 0.9, color = "black", size = 1, linetype = 'solid') +
    xlim(0.67, 1) +
    scale_y_discrete(labels = method_labels)+
    ggtitle(rho_labels[[as.character(rho_val)]]) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  width_plot <- ggplot(df_rho, aes(x = Width, y = Method, fill = Method)) +
    highlight_rect +
    geom_boxplot(outlier.size = 1.2) +
    xlim(0, 15) +
    scale_y_discrete(labels = method_labels)+
    ggtitle(rho_labels[[as.character(rho_val)]]) +
    scale_fill_manual(values = custom_colors) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  # Remove x-axis text/ticks unless it's the last row
  if (hide_x_axis) {
    coverage_plot <- coverage_plot +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    width_plot <- width_plot +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  
  list(coverage = coverage_plot, width = width_plot)
}






# Create plots, only show x-axis on the bottom row
plot_pairs <- lapply(seq_along(rho_vals), function(i) {
  generate_plots(rho_vals[i], hide_x_axis = i != length(rho_vals))
})

# Extract coverage and width plots
coverage_plots <- lapply(plot_pairs, `[[`, "coverage")
width_plots <- lapply(plot_pairs, `[[`, "width")

# Arrange with shared guides and synced layout
left_col <- wrap_plots(coverage_plots, ncol = 1, guides = "collect")
right_col <- wrap_plots(width_plots, ncol = 1, guides = "collect")

# Column titles
title_cov <- wrap_elements(full = textGrob("Coverage", gp = gpar(fontsize = 14, fontface = "bold")))
title_wid <- wrap_elements(full = textGrob("Width", gp = gpar(fontsize = 14, fontface = "bold")))

# Final combined layout
final_plot <- (title_cov | title_wid) / (left_col | right_col) +
  plot_layout(heights = c(0.06, 1))

# Display
print(final_plot)
