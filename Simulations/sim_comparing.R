#Simulations about the coverage of the CW_rho intervals in comparison with other approaches such as Lei et al.

set.seed(123)
rho_vals <- c(1, 0.5, 0, -0.5, -1)
d = 1
n = 2000
n_iter = 50
data_generation = 'Synthetic Gaussian'



result_CW_rho <- result_CW_rho_with_CI <-result_CW_rho_0.25 <- list()
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
    
    
    
    
    CW_rho <- CW_rho_intervals(X = data[, 1:d], Y = data$Y_obs, treatment = data$treatment, 
                                     new_points = new_points, rho = rho_true, 
                                     weighted_conformal = TRUE, add_confidence_intervals = FALSE)
    if(rho_true == -1) {CW_rho_with_CI <-CW_rho_0.25 <- CW_rho}else{
      CW_rho_with_CI <- CW_rho_intervals(X = data[, 1:d], Y = data$Y_obs, treatment = data$treatment, 
                                               new_points = new_points, rho = rho_true, 
                                               weighted_conformal = TRUE, add_confidence_intervals = TRUE)
      
      CW_rho_0.25 <- CW_rho_intervals(X = data[, 1:d], Y = data$Y_obs, treatment = data$treatment, 
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
    
    
    coverage_CW_rho <- estimate_coverage(X_new, Y_new, 
                                        l = CW_rho$lower, u = CW_rho$upper)
    coverage_CW_rho_with_CI <- estimate_coverage(X_new,  Y_new, 
                                                l = CW_rho_with_CI$lower, 
                                                u = CW_rho_with_CI$upper)
    coverage_CW_rho_0.25 <- estimate_coverage(X_new,  Y_new, 
                                             l = CW_rho_0.25$lower, 
                                             u = CW_rho_0.25$upper)
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
    
    length_CW_rho <- estimate_average_length(CW_rho$lower, CW_rho$upper)
    length_CW_rho_with_CI <- estimate_average_length(CW_rho_with_CI$lower, CW_rho_with_CI$upper)
    length_CW_rho_0.25 <- estimate_average_length(CW_rho_0.25$lower, CW_rho_0.25$upper)
    length_Lei_inexact <- estimate_average_length(Lei_inexact$lower, Lei_inexact$upper)
    length_Lei_exact <- estimate_average_length(Lei_exact$lower, Lei_exact$upper)
    length_jonkers <- estimate_average_length(jonkers$lower, jonkers$upper)
    length_alaa <- estimate_average_length(alaa$lower, alaa$upper)
    
    result_CW_rho[[k]] <- c(coverage_CW_rho, length_CW_rho)
    result_CW_rho_with_CI[[k]] <- c(coverage_CW_rho_with_CI, length_CW_rho_with_CI)
    result_CW_rho_0.25[[k]] <- c(coverage_CW_rho_0.25, length_CW_rho_0.25)
    result_lei_inexact[[k]] <- c(coverage_ITE_Lei_inexact, length_Lei_inexact)
    result_lei_exact[[k]] <- c(coverage_ITE_Lei_exact, length_Lei_exact)
    result_jonkers[[k]] <- c(coverage_ITE_jonkers, length_jonkers)
    result_alaa[[k]] <- c(coverage_ITE_alaa, length_alaa)
    
    cat('Iteration:', k, ' and rho = ', rho_true ,  '\n', 'CW_rho:', result_CW_rho[[k]], '\n', 'CW_rho_with_CI:', result_CW_rho_with_CI[[k]], '\n', 'CW_rho_0.25:', result_CW_rho_0.25[[k]], '\n',
        'Lei_inexact:', result_lei_inexact[[k]], '\n', 'Lei_exact:', result_lei_exact[[k]], '\n', 'Jonkers:', result_jonkers[[k]], '\n', 'Alaa:', result_alaa[[k]], '\n')
  }
  
  # Combine results
  methods <- list(
    CW_rho = result_CW_rho,
    CW_rho_with_CI = result_CW_rho_with_CI,
    CW_rho_0.25 = result_CW_rho_0.25,
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

################now run it for d=1 and d=15, save the files and then combine them
#write.csv(all_results, file = "sim1_results_d1.csv", row.names = TRUE)
#write.csv(all_results, file = "sim1_results_d15.csv", row.names = TRUE)
#write.csv(all_results, file = "sim1_results_IHDP.csv", row.names = TRUE)

#r1 = read.csv("sim1_results_d1.csv", row.names = 1);r1$d <- 1;
#r2 = read.csv("sim1_results_d15.csv", row.names = 1);r2$d <- 15
#all_results <- rbind(r1, r2)
#write.csv(all_results, file = "sim1_results_final_combined.csv", row.names = TRUE)



#all_results <- read.csv("sim1_results_final_combined.csv", row.names = 1)
#all_results_IHDP <- read.csv("sim1_results_IHDP.csv", row.names = 1)








################################## PLOTTING ##################################
library(ggplot2)
library(dplyr)
library(patchwork)
library(grid)
library(ggpattern)



generate_plots_all_tohether_final = function(all_results){
  
  
  rho_vals <- c(-1, -0.5, 0, 0.5, 1)
  
  all_results$Method <- factor(all_results$Method,
                               levels = c("Lei_exact", "Lei_inexact", "DR", "CMC", "D_rho", "D_rho_0.25", "D_rho_with_CI"))
  
  method_labels <- c(
    "D_rho"        = expression(CW(rho)),
    "D_rho_with_CI" = expression(CW(rho)*"+CI"),
    "D_rho_0.25"   = expression(CW("misspec. "*rho)),
    "CMC"          = "CMC",
    "DR"           = "DR",
    "Lei_inexact"  = "Lei (inexact)",
    "Lei_exact"    = "Lei (naive)"
  )
  
  
  rho_labels <- list(
    `1` = expression(rho == 1),
    `0.5` = expression(rho == 0.5),
    `0` = expression(rho == 0),
    `-0.5` = expression(rho == -0.5),
    `-1` = expression(rho == -1)
  )
  
  custom_colors <- c(
    "D_rho" = "#e41a1c",          # strong red
    "D_rho_with_CI" = "#d62728",  # medium red
    "D_rho_0.25" = "#f6615f",     # light red
    "CMC" = "gray40",
    "DR" = "gray60",
    "Lei_inexact" = "gray80",
    "Lei_exact" = "gray99"
  )
  
  
  generate_plots_overlay <- function(rho_val, hide_x_axis = TRUE) {
    df_rho <- filter(all_results, Rho == paste0("rho = ", rho_val))
    show_legend <- rho_val == 0
    method_levels <- levels(df_rho$Method)
    y_indices <- seq_along(method_levels)
    
    # Grey background for first three rows (our proposed methods)
    # Light blue background for every second row (striping)
    stripes <- lapply(y_indices[y_indices %% 2 == 0], function(i) {
      annotate("rect",
               xmin = -Inf, xmax = Inf,
               ymin = i - 0.5, ymax = i + 0.5,
               fill = "lightblue", alpha = 0.3)
    })
    
    
    
    
    df_rho$pattern <- factor(df_rho$d, levels = c(15, 1))
    pattern_vals <- c(`1` = "circle", `15` = "none")
    
    coverage_plot <- ggplot(df_rho, aes(
      x = Coverage,
      y = Method,
      fill = Method,
      pattern = pattern
    )) +
      stripes +
      ggpattern::geom_boxplot_pattern(
        position = position_dodge(width = 0.6),
        width = 0.5, outlier.size = 0.8,
        color = "black", linewidth = 0.4,
        pattern_fill = "gray20",
        pattern_density = 0.7,
        pattern_spacing = 0.025,
        pattern_key_scale_factor = 0.5
      ) +
      geom_vline(xintercept = 0.9, color = "black", size = 1) +
      xlim(0.67, 1) +
      scale_y_discrete(labels = method_labels) +
      scale_fill_manual(values = custom_colors, guide = "none") +
      scale_pattern_manual(
        values = pattern_vals,
        breaks = c("1", "15"),
        labels = c("1" = "d = 1", "15" = "d = 15"),
        name = if (show_legend) "Dimension" else NULL
      ) +
      theme_minimal(base_size = 11) +
      guides() +
      theme(
        legend.position = if (show_legend) "right" else "none",
        legend.box = "vertical",
        legend.box.just = "center",
        legend.margin = ggplot2::margin(t = 0, r = 0, b = 500, l = 0, unit = "pt"),
        axis.title = element_blank(),
        plot.title = element_blank()
      )
    
    
    
    
    width_plot <- ggplot(df_rho, aes(
      x = Width,
      y = Method,
      fill = Method,
      pattern = pattern
    )) +
      stripes +
      ggpattern::geom_boxplot_pattern(
        position = position_dodge(width = 0.6),
        width = 0.5, outlier.size = 0.8,
        color = "black", linewidth = 0.4,
        pattern_fill = "gray20",
        pattern_density = 0.7,
        pattern_spacing = 0.025,
        pattern_key_scale_factor = 0.5
      ) +
      xlim(0, 15) +
      scale_y_discrete(labels = method_labels) +
      scale_fill_manual(values = custom_colors, guide = "none") +
      scale_pattern_manual(
        values = pattern_vals,
        breaks = c("1", "15"),
        labels = c("1" = "d = 1", "15" = "d = 15"),
        name = NULL
      ) +
      theme_minimal(base_size = 11) +
      guides() +
      theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank()
      )
    
    if (hide_x_axis) {
      coverage_plot <- coverage_plot +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      width_plot <- width_plot +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    
    list(coverage = coverage_plot, width = width_plot)
  }
  
  # Add x-axis to the top row of plots
  # Re-generate the first row of plots without hiding x-axis
  plot_pairs <- lapply(seq_along(rho_vals), function(i) {
    generate_plots_overlay(rho_vals[i], hide_x_axis = !(i %in% c(1, length(rho_vals))))
  })
  
  
  coverage_plots <- lapply(plot_pairs, `[[`, "coverage")
  width_plots <- lapply(plot_pairs, `[[`, "width")
  
  rho_labels_grobs <- lapply(rho_vals, function(rho_val) {
    label_expr <- rho_labels[[as.character(rho_val)]]
    ggplot() + 
      annotate("text", x = 1, y = 1, label = label_expr, parse = TRUE, size = 5) +
      theme_void()
  })
  
  rows_with_labels <- Map(function(lbl, cov, wid) {
    lbl | cov | wid
  }, lbl = rho_labels_grobs, cov = coverage_plots, wid = width_plots)
  
  title_rho <- ggplot() + theme_void()
  title_cov <- ggplot() +
    annotate("text", x = 1, y = 1, label = "Coverage", size = 5, fontface = "bold") +
    theme_void()
  title_wid <- ggplot() +
    annotate("text", x = 1, y = 1, label = "Width", size = 5, fontface = "bold") +
    theme_void()
  
  header <- title_rho | title_cov | title_wid
  
  row_height <- 1
  row_heights <- rep(row_height, length(rows_with_labels))
  header_height <- 0.06
  
  full_plot <- wrap_plots(
    list(header, rows_with_labels[[1]], rows_with_labels[[2]],
         rows_with_labels[[3]], rows_with_labels[[4]], rows_with_labels[[5]]),
    ncol = 1
  )
  
  layout_matrix <- patchwork::area(
    t = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    l = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    b = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    r = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)
  )
  
  rows_split <- lapply(seq_along(rows_with_labels), function(i) {
    list(rho_labels_grobs[[i]], coverage_plots[[i]], width_plots[[i]])
  })
  
  plots_flat <- c(
    list(title_rho, title_cov, title_wid),
    do.call(c, rows_split)
  )
  
  final_plot <- wrap_plots(plots_flat, design = layout_matrix) +
    plot_layout(widths = c(1.5, 5, 5), heights = c(0.3, rep(1, 5)))


# Display
print(final_plot)
}



generate_plots_all_tohether_final(all_results)
ggsave("sim1_results_SIMULATED.pdf", plot = final_plot, width = 8, height = 9, dpi = 500)





generate_plots_IHDP = function(all_results_IHDP){
  
  
  
  rho_vals <- c(-1, -0.5, 0, 0.5, 1)
  
  all_results$Method <- factor(all_results$Method,
                               levels = c("Lei_exact", "Lei_inexact", "DR", "CMC", "CW_rho", "CW_rho_0.25", "CW_rho_with_CI"))
  
  method_labels <- c(
    "CW_rho"        = expression(CW(rho)),
    "CW_rho_with_CI" = expression(CW(rho)*"+CI"),
    "CW_rho_0.25"   = expression(CW("misspec. "*rho)),
    "CMC"          = "CMC",
    "DR"           = "DR",
    "Lei_inexact"  = "Lei (inexact)",
    "Lei_exact"    = "Lei (naive)"
  )
  
  rho_labels <- list(
    `1` = expression(rho == 1),
    `0.5` = expression(rho == 0.5),
    `0` = expression(rho == 0),
    `-0.5` = expression(rho == -0.5),
    `-1` = expression(rho == -1)
  )
  
  custom_colors <- c(
    "CW_rho" = "#e41a1c",          # strong red
    "CW_rho_with_CI" = "#d62728",  # medium red
    "CW_rho_0.25" = "#f6615f",     # light red
    "CMC" = "gray40",
    "DR" = "gray60",
    "Lei_inexact" = "gray80",
    "Lei_exact" = "gray99"
  )
  
  
  generate_plots <- function(rho_val, hide_x_axis = TRUE) {
    df_rho <- filter(all_results, Rho == paste0("rho = ", rho_val))
    method_levels <- levels(df_rho$Method)
    y_indices <- seq_along(method_levels)
    
    blue_stripes <- lapply(y_indices[y_indices %% 2 == 0], function(i) {
      annotate("rect", xmin = -Inf, xmax = Inf,
               ymin = i - 0.5, ymax = i + 0.5,
               fill = "lightblue", alpha = 0.3)
    })
    
    coverage_plot <- ggplot(df_rho, aes(x = Coverage, y = Method, fill = Method)) +
      blue_stripes +
      geom_boxplot(outlier.size = 1.2) +
      geom_vline(xintercept = 0.9, color = "black", size = 1) +
      xlim(0.67, 1) +
      scale_y_discrete(labels = method_labels) +
      scale_fill_manual(values = custom_colors) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        axis.title = element_blank(),
        plot.title = element_blank()
      )
    
    width_plot <- ggplot(df_rho, aes(x = Width, y = Method, fill = Method)) +
      blue_stripes +
      geom_boxplot(outlier.size = 1.2) +
      xlim(0, 15) +
      scale_y_discrete(labels = method_labels) +
      scale_fill_manual(values = custom_colors) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_blank()
      )
    
    if (hide_x_axis) {
      coverage_plot <- coverage_plot +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      width_plot <- width_plot +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    
    list(coverage = coverage_plot, width = width_plot)
  }
  
  plot_pairs <- lapply(seq_along(rho_vals), function(i) {
    generate_plots(rho_vals[i], hide_x_axis = !(i %in% c(1, length(rho_vals))))
  })
  
  coverage_plots <- lapply(plot_pairs, `[[`, "coverage")
  width_plots <- lapply(plot_pairs, `[[`, "width")
  
  rho_labels_grobs <- lapply(rho_vals, function(rho_val) {
    label_expr <- rho_labels[[as.character(rho_val)]]
    ggplot() +
      annotate("text", x = 1, y = 1, label = label_expr, parse = TRUE, size = 5) +
      theme_void()
  })
  
  rows_with_labels <- Map(function(lbl, cov, wid) {
    lbl | cov | wid
  }, lbl = rho_labels_grobs, cov = coverage_plots, wid = width_plots)
  
  title_rho <- ggplot() + theme_void()
  title_cov <- ggplot() +
    annotate("text", x = 1, y = 1, label = "Coverage", size = 5, fontface = "bold") +
    theme_void()
  title_wid <- ggplot() +
    annotate("text", x = 1, y = 1, label = "Width", size = 5, fontface = "bold") +
    theme_void()
  
  header <- title_rho | title_cov | title_wid
  
  full_plot <- wrap_plots(
    list(header, rows_with_labels[[1]], rows_with_labels[[2]],
         rows_with_labels[[3]], rows_with_labels[[4]], rows_with_labels[[5]]),
    ncol = 1
  )
  
  layout_matrix <- patchwork::area(
    t = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    l = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3),
    b = c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6),
    r = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)
  )
  
  rows_split <- lapply(seq_along(rows_with_labels), function(i) {
    list(rho_labels_grobs[[i]], coverage_plots[[i]], width_plots[[i]])
  })
  
  plots_flat <- c(
    list(title_rho, title_cov, title_wid),
    do.call(c, rows_split)
  )
  
  final_plot <- wrap_plots(plots_flat, design = layout_matrix) +
    plot_layout(widths = c(1.5, 5, 5), heights = c(0.3, rep(1, 5)))
  
  # Display
  print(final_plot)
}


#generate_plots_IHDP(all_results)
#ggsave("sim1_results_IHDP.pdf", plot = final_plot, width = 8, height = 9, dpi = 500)






