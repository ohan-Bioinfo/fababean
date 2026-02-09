# First, create the CV summary file
cat > cv_summary.txt << EOF
K CV_Error
2 0.88910
3 1.07480
4 1.28890
5 1.30272
6 1.63901
7 1.70065
8 1.36837
EOF

# Now create the R visualization script
cat > visualize_cv.R << 'EOF'
#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)
library(dplyr)
library(gridExtra)

# Create the data
cv_data <- data.frame(
  K = 2:8,
  CV_Error = c(0.88910, 1.07480, 1.28890, 1.30272, 1.63901, 1.70065, 1.36837)
)

# Find best K
best_k <- cv_data$K[which.min(cv_data$CV_Error)]
best_cv_error <- min(cv_data$CV_Error)

cat("=== ADMIXTURE CROSS-VALIDATION RESULTS ===\n")
cat("Best K:", best_k, "(CV error =", best_cv_error, ")\n\n")
print(cv_data)

# Create the main CV error plot
p1 <- ggplot(cv_data, aes(x = K, y = CV_Error)) +
  geom_line(color = "navy", size = 1.2, alpha = 0.7) +
  geom_point(aes(color = K == best_k), size = 4) +
  geom_point(data = cv_data[cv_data$K == best_k, ], 
             aes(x = K, y = CV_Error), 
             color = "red", size = 6, shape = 1, stroke = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkblue"), guide = "none") +
  labs(
    title = paste("ADMIXTURE Cross-Validation Error\nBest K =", best_k),
    subtitle = paste("Lowest CV error =", best_cv_error),
    x = "Number of populations (K)",
    y = "Cross-validation error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95")
  ) +
  scale_x_continuous(breaks = 2:8) +
  annotate("text", x = best_k, y = best_cv_error, 
           label = paste("Best K =", best_k), 
           vjust = -2, color = "red", fontface = "bold", size = 5)

# Create a bar plot version
p2 <- ggplot(cv_data, aes(x = factor(K), y = CV_Error, fill = K == best_k)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = round(CV_Error, 4)), vjust = -0.5, fontface = "bold") +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"), guide = "none") +
  labs(
    title = "CV Error by K Value",
    x = "Number of populations (K)",
    y = "Cross-validation error"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "grey90")
  ) +
  ylim(0, max(cv_data$CV_Error) * 1.1)

# Create comparison table
cv_table <- cv_data %>%
  mutate(
    CV_Error = round(CV_Error, 5),
    Best = ifelse(K == best_k, "✓", ""),
    Rank = rank(CV_Error)
  ) %>%
  arrange(Rank)

table_plot <- tableGrob(cv_table[, c("K", "CV_Error", "Best")], 
                       rows = NULL,
                       theme = ttheme_minimal(
                         base_size = 12,
                         padding = unit(c(8, 4), "mm")
                       ))

# Combine plots
combined_plot <- grid.arrange(
  p1, p2, table_plot,
  ncol = 2,
  nrow = 2,
  layout_matrix = rbind(c(1, 1), c(2, 3)),
  top = textGrob("ADMIXTURE Cross-Validation Analysis", 
                 gp = gpar(fontsize = 18, fontface = "bold"))
)

# Save plots
ggsave("cv_analysis_comprehensive.png", combined_plot, width = 14, height = 10, dpi = 300)
ggsave("cv_error_plot.png", p1, width = 10, height = 8, dpi = 300)
ggsave("cv_bar_plot.png", p2, width = 8, height = 6, dpi = 300)

# Create individual admixture plots for key K values
plot_admixture <- function(K) {
  q_file <- paste0("Faba_chrOnly_pruned_numeric.", K, ".Q")
  if (!file.exists(q_file)) {
    cat("File not found:", q_file, "\n")
    return(NULL)
  }
  
  q_matrix <- read.table(q_file)
  
  # Read sample names
  fam_data <- read.table("Faba_chrOnly_pruned_numeric.fam", header = FALSE)
  sample_names <- fam_data$V1
  
  # Prepare data
  plot_data <- data.frame(
    Sample = rep(sample_names, each = K),
    Population = factor(rep(paste0("Pop", 1:K), times = length(sample_names))),
    Proportion = as.vector(t(as.matrix(q_matrix)))
  )
  
  p <- ggplot(plot_data, aes(x = factor(Sample, levels = unique(Sample)), 
                            y = Proportion, fill = Population)) +
    geom_bar(stat = "identity", width = 1) +
    labs(
      title = paste("Admixture Proportions - K =", K),
      x = "Samples",
      y = "Ancestry Proportion"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  if (K == best_k) {
    p <- p + labs(title = paste("Admixture Proportions - K =", K, "(BEST)")) +
      theme(plot.title = element_text(color = "red"))
  }
  
  return(p)
}

# Create directory for plots
dir.create("admixture_plots", showWarnings = FALSE)

# Plot admixture for best K and a few others
key_ks <- c(2, 3, 4, best_k)
key_ks <- unique(key_ks)

for (k in key_ks) {
  cat("Creating admixture plot for K =", k, "\n")
  p <- plot_admixture(k)
  if (!is.null(p)) {
    ggsave(paste0("admixture_plots/admixture_K", k, ".png"), p, 
           width = 14, height = 6, dpi = 300)
  }
}

# Create results summary
results_summary <- data.frame(
  Parameter = c("Best K", "CV Error", "Total Samples", "Total SNPs"),
  Value = c(best_k, 
            best_cv_error,
            nrow(read.table("Faba_chrOnly_pruned_numeric.fam", header = FALSE)),
            nrow(read.table("Faba_chrOnly_pruned_numeric.bim", header = FALSE)))
)

write.csv(cv_data, "cv_detailed_results.csv", row.names = FALSE)
write.csv(results_summary, "analysis_summary.csv", row.names = FALSE)

cat("\n" + rep("=", 60) + "\n")
cat("ANALYSIS COMPLETE\n")
cat(rep("=", 60) + "\n")
cat("Best K:", best_k, "with CV error =", best_cv_error, "\n")
cat("\nKey files generated:\n")
cat("- cv_analysis_comprehensive.png: Main analysis plot\n")
cat("- cv_error_plot.png: CV error line plot\n")
cat("- cv_bar_plot.png: CV error bar plot\n")
cat("- admixture_plots/: Directory with admixture bar plots\n")
cat("- cv_detailed_results.csv: Detailed CV results\n")
cat("- analysis_summary.csv: Analysis summary\n")
cat(rep("=", 60) + "\n")
EOF

# Run the visualization
chmod +x visualize_cv.R
/usr/bin/Rscript visualize_cv.R

# Create a quick text-based visualization
echo "=== ADMIXTURE CV ERROR SUMMARY ==="
echo "K    CV Error    Status"
echo "----------------------------"
for k in 2 3 4 5 6 7 8; do
  cv=$(grep "CV error (K=$k):" log_K${k}.out | awk '{print $4}')
  if [ "$k" -eq 2 ]; then
    echo "$k    $cv    ★ BEST"
  else
    echo "$k    $cv"
  fi
done

# Create a simple ASCII plot
echo ""
echo "=== CV ERROR TREND (ASCII) ==="
max_cv=1.70
for k in 2 3 4 5 6 7 8; do
  cv=$(grep "CV error (K=$k):" log_K${k}.out | awk '{print $4}')
  bars=$(echo "scale=0; ($cv / $max_cv) * 50" | bc)
  if [ "$k" -eq 2 ]; then
    marker="★"
  else
    marker=" "
  fi
  printf "K=%d %s |%s%*s (%.5f)\n" $k "$marker" "$(printf '%*s' $bars | tr ' ' '=')" $((50 - bars)) "" $cv
done