#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(gridExtra)
  library(grid)
})

load_cv_data <- function() {
  if (file.exists("cv_summary.txt")) {
    cv <- tryCatch(
      read.table("cv_summary.txt", header = TRUE),
      error = function(e) NULL
    )
    if (!is.null(cv) && nrow(cv) > 0 && all(c("K", "CV_Error") %in% colnames(cv))) {
      return(cv)
    }
  }

  logs <- list.files(pattern = "^log_K[0-9]+\\.out$")
  rows <- data.frame(K = integer(), CV_Error = numeric())
  if (length(logs) == 0) {
    stop("No cv_summary.txt or log_K*.out files found.")
  }

  for (log_file in logs) {
    lines <- readLines(log_file, warn = FALSE)
    cv_line <- grep("CV error \\(K=", lines, value = TRUE)
    if (length(cv_line) == 0) next
    k_value <- as.integer(sub(".*CV error \\(K=([0-9]+)\\):.*", "\\1", cv_line[1]))
    cv_value <- as.numeric(sub(".*CV error \\(K=[0-9]+\\):\\s*([0-9.]+).*", "\\1", cv_line[1]))
    rows <- rbind(rows, data.frame(K = k_value, CV_Error = cv_value))
  }

  rows %>% arrange(K)
}

plot_admixture <- function(k_value, sample_names, best_k) {
  q_file <- paste0("Faba_chrOnly_pruned_numeric.", k_value, ".Q")
  if (!file.exists(q_file)) return(NULL)

  q_matrix <- read.table(q_file)
  if (nrow(q_matrix) != length(sample_names)) return(NULL)

  plot_data <- data.frame(
    Sample = rep(sample_names, each = k_value),
    Population = factor(rep(paste0("Pop", seq_len(k_value)), times = length(sample_names))),
    Proportion = as.vector(t(as.matrix(q_matrix)))
  )

  title_text <- paste("Admixture Proportions - K =", k_value)
  if (k_value == best_k) title_text <- paste(title_text, "(BEST)")

  ggplot(plot_data, aes(x = factor(Sample, levels = unique(Sample)), y = Proportion, fill = Population)) +
    geom_bar(stat = "identity", width = 1) +
    labs(title = title_text, x = "Samples", y = "Ancestry Proportion") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

cv_data <- load_cv_data() %>% arrange(K)
if (nrow(cv_data) == 0) stop("No valid CV data found.")

best_k <- cv_data$K[which.min(cv_data$CV_Error)]
best_cv <- min(cv_data$CV_Error)

if (!dir.exists("plots")) dir.create("plots", recursive = TRUE)
if (!dir.exists("plots/Figure")) dir.create("plots/Figure", recursive = TRUE)

cat("=== ADMIXTURE CROSS-VALIDATION RESULTS ===\n")
cat(sprintf("Best K: %d (CV error = %.5f)\n\n", best_k, best_cv))
print(cv_data)

p_line <- ggplot(cv_data, aes(x = K, y = CV_Error)) +
  geom_line(color = "navy", linewidth = 1.1) +
  geom_point(aes(color = K == best_k), size = 3.5) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "darkblue"), guide = "none") +
  labs(
    title = paste("ADMIXTURE CV Error (Best K =", best_k, ")"),
    x = "Number of populations (K)",
    y = "Cross-validation error"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  scale_x_continuous(breaks = sort(unique(cv_data$K)))

p_bar <- ggplot(cv_data, aes(x = factor(K), y = CV_Error, fill = K == best_k)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.4f", CV_Error)), vjust = -0.4, size = 3.3) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "steelblue"), guide = "none") +
  labs(title = "CV Error by K", x = "K", y = "CV Error") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

cv_table <- cv_data %>%
  mutate(CV_Error = round(CV_Error, 5), Best = ifelse(K == best_k, "yes", ""))
table_plot <- tableGrob(cv_table, rows = NULL, theme = ttheme_minimal(base_size = 10))

combined <- grid.arrange(
  p_line, p_bar, table_plot,
  ncol = 2,
  nrow = 2,
  layout_matrix = rbind(c(1, 1), c(2, 3)),
  top = textGrob("ADMIXTURE Cross-Validation Summary", gp = gpar(fontsize = 16, fontface = "bold"))
)

ggsave("plots/cv_error_highlighted.png", combined, width = 14, height = 9, dpi = 300)
ggsave("plots/cv_error_highlighted.pdf", combined, width = 14, height = 9)

if (file.exists("Faba_chrOnly_pruned_numeric.fam")) {
  fam <- read.table("Faba_chrOnly_pruned_numeric.fam", header = FALSE)
  sample_names <- fam$V1
  selected_k <- unique(c(best_k, 2, 3, 4))
  for (k_value in selected_k) {
    p <- plot_admixture(k_value, sample_names, best_k)
    if (is.null(p)) next
    ggsave(
      paste0("plots/Figure/Admixture_K", k_value, ".png"),
      p,
      width = 14,
      height = 6,
      dpi = 300
    )
    ggsave(
      paste0("plots/Figure/Admixture_K", k_value, ".pdf"),
      p,
      width = 14,
      height = 6
    )
  }
}

write.csv(cv_data, "plots/admixture_detailed_summary.csv", row.names = FALSE)

cat("\n============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat(sprintf("Best K: %d with CV error %.5f\n", best_k, best_cv))
cat("============================================================\n")
