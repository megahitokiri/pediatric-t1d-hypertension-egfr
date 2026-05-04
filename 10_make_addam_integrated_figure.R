#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(cowplot)
  library(scales)
})

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && nzchar(x[1])) x else y
}

script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) %||% getwd()
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
if (basename(script_dir) != "reproducible_pipeline") {
  script_dir <- "/Users/jose/Documents/PHD experimental Medicine/HTA Paper/reproducible_pipeline"
}

tables_dir <- file.path(script_dir, "results", "addam_genetics", "tables")
fig_dir <- file.path(script_dir, "results", "paper_figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.08), hjust = 0),
      plot.subtitle = element_text(color = "#4d5a5f", size = rel(0.9), hjust = 0),
      axis.title = element_text(face = "bold", color = "#1a242d"),
      axis.text = element_text(color = "#3d454a"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(7, 12, 7, 12)
    )
}

domain_colors <- c(
  Clinical = "#2B83BA",
  Immune = "#D64745",
  Genetic = "#8E6BBE",
  Anthropometry = "#4DAF4A"
)

feature_labels <- c(
  X96GAD = "GAD65 autoantibody",
  age_BP = "Age at BP assessment",
  age_atb = "Age at autoantibody assessment",
  X6.29840255_A1 = "Top SNP dosage",
  cor_insulin = "Insulin dose",
  Cluster_diabete_history = "Diabetes duration/history",
  bmiz = "BMI z-score",
  GRS = "T1D GRS/PRS",
  IA2 = "IA-2 autoantibody",
  Cluster_sex_bin = "Sex",
  Cluster_A1c = "HbA1c cluster"
)

feature_domain <- function(x) {
  case_when(
    x %in% c("X96GAD", "IA2", "ZnT8") ~ "Immune",
    x %in% c("GRS", "X6.29840255_A1") ~ "Genetic",
    x %in% c("bmiz", "bmi", "bmip") ~ "Anthropometry",
    TRUE ~ "Clinical"
  )
}

perf <- read_csv(file.path(tables_dir, "Table_ADDAM_ML_LOOCV_Performance.csv"), show_col_types = FALSE) %>%
  mutate(
    auc = as.numeric(auc),
    model = factor(model, levels = rev(model[order(auc)]))
  )

p_auc <- ggplot(perf, aes(x = auc, y = model)) +
  geom_segment(aes(x = 0.50, xend = auc, yend = model), color = "#d7e0e6", linewidth = 1.1) +
  geom_vline(xintercept = 0.50, linetype = "dashed", color = "#7a8a8f") +
  geom_point(color = "#2B83BA", size = 3.2) +
  geom_text(aes(label = sprintf("%.3f", auc)), hjust = -0.25, size = 3.4) +
  scale_x_continuous(limits = c(0.50, 0.72), breaks = c(0.50, 0.60, 0.70)) +
  labs(
    title = "A. Model discrimination",
    subtitle = "Leave-one-out cross-validation in ADDAM",
    x = "AUC",
    y = NULL
  ) +
  theme_paper()

incremental <- read_csv(file.path(tables_dir, "Table_ADDAM_Incremental_LOOCV.csv"), show_col_types = FALSE) %>%
  mutate(
    delta_auc_vs_clinical = as.numeric(delta_auc_vs_clinical),
    model_set = factor(model_set, levels = rev(model_set[order(delta_auc_vs_clinical)]))
  )

p_incremental <- ggplot(incremental, aes(x = delta_auc_vs_clinical, y = model_set)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#7a8a8f") +
  geom_segment(aes(x = 0, xend = delta_auc_vs_clinical, yend = model_set), color = "#e5d4ad", linewidth = 1.1) +
  geom_point(color = "#D9822B", size = 3.2) +
  geom_text(aes(label = sprintf("%+.3f", delta_auc_vs_clinical)), hjust = ifelse(incremental$delta_auc_vs_clinical >= 0, -0.25, 1.15), size = 3.2) +
  scale_x_continuous(limits = c(-0.03, 0.012), breaks = c(-0.03, -0.02, -0.01, 0, 0.01)) +
  labs(
    title = "B. Added value beyond clinical core",
    subtitle = "PRS and ancestry did not materially improve AUC",
    x = "Delta AUC vs clinical core",
    y = NULL
  ) +
  theme_paper()

importance <- read_csv(file.path(tables_dir, "Table_ADDAM_ML_FeatureImportance.csv"), show_col_types = FALSE) %>%
  slice_head(n = 10) %>%
  mutate(
    label = recode(feature, !!!feature_labels),
    domain = factor(feature_domain(feature), levels = names(domain_colors)),
    label = factor(label, levels = rev(label[order(median_importance)]))
  )

p_importance <- ggplot(importance, aes(x = median_importance, y = label, color = domain)) +
  geom_errorbar(aes(xmin = ci_2_5, xmax = ci_97_5), width = 0, linewidth = 1.1, alpha = 0.85, orientation = "y") +
  geom_point(size = 3.3) +
  scale_color_manual(values = domain_colors, drop = FALSE) +
  scale_x_continuous(limits = c(0, 105), breaks = c(0, 25, 50, 75, 100)) +
  labs(
    title = "C. Biologic feature-importance profile",
    subtitle = "Immune and age-related variables carry the strongest signal",
    x = "Median bootstrap importance (95% interval)",
    y = NULL,
    color = "Domain"
  ) +
  theme_paper()

prs <- read_csv(file.path(tables_dir, "Table_ADDAM_PRS_Univariate.csv"), show_col_types = FALSE) %>%
  filter(feature %in% c("IA2", "ZnT8", "bmiz", "duration_diabetes_per_year", "GRS", "age_atb_per_year", "age_BP_per_year", "X96GAD")) %>%
  mutate(
    label = recode(
      feature,
      IA2 = "IA-2 autoantibody",
      ZnT8 = "ZnT8 autoantibody",
      bmiz = "BMI z-score",
      duration_diabetes_per_year = "Diabetes duration (per year)",
      GRS = "T1D GRS/PRS",
      age_atb_per_year = "Age at autoantibody assessment",
      age_BP_per_year = "Age at BP assessment",
      X96GAD = "GAD65 autoantibody"
    ),
    domain = factor(feature_domain(gsub("_per_year", "", feature)), levels = names(domain_colors)),
    label = factor(label, levels = rev(label[order(odds_ratio)]))
  )

p_or <- ggplot(prs, aes(x = odds_ratio, y = label, color = domain)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#7a8a8f") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 1.1, alpha = 0.85, orientation = "y") +
  geom_point(size = 3.3) +
  scale_color_manual(values = domain_colors, drop = FALSE) +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30), labels = label_number()) +
  labs(
    title = "D. Selected univariate associations",
    subtitle = "PRS was contextual, not a stand-alone BP discriminator",
    x = "Odds ratio (log scale)",
    y = NULL,
    color = "Domain"
  ) +
  theme_paper()

top <- plot_grid(p_auc, p_incremental, labels = NULL, ncol = 2, rel_widths = c(1, 1.08))
bottom <- plot_grid(p_importance, p_or, labels = NULL, ncol = 2, rel_widths = c(1.05, 1))
body <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1.08))

title <- ggdraw() +
  draw_label(
    "Figure 4. ADDAM deep phenotyping links elevated BP to immune and age-related risk signals",
    x = 0.01, y = 0.72, hjust = 0, fontface = "bold", size = 18, color = "#101820"
  ) +
  draw_label(
    "Clinical, autoantibody, anthropometric, ancestry, SNP dosage, and T1D genetic risk score features were analyzed as a mechanistic extension of the UDDB population signal.",
    x = 0.01, y = 0.28, hjust = 0, size = 11, color = "#4d5a5f"
  )

footer <- ggdraw() +
  draw_label(
    "ADDAM n=288; elevated BP n=48. Model performance reflects leave-one-out cross-validation. Colors: clinical blue; immune red; genetic purple; anthropometry green.",
    x = 0.01, y = 0.5, hjust = 0, size = 9.5, color = "#5a656b"
  )

fig <- plot_grid(title, body, footer, ncol = 1, rel_heights = c(0.12, 1, 0.06))

ggsave(
  filename = file.path(fig_dir, "Figure4_ADDAM_Integrated_Model_Biology.png"),
  plot = fig,
  width = 13,
  height = 8.7,
  dpi = 300,
  bg = "white"
)
