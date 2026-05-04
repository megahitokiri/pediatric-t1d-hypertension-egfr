script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("dplyr", "readr", "ggplot2", "scales", "patchwork", "tidyr", "forcats"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(patchwork))

dirs <- ensure_pipeline_dirs(script_dir)
paper_figures_dir <- file.path(dirs$results_dir, "paper_figures")
supplemental_figures_dir <- file.path(dirs$results_dir, "supplemental_figures")
dir.create(paper_figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supplemental_figures_dir, recursive = TRUE, showWarnings = FALSE)

save_paper_plot <- function(plot, filename, width, height, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(paper_figures_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

copy_if_exists <- function(from, to_dir, to_name = basename(from)) {
  if (file.exists(from)) {
    file.copy(from, file.path(to_dir, to_name), overwrite = TRUE)
    return(invisible(TRUE))
  }
  invisible(FALSE)
}

fmt_n_pct <- function(n, pct) {
  paste0(n, " (", sprintf("%.1f", pct), "%)")
}

table1 <- readr::read_csv(file.path(dirs$tables_dir, "Table1_T1D_T1DImputed_U26.csv"), show_col_types = FALSE)
table2_class <- readr::read_csv(file.path(dirs$tables_dir, "Table2_HeightInformed_Pediatric_BP_Classification.csv"), show_col_types = FALSE)
table7_summary <- readr::read_csv(file.path(dirs$tables_dir, "Table7_PediatricHTN_summary.csv"), show_col_types = FALSE)
table8_qc <- readr::read_csv(file.path(dirs$tables_dir, "Table8_eGFR_HeightMatching_QC.csv"), show_col_types = FALSE)
model_df <- load_derived_rds("height_informed_pediatric_bp_model_df", dirs)
creat_labs_u25 <- load_derived_rds("creat_labs_u25", dirs)

addam_qc_path <- file.path(dirs$results_dir, "addam_genetics", "tables", "Table_ADDAM_Cohort_QC.csv")
addam_ml_path <- file.path(dirs$results_dir, "addam_genetics", "tables", "Table_ADDAM_ML_LOOCV_Performance.csv")
addam_imp_path <- file.path(dirs$results_dir, "addam_genetics", "tables", "Table_ADDAM_ML_FeatureImportance.csv")
addam_inc_path <- file.path(dirs$results_dir, "addam_genetics", "tables", "Table_ADDAM_Incremental_LOOCV.csv")
addam_univ_path <- file.path(dirs$results_dir, "addam_genetics", "tables", "Table_ADDAM_PRS_Univariate.csv")
addam_qc <- if (file.exists(addam_qc_path)) readr::read_csv(addam_qc_path, show_col_types = FALSE) else NULL
addam_ml <- if (file.exists(addam_ml_path)) readr::read_csv(addam_ml_path, show_col_types = FALSE) else NULL
addam_imp <- if (file.exists(addam_imp_path)) readr::read_csv(addam_imp_path, show_col_types = FALSE) else NULL
addam_inc <- if (file.exists(addam_inc_path)) readr::read_csv(addam_inc_path, show_col_types = FALSE) else NULL
addam_univ <- if (file.exists(addam_univ_path)) readr::read_csv(addam_univ_path, show_col_types = FALSE) else NULL

pretty_addam_feature <- function(x) {
  x <- gsub("^X96GAD$", "GAD65 autoantibody", x)
  x <- gsub("^IA2$", "IA-2 autoantibody", x)
  x <- gsub("^ZnT8$", "ZnT8 autoantibody", x)
  x <- gsub("^GRS$", "T1D GRS/PRS", x)
  x <- gsub("^age_BP$", "Age at BP assessment", x)
  x <- gsub("^age_BP_per_year$", "Age at BP assessment (per year)", x)
  x <- gsub("^age_atb$", "Age at autoantibody assessment", x)
  x <- gsub("^age_atb_per_year$", "Age at autoantibody assessment (per year)", x)
  x <- gsub("^duration_diabetes_per_year$", "Diabetes duration (per year)", x)
  x <- gsub("^Cluster_A1c$", "HbA1c cluster", x)
  x <- gsub("^Cluster_sex_bin$", "Sex", x)
  x <- gsub("^cluster_shannon$", "Ancestry diversity", x)
  x <- gsub("^cluster_family_12_degree_bin$", "Family history", x)
  x <- gsub("^cor_insulin$", "Insulin dose", x)
  x <- gsub("^imputed_C_pep$", "C-peptide", x)
  x <- gsub("^bmiz$", "BMI z-score", x)
  x <- gsub("^bmip$", "BMI percentile", x)
  x <- gsub("^bmi$", "BMI", x)
  x <- gsub("^Cluster_diabete_history$", "Diabetes duration/history", x)
  x <- gsub("^cluster_autoimmune_disease$", "Autoimmune disease", x)
  x
}

addam_domain <- function(x) {
  dplyr::case_when(
    x %in% c("X96GAD", "IA2", "ZnT8", "age_atb", "age_atb_per_year", "cluster_autoimmune_disease") ~ "Immune",
    grepl("GRS|^X[0-9]|_A[123]$", x) ~ "Genetic",
    x %in% c("bmiz", "bmip", "bmi") ~ "Anthropometry",
    x %in% c("age_BP", "age_BP_per_year", "Cluster_A1c", "cor_insulin", "imputed_C_pep", "duration_diabetes_per_year", "Cluster_diabete_history", "Cluster_sex_bin", "cluster_family_12_degree_bin", "cluster_shannon") ~ "Clinical",
    TRUE ~ "Other"
  )
}

cohort_n <- table1 |> filter(Characteristic == "N") |> pull(Value) |> as.character()
followup_median <- table1 |> filter(Characteristic == "Follow-up years, median [Q1, Q3]") |> pull(Value) |> as.character()
bp_days_median <- table1 |> filter(Characteristic == "Unique BP days, median [Q1, Q3]") |> pull(Value) |> as.character()
classifiable_n <- table7_summary$N_with_height_informed_peds_bp_classification[1]
peds_htn_n <- table7_summary$N_peds_htn_3plus_days[1]
peds_htn_pct <- table7_summary$pct_peds_htn_3plus_days[1]
height_matched_patients <- table8_qc |> filter(grepl("Patients with valid height within", Metric)) |> pull(Value)
height_matched_labs <- table8_qc |> filter(grepl("Creatinine labs with valid height within", Metric)) |> pull(Value)
fmt_count <- function(x) scales::comma(suppressWarnings(as.numeric(x)), accuracy = 1)

addam_label <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) {
  paste0(
    "ADDAM biobank\nn=", addam_qc$n_participants[1],
    "; elevated BP ", addam_qc$elevated_bp_n[1],
    " (", sprintf("%.1f", addam_qc$elevated_bp_percent[1]), "%)"
  )
} else {
  "ADDAM biobank\npending local results"
}

addam_auc <- if (!is.null(addam_ml) && nrow(addam_ml) > 0) {
  best <- addam_ml |> arrange(desc(auc)) |> slice(1)
  paste0(best$model[1], " LOOCV AUC ", sprintf("%.3f", best$auc[1]))
} else {
  "deep phenotype/genetics arm"
}

addam_feature_n <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) addam_qc$feature_columns_n[1] else NA_integer_
addam_snp_n <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) addam_qc$snp_dosage_columns_n[1] else NA_integer_
addam_prs_n <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) addam_qc$prs_available_n[1] else NA_integer_
addam_n <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) addam_qc$n_participants[1] else NA_integer_
addam_elevated_n <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) addam_qc$elevated_bp_n[1] else NA_integer_
addam_elevated_pct <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) addam_qc$elevated_bp_percent[1] else NA_real_
addam_best_auc <- if (!is.null(addam_ml) && nrow(addam_ml) > 0) {
  best <- addam_ml |> arrange(desc(auc)) |> slice(1)
  sprintf("%.3f", best$auc[1])
} else {
  "pending"
}

fig1_rapid_df <- model_df |>
  filter(!is.na(peds_htn_3plus_days), !is.na(egfr_slope_per_year)) |>
  mutate(group = ifelse(peds_htn_3plus_days == 1, "Sustained HTN", "No sustained HTN")) |>
  group_by(group) |>
  summarise(
    rapid_pct = mean(egfr_slope_per_year < -3, na.rm = TRUE) * 100,
    rapid_n = sum(egfr_slope_per_year < -3, na.rm = TRUE),
    denom = n(),
    .groups = "drop"
  ) |>
  mutate(
    group = factor(group, levels = c("No sustained HTN", "Sustained HTN")),
    y = ifelse(group == "No sustained HTN", 3.18, 2.88),
    x0 = 1.35,
    x1 = x0 + pmin(rapid_pct, 55) / 55 * 2.2,
    label = paste0(sprintf("%.1f", rapid_pct), "%"),
    note = paste0(rapid_n, "/", denom),
    color_hex = ifelse(group == "Sustained HTN", "#c43c39", "#2d82b7")
  )

fig1_imp_df <- if (!is.null(addam_imp) && nrow(addam_imp) > 0) {
  addam_imp |>
    slice_max(order_by = median_importance, n = 5) |>
    arrange(desc(median_importance)) |>
    mutate(
      feature_label = dplyr::case_when(
        grepl("^X[0-9].*_A[123]$", feature) ~ "Top SNP dosage",
        feature == "age_atb" ~ "Age at autoantibody",
        TRUE ~ pretty_addam_feature(feature)
      ),
      feature_label = ifelse(nchar(feature_label) > 27, paste0(substr(feature_label, 1, 25), "..."), feature_label),
      domain = addam_domain(feature),
      color_hex = dplyr::case_when(
        domain == "Clinical" ~ "#2d82b7",
        domain == "Immune" ~ "#c43c39",
        domain == "Genetic" ~ "#8c6bb1",
        domain == "Anthropometry" ~ "#4daf4a",
        TRUE ~ "#7f8c8d"
      ),
      y = seq(3.20, 2.60, length.out = n()),
      x0 = 7.62,
      x1 = x0 + median_importance / max(median_importance, na.rm = TRUE) * 1.85
    )
} else {
  tibble::tibble(feature_label = character(), domain = character(), color_hex = character(), y = numeric(), x0 = numeric(), x1 = numeric(), median_importance = numeric())
}

fig1_cards <- tibble::tibble(
  id = c("uddb", "addam"),
  xmin = c(0.45, 6.35),
  xmax = c(5.65, 11.55),
  ymin = c(1.85, 1.85),
  ymax = c(6.25, 6.25),
  fill = c("#fff5ec", "#edf7fd"),
  border = c("#d5602a", "#2d82b7"),
  title = c("Utah Diabetes Database (UDDB)", "ADDAM Biobank, Quebec"),
  subtitle = c("Population EHR: BP burden and kidney trajectories", "Deep phenotype: immune, clinical, and genetic risk"),
  accent = c("#d5602a", "#2d82b7")
)

fig1_tiles <- tibble::tibble(
  x = c(1.25, 2.95, 4.65, 7.05, 8.70, 10.35),
  y = c(4.70, 4.70, 4.70, 4.70, 4.70, 4.70),
  value = c(
    fmt_count(cohort_n),
    fmt_count(classifiable_n),
    paste0(fmt_count(peds_htn_n), "\n(", sprintf("%.1f", peds_htn_pct), "%)"),
    fmt_count(addam_n),
    paste0(fmt_count(addam_elevated_n), "\n(", sprintf("%.1f", addam_elevated_pct), "%)"),
    paste0("AUC\n", addam_best_auc)
  ),
  label = c(
    "pediatric T1D",
    "classifiable BP",
    ">=3 HTN encounter-days",
    "deeply phenotyped",
    "elevated BP",
    "best LOOCV"
  ),
  fill = c(rep("#fde7d4", 3), rep("#dceef9", 3)),
  border = c(rep("#e09463", 3), rep("#6ba9cf", 3))
)

fig1 <- ggplot() +
  geom_rect(
    data = fig1_cards,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, color = border),
    linewidth = 0.85
  ) +
  geom_rect(
    data = fig1_tiles,
    aes(xmin = x - 0.62, xmax = x + 0.62, ymin = y - 0.48, ymax = y + 0.48, fill = fill, color = border),
    linewidth = 0.4
  ) +
  geom_text(
    data = fig1_cards,
    aes(x = xmin + 0.32, y = ymax - 0.45, label = title),
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4.2,
    color = "#17202a"
  ) +
  geom_text(
    data = fig1_cards,
    aes(x = xmin + 0.32, y = ymax - 0.86, label = subtitle, color = accent),
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 2.8
  ) +
  geom_text(
    data = fig1_tiles,
    aes(x = x, y = y + 0.08, label = value),
    hjust = 0.5,
    vjust = 0.5,
    fontface = "bold",
    size = 3.55,
    lineheight = 0.85,
    color = "#17202a"
  ) +
  geom_text(
    data = fig1_tiles,
    aes(x = x, y = y - 0.32, label = label),
    hjust = 0.5,
    vjust = 0.5,
    size = 2.35,
    color = "#465760"
  ) +
  annotate("text", x = 0.80, y = 3.78, label = "Kidney signal", hjust = 0, fontface = "bold", size = 3.2, color = "#17202a") +
  annotate("text", x = 0.80, y = 3.50, label = "Rapid CKiD U25 eGFR decline (< -3/year)", hjust = 0, size = 2.55, color = "#465760") +
  geom_segment(data = fig1_rapid_df, aes(x = x0, xend = x1, y = y, yend = y), linewidth = 4.8, color = "#e6d3c4", lineend = "round") +
  geom_point(data = fig1_rapid_df, aes(x = x1, y = y, color = color_hex), size = 3.4) +
  geom_text(data = fig1_rapid_df, aes(x = x0 - 0.08, y = y, label = group), hjust = 1, vjust = 0.5, size = 2.45, color = "#37474f") +
  geom_text(data = fig1_rapid_df, aes(x = x1 + 0.18, y = y, label = paste(label, note, sep = "  ")), hjust = 0, vjust = 0.5, size = 2.45, color = "#37474f") +
  annotate(
    "text",
    x = 0.80,
    y = 2.28,
    label = paste0("CKiD U25 eGFR: ", fmt_count(height_matched_patients), " patients; ", fmt_count(height_matched_labs), " height-matched creatinine labs"),
    hjust = 0,
    size = 2.45,
    color = "#465760"
  ) +
  annotate(
    "text",
    x = 0.80,
    y = 2.05,
    label = "Historical UDDB database build retained as Supplemental Figure 1",
    hjust = 0,
    size = 2.25,
    color = "#607178"
  ) +
  annotate("text", x = 6.72, y = 3.78, label = "Biologic predictor signal", hjust = 0, fontface = "bold", size = 3.2, color = "#17202a") +
  annotate("text", x = 6.72, y = 3.50, label = "Top ADDAM feature-importance signals", hjust = 0, size = 2.55, color = "#465760") +
  geom_segment(data = fig1_imp_df, aes(x = x0, xend = x1, y = y, yend = y, color = color_hex), linewidth = 4.2, alpha = 0.82, lineend = "round") +
  geom_point(data = fig1_imp_df, aes(x = x1, y = y, color = color_hex), size = 2.8) +
  geom_text(data = fig1_imp_df, aes(x = x0 - 0.08, y = y, label = feature_label), hjust = 1, vjust = 0.5, size = 2.35, color = "#37474f") +
  annotate(
    "text",
    x = 6.72,
    y = 2.28,
    label = paste0(fmt_count(addam_feature_n), " features, ", fmt_count(addam_snp_n), " SNP dosage columns, GRS n=", fmt_count(addam_prs_n)),
    hjust = 0,
    size = 2.45,
    color = "#465760"
  ) +
  annotate(
    "text",
    x = 6.72,
    y = 2.05,
    label = "PRS/genetic information is contextual, not a stand-alone BP discriminator",
    hjust = 0,
    size = 2.25,
    color = "#607178"
  ) +
  annotate("rect", xmin = 0.85, xmax = 11.15, ymin = 0.35, ymax = 1.35, fill = "#223747", color = "#223747") +
  annotate("text", x = 1.20, y = 1.10, label = "Central inference", hjust = 0, fontface = "bold", size = 3.4, color = "white") +
  annotate(
    "text",
    x = 1.20,
    y = 0.70,
    label = "Elevated BP in pediatric T1D marks heterogeneous kidney-risk trajectories\nand an immune/age-related risk profile not captured by HbA1c, BMI, or PRS alone.",
    hjust = 0,
    size = 2.45,
    lineheight = 0.95,
    color = "white"
  ) +
  annotate("text", x = 6.0, y = 1.62, label = "one integrated paper", fontface = "bold", size = 2.75, color = "#607178") +
  scale_fill_identity() +
  scale_color_identity() +
  coord_cartesian(xlim = c(0, 12), ylim = c(0.15, 6.65), clip = "off") +
  labs(
    title = "Figure 1. Elevated BP in pediatric T1D links kidney-risk trajectories with immune/age signals",
    subtitle = "UDDB provides population EHR evidence; ADDAM tests whether deep phenotype and genetics explain elevated BP risk"
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "#4d5656", margin = margin(b = 12)),
    plot.margin = margin(16, 18, 12, 18)
  )

classification_long <- table2_class |>
  mutate(group_label = paste(SEX_CD, age_band_method)) |>
  select(group_label, method_label, normal_days, elevated_days, stage1_days, stage2_days) |>
  tidyr::pivot_longer(
    cols = c(normal_days, elevated_days, stage1_days, stage2_days),
    names_to = "category",
    values_to = "days"
  ) |>
  mutate(
    category = factor(
      category,
      levels = c("normal_days", "elevated_days", "stage1_days", "stage2_days"),
      labels = c("Normal", "Elevated", "Stage 1 HTN", "Stage 2 HTN")
    ),
    group_label = factor(group_label, levels = unique(table2_class |> mutate(group_label = paste(SEX_CD, age_band_method)) |> pull(group_label)))
  )

fig2 <- ggplot(classification_long, aes(x = group_label, y = days, fill = category)) +
  geom_col(position = "fill", width = 0.72, color = "white", linewidth = 0.2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(values = c("Normal" = "#7fbf7b", "Elevated" = "#fdd66c", "Stage 1 HTN" = "#f28e5c", "Stage 2 HTN" = "#c43c39")) +
  labs(
    title = "Figure 2. Height-informed pediatric BP classification",
    subtitle = "Age/sex/height percentiles for ages 1-12 years; fixed adolescent thresholds for ages 13-17 years",
    x = "Sex and age band",
    y = "Proportion of classifiable BP days",
    fill = "BP category"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.text.x = element_text(angle = 35, hjust = 1),
    legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

traj_by_peds_htn <- creat_labs_u25 |>
  left_join(model_df |> select(Patient_ID, peds_htn_3plus_days), by = "Patient_ID") |>
  mutate(
    peds_htn_group = case_when(
      peds_htn_3plus_days == 1 ~ "Pediatric HTN >=3 BP days",
      peds_htn_3plus_days == 0 ~ "No pediatric HTN",
      TRUE ~ NA_character_
    ),
    year_bin = floor(years_since_dx),
    eGFR_u25_cap = pmin(eGFR_u25, 150)
  ) |>
  filter(!is.na(peds_htn_group), is.finite(year_bin), year_bin >= 0, year_bin <= 15)

traj_summary <- traj_by_peds_htn |>
  group_by(peds_htn_group, year_bin) |>
  summarise(
    eGFR_median = median(eGFR_u25_cap, na.rm = TRUE),
    eGFR_Q1 = quantile(eGFR_u25_cap, 0.25, na.rm = TRUE),
    eGFR_Q3 = quantile(eGFR_u25_cap, 0.75, na.rm = TRUE),
    n_patients = n_distinct(Patient_ID),
    .groups = "drop"
  ) |>
  filter(n_patients >= 5)

fig3 <- ggplot(traj_summary, aes(x = year_bin, y = eGFR_median, color = peds_htn_group, fill = peds_htn_group)) +
  geom_ribbon(aes(ymin = eGFR_Q1, ymax = eGFR_Q3), alpha = 0.16, color = NA) +
  geom_line(linewidth = 1.1) +
  geom_point(aes(size = n_patients), alpha = 0.85) +
  geom_hline(yintercept = 90, linetype = "dashed", color = "#555555") +
  scale_fill_manual(values = c("No pediatric HTN" = "#2878b5", "Pediatric HTN >=3 BP days" = "#c43c39")) +
  scale_color_manual(values = c("No pediatric HTN" = "#2878b5", "Pediatric HTN >=3 BP days" = "#c43c39")) +
  scale_size_continuous(range = c(1.5, 4.5)) +
  coord_cartesian(ylim = c(40, 150)) +
  labs(
    title = "Figure 3. CKiD U25 eGFR trajectories by pediatric BP phenotype",
    subtitle = "Height-based CKiD U25 eGFR; points show supported year bins with >=5 patients",
    x = "Years since T1D diagnosis",
    y = "Median eGFR (IQR), capped at 150",
    color = "",
    fill = "",
    size = "Patients"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    legend.position = "bottom"
  )

slope_df <- model_df |>
  filter(!is.na(peds_htn_3plus_days), !is.na(egfr_slope_per_year)) |>
  mutate(
    peds_htn_group = factor(
      ifelse(peds_htn_3plus_days == 1, "Pediatric HTN >=3 BP days", "No pediatric HTN"),
      levels = c("No pediatric HTN", "Pediatric HTN >=3 BP days")
    ),
    rapid_decline = egfr_slope_per_year < -3
  )

rapid_summary <- slope_df |>
  group_by(peds_htn_group) |>
  summarise(
    rapid_decline_pct = mean(rapid_decline, na.rm = TRUE) * 100,
    rapid_decline_n = sum(rapid_decline, na.rm = TRUE),
    slope_n = n(),
    label = paste0(sprintf("%.1f", rapid_decline_pct), "%\n", rapid_decline_n, "/", slope_n),
    .groups = "drop"
  )

fig4a <- ggplot(rapid_summary, aes(x = peds_htn_group, y = rapid_decline_pct, fill = peds_htn_group)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = label), vjust = -0.2, size = 3.5, lineheight = 0.9) +
  scale_fill_manual(values = c("No pediatric HTN" = "#2878b5", "Pediatric HTN >=3 BP days" = "#c43c39")) +
  coord_cartesian(ylim = c(0, max(rapid_summary$rapid_decline_pct, na.rm = TRUE) * 1.25)) +
  labs(
    title = "A. Rapid eGFR decline",
    x = "",
    y = "Patients with slope < -3 (%)",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

fig4b <- ggplot(slope_df, aes(x = peds_htn_group, y = egfr_slope_per_year, fill = peds_htn_group)) +
  geom_hline(yintercept = -3, linetype = "dashed", color = "#555555") +
  geom_hline(yintercept = 0, linetype = "solid", color = "#999999", linewidth = 0.3) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, alpha = 0.25, size = 1.2) +
  scale_fill_manual(values = c("No pediatric HTN" = "#2878b5", "Pediatric HTN >=3 BP days" = "#c43c39")) +
  labs(
    title = "B. CKiD U25 eGFR slope distribution",
    x = "",
    y = "eGFR slope (mL/min/1.73m2/year)",
    subtitle = "Dashed line = rapid decline threshold (< -3)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  )

fig4 <- (fig4a | fig4b) +
  plot_annotation(
    title = "Figure 4. Height-informed pediatric BP phenotype and kidney decline",
    subtitle = "Phenotype defined as stage 1 or stage 2 hypertension on >=3 unique pediatric BP days (encounter dates) after same-day readings were collapsed",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10.5, color = "#4d5656")
    )
  )

if (!is.null(addam_ml) && !is.null(addam_inc) && !is.null(addam_imp) && !is.null(addam_univ)) {
  ml_plot_df <- addam_ml |>
    mutate(
      model = forcats::fct_reorder(model, auc),
      auc_label = sprintf("%.3f", auc)
    )

  fig5a_left <- ggplot(ml_plot_df, aes(x = auc, y = model)) +
    geom_segment(aes(x = 0.5, xend = auc, yend = model), linewidth = 1.1, color = "#d7dde5") +
    geom_point(size = 3.4, color = "#2c7fb8") +
    geom_text(aes(label = auc_label), nudge_x = 0.012, size = 3.4, hjust = 0) +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "#7f8c8d") +
    coord_cartesian(xlim = c(0.48, 0.75), clip = "off") +
    labs(
      title = "A. Model discrimination across families",
      x = "AUC",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  inc_plot_df <- addam_inc |>
    mutate(
      model_set = dplyr::recode(
        model_set,
        "Clinical core" = "Clinical core",
        "Clinical + autoantibodies" = "Clinical +\nautoantibodies",
        "Clinical + autoantibodies + PRS" = "Clinical + autoantibodies\n+ PRS",
        "Clinical + autoantibodies + ancestry" = "Clinical + autoantibodies\n+ ancestry",
        "Clinical + autoantibodies + ancestry + PRS" = "Clinical + autoantibodies\n+ ancestry + PRS"
      ),
      model_set = factor(
        model_set,
        levels = c(
          "Clinical core",
          "Clinical +\nautoantibodies",
          "Clinical + autoantibodies\n+ PRS",
          "Clinical + autoantibodies\n+ ancestry",
          "Clinical + autoantibodies\n+ ancestry + PRS"
        )
      ),
      delta_label = sprintf("%+.003f", delta_auc_vs_clinical)
    )

  fig5a_right <- ggplot(inc_plot_df, aes(x = delta_auc_vs_clinical, y = model_set)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "#7f8c8d") +
    geom_segment(aes(x = 0, xend = delta_auc_vs_clinical, yend = model_set), linewidth = 1.1, color = "#e7d8b9") +
    geom_point(size = 3.4, color = "#d17c2f") +
    geom_text(aes(label = delta_label), nudge_x = 0.004, size = 3.2, hjust = 0) +
    coord_cartesian(xlim = c(min(-0.03, min(inc_plot_df$delta_auc_vs_clinical, na.rm = TRUE) - 0.005), 0.02), clip = "off") +
    labs(
      title = "B. Added value beyond the clinical core",
      x = "Delta AUC vs clinical core",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank()
    )

  fig5a <- (fig5a_left / fig5a_right) +
    plot_annotation(
      title = "Figure 5A. ADDAM model performance within the integrated paper",
      subtitle = "ADDAM extends the UDDB clinical story by testing whether deeper phenotyping improves elevated-BP discrimination",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10.5, color = "#4d5656")
      )
    )

  imp_plot_df <- addam_imp |>
    slice_max(order_by = median_importance, n = 10) |>
    mutate(
      feature_label = pretty_addam_feature(feature),
      domain = addam_domain(feature),
      feature_label = forcats::fct_reorder(feature_label, median_importance)
    )

  fig5b_left <- ggplot(imp_plot_df, aes(x = median_importance, y = feature_label, color = domain)) +
    geom_segment(aes(x = ci_2_5, xend = ci_97_5, yend = feature_label), linewidth = 1.1, alpha = 0.9) +
    geom_point(size = 3) +
    scale_color_manual(values = c(Clinical = "#2c7fb8", Immune = "#c43c39", Anthropometry = "#4daf4a", Genetic = "#8c6bb1", Other = "#7f8c8d")) +
    labs(
      title = "A. Bootstrap feature-importance profile",
      x = "Median bootstrap importance (95% interval)",
      y = NULL,
      color = "Domain"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom",
      panel.grid.major.y = element_blank()
    )

  univ_keep <- c("age_BP_per_year", "IA2", "age_atb_per_year", "GRS", "bmiz", "ZnT8", "X96GAD", "duration_diabetes_per_year")
  univ_plot_df <- addam_univ |>
    filter(feature %in% univ_keep, !is.na(odds_ratio), !is.na(ci_low), !is.na(ci_high)) |>
    mutate(
      feature_label = pretty_addam_feature(feature),
      domain = addam_domain(feature),
      feature_label = forcats::fct_reorder(feature_label, odds_ratio)
    )

  fig5b_right <- ggplot(univ_plot_df, aes(x = odds_ratio, y = feature_label, color = domain)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "#7f8c8d") +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = feature_label), linewidth = 1.1, alpha = 0.9) +
    geom_point(size = 3) +
    scale_x_log10() +
    scale_color_manual(values = c(Clinical = "#2c7fb8", Immune = "#c43c39", Anthropometry = "#4daf4a", Genetic = "#8c6bb1", Other = "#7f8c8d")) +
    labs(
      title = "B. Selected univariate associations with elevated BP",
      x = "Odds ratio (log scale)",
      y = NULL,
      color = "Domain"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      panel.grid.major.y = element_blank()
    )

  fig5b <- (fig5b_left / fig5b_right) +
    plot_annotation(
      title = "Figure 5B. ADDAM biologic and genetic signal profile",
      subtitle = "Immune and age-related signals dominate; BMI and T1D GRS/PRS are contextual rather than stand-alone discriminators",
      theme = theme(
        plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10.5, color = "#4d5656")
      )
    )

  save_paper_plot(fig5a, "Figure5A_ADDAM_ROC_ModelPerformance.jpeg", width = 9.5, height = 8.2)
  save_paper_plot(fig5b, "Figure5B_ADDAM_FeatureImportance.jpeg", width = 10, height = 9.6)
}

save_paper_plot(fig1, "Figure1_Integrated_UDDB_ADDAM_Cohorts.png", width = 12, height = 5.8)
save_paper_plot(fig2, "Figure2_HeightInformed_Pediatric_BP_Classification.png", width = 8.5, height = 5.2)
save_paper_plot(fig3, "Figure3_CKiD_U25_eGFR_by_Pediatric_BP_Phenotype.png", width = 8, height = 5.2)
save_paper_plot(fig4, "Figure4_ClinicalBridge_PediatricBP_AdultKidney.png", width = 10, height = 5.2)

addam_fig_dir <- file.path(dirs$results_dir, "addam_genetics", "figures")
copy_if_exists(file.path(addam_fig_dir, "7_8_PRC.jpeg"), supplemental_figures_dir, "Supplemental_ADDAM_PrecisionRecall.jpeg")
copy_if_exists(file.path(addam_fig_dir, "7_8_PRC_band_Xgboost.jpeg"), supplemental_figures_dir, "Supplemental_ADDAM_XGBoost_PR_Band.jpeg")
copy_if_exists(file.path(addam_fig_dir, "7_8_PRC_band_Random Forestx.jpeg"), supplemental_figures_dir, "Supplemental_ADDAM_RF_PR_Band.jpeg")

supplemental_map <- c(
  "Figure1_AdultSlopeDensity.png",
  "Figure2_AdultSlopeBoxplot.png",
  "Figure3_eGFR_Trajectories_Subset.png",
  "Figure5A_spaghetti.png",
  "Figure5B_median_IQR.png",
  "Figure5C_trajectory_categories.png",
  "Figure6A_BP_spaghetti.png",
  "Figure6B_BP_median.png",
  "Figure6C_BP_trajectory.png",
  "Figure7A_AdultHTN_by_HeightInformedPedsHTN.png",
  "Figure7B_eGFRSlope_by_HeightInformedPedsHTN.png"
)

for (fname in supplemental_map) {
  copy_if_exists(file.path(dirs$figures_dir, fname), supplemental_figures_dir)
}

legacy_omit <- c(
  "Figure7A_AdultHTN_by_PedsPercentileHTN.png",
  "Figure7B_eGFRSlope_by_PedsPercentileHTN.png"
)

catalog <- c(
  "# Figure Catalog and Paper Story",
  "",
  "## Main Paper Figures",
  "",
  "1. `Figure1_Integrated_UDDB_ADDAM_Cohorts.png` - central graphical summary showing why UDDB and ADDAM are one paper: population kidney-risk signal plus deep phenotype/genetic risk profiling.",
  "2. `Figure2_HeightInformed_Pediatric_BP_Classification.png` - primary pediatric BP phenotype figure. This replaces the raw median/IQR BP table as the clinical BP result.",
  "3. `Figure3_CKiD_U25_eGFR_by_Pediatric_BP_Phenotype.png` - height-based CKiD U25 kidney trajectories by height-informed pediatric BP phenotype.",
  "4. `Figure4_ClinicalBridge_PediatricBP_AdultKidney.png` - links the height-informed pediatric phenotype to patient-level kidney decline. Here, `>=3 BP days` means >=3 unique encounter dates/patient-days after same-day measurements were collapsed.",
  "5. `Figure5A_ADDAM_ROC_ModelPerformance.jpeg` - ADDAM discrimination and incremental-value figure, positioned as a complementary extension of the same paper rather than a separate study.",
  "6. `Figure5B_ADDAM_FeatureImportance.jpeg` - ADDAM biologic-signal figure showing which immune, clinical, anthropometric, and genetic variables carry the strongest signal.",
  "",
  "## Supplemental Figures",
  "",
  "Use `results/supplemental_figures/` for adult-only BP/slope plots, raw eGFR spaghetti/median trajectories, SBP-tertile sensitivity plots, and ADDAM PR curves.",
  "Keep the historical UDDB database construction diagram as Supplemental Figure 1 rather than repeating those database-freeze details in the main figure set.",
  "",
  "## Figures Removed From the Main Story",
  "",
  paste0("- `", legacy_omit, "` - older duplicate labels using the percentile wording; replaced by height-informed terminology."),
  "- SBP tertile figures should not lead the paper because height-informed pediatric BP classification is more clinically defensible.",
  "- Raw BP median/IQR plots/tables are descriptive and should be supplemental unless a reviewer asks for them in the main text.",
  "",
  "## Suggested Narrative",
  "",
  "Lead with height-informed BP classification rather than raw BP distributions. Then connect that phenotype to height-based CKiD U25 eGFR trajectories and slope-defined kidney decline. Close the paper with ADDAM as a mechanistic and predictive extension of the same central story, not as a separate paper."
)

writeLines(catalog, file.path(dirs$results_dir, "Figure_Catalog.md"), useBytes = TRUE)
writeLines(catalog, file.path(dirs$publication_dir, "Figure_Strategy.md"), useBytes = TRUE)
message("Curated paper figures in: ", paper_figures_dir)
message("Curated supplemental figures in: ", supplemental_figures_dir)
