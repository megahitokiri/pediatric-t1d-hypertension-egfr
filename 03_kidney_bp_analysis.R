script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("dplyr", "readr", "ggplot2"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

dirs <- ensure_pipeline_dirs(script_dir)

kidney_bp_analysis <- load_derived_rds("kidney_bp_analysis", dirs) |>
  mutate(
    HTN_group = case_when(
      is.na(HTN_2plus_days) ~ "No adult BP phenotype",
      HTN_2plus_days ~ "Adult HTN (>=2 days)",
      TRUE ~ "No adult HTN"
    ),
    adultBP_available = !is.na(n_BP_days_adult) & n_BP_days_adult > 0
  )

creat_labs_u25 <- load_derived_rds("creat_labs_u25", dirs)

table3_all <- kidney_bp_analysis |>
  group_by(HTN_group) |>
  summarise(
    N = n(),
    N_with_creat = sum(!is.na(last_eGFR)),
    N_with_slope = sum(!is.na(egfr_slope_per_year)),
    last_eGFR = fmt_med_q(last_eGFR),
    slope = fmt_med_q(egfr_slope_per_year),
    pct_any_eGFR_lt90 = round(mean(any_eGFR_lt90, na.rm = TRUE) * 100, 1),
    pct_any_eGFR_lt60 = round(mean(any_eGFR_lt60, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

table3_adultBPonly <- kidney_bp_analysis |>
  filter(adultBP_available) |>
  mutate(HTN_group2 = ifelse(HTN_2plus_days, "Adult HTN (>=2 days)", "No adult HTN")) |>
  group_by(HTN_group2) |>
  summarise(
    N = n(),
    N_with_creat = sum(!is.na(last_eGFR)),
    N_with_slope = sum(!is.na(egfr_slope_per_year)),
    last_eGFR = fmt_med_q(last_eGFR),
    slope = fmt_med_q(egfr_slope_per_year),
    pct_any_eGFR_lt90 = round(mean(any_eGFR_lt90, na.rm = TRUE) * 100, 1),
    pct_any_eGFR_lt60 = round(mean(any_eGFR_lt60, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

slope_df <- kidney_bp_analysis |>
  filter(!is.na(egfr_slope_per_year)) |>
  mutate(
    age_group = ifelse(AgeAtLastVisit < 18, "Pediatric (<18y)", "Adult (>=18y)"),
    HTN_group_simple = case_when(
      is.na(HTN_2plus_days) ~ "No adult BP phenotype",
      HTN_2plus_days ~ "Adult HTN",
      TRUE ~ "No adult HTN"
    )
  )

table4A <- slope_df |>
  group_by(age_group) |>
  summarise(
    N = n(),
    slope = fmt_med_q(egfr_slope_per_year),
    pct_decline_gt3 = round(mean(egfr_slope_per_year < -3, na.rm = TRUE) * 100, 1),
    pct_decline_gt5 = round(mean(egfr_slope_per_year < -5, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

adult_slope_test <- slope_df |>
  filter(age_group == "Adult (>=18y)", HTN_group_simple %in% c("Adult HTN", "No adult HTN"))

table4B <- adult_slope_test |>
  group_by(HTN_group_simple) |>
  summarise(
    N = n(),
    slope = fmt_med_q(egfr_slope_per_year),
    pct_decline_gt3 = round(mean(egfr_slope_per_year < -3, na.rm = TRUE) * 100, 1),
    pct_decline_gt5 = round(mean(egfr_slope_per_year < -5, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )

table4C <- if (nrow(adult_slope_test) > 1 && length(unique(adult_slope_test$HTN_group_simple)) == 2) {
  wilcox_slope <- stats::wilcox.test(egfr_slope_per_year ~ HTN_group_simple, data = adult_slope_test, exact = FALSE)
  rapid_decline <- table(adult_slope_test$HTN_group_simple, adult_slope_test$egfr_slope_per_year < -5)
  fisher_rapid <- stats::fisher.test(rapid_decline)
  tibble::tibble(
    Test = c("Wilcoxon slope (Adult HTN vs No adult HTN)", "Fisher rapid decline (>5)"),
    p_value = c(wilcox_slope$p.value, fisher_rapid$p.value)
  )
} else {
  tibble::tibble(Test = "Adult slope comparison", p_value = NA_real_)
}

write_table_csv(table3_all, "Table3_all.csv", dirs)
write_table_csv(table3_adultBPonly, "Table3_adultBPonly.csv", dirs)
write_table_csv(table4A, "Table4A_slopes_by_age.csv", dirs)
write_table_csv(table4B, "Table4B_adult_slopes_by_HTN.csv", dirs)
write_table_csv(table4C, "Table4C_tests.csv", dirs)
write_table_csv(table3_all, "Table3_Kidney_by_HTN.csv", dirs)
write_table_csv(table3_adultBPonly, "Table3_Kidney_by_HTN_adultBPonly.csv", dirs)

if (nrow(adult_slope_test) > 0) {
  p1 <- ggplot(adult_slope_test, aes(x = egfr_slope_per_year, fill = HTN_group_simple)) +
    geom_density(alpha = 0.35) +
    geom_vline(xintercept = -3, linetype = "dashed") +
    geom_vline(xintercept = -5, linetype = "dotted") +
    labs(
      title = "Adult eGFR slopes by hypertension phenotype",
      x = "eGFR slope (mL/min/1.73m2/year)",
      y = "Density",
      fill = "Group"
    ) +
    theme_minimal()

  p2 <- ggplot(adult_slope_test, aes(x = HTN_group_simple, y = egfr_slope_per_year)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.25) +
    geom_hline(yintercept = -3, linetype = "dashed") +
    geom_hline(yintercept = -5, linetype = "dotted") +
    labs(
      title = "Adult eGFR slopes (distribution) by hypertension phenotype",
      x = "",
      y = "eGFR slope (mL/min/1.73m2/year)"
    ) +
    theme_minimal()

  traj_df <- creat_labs_u25 |>
    left_join(
      kidney_bp_analysis |> select(Patient_ID, HTN_group, adultBP_available),
      by = "Patient_ID"
    ) |>
    filter(adultBP_available)

  set.seed(1)
  example_ids <- traj_df |>
    distinct(Patient_ID, HTN_group) |>
    group_by(HTN_group) |>
    slice_sample(prop = 1) |>
    slice_head(n = 60) |>
    ungroup() |>
    pull(Patient_ID)

  traj_small <- traj_df |> filter(Patient_ID %in% example_ids)

  p3 <- ggplot(traj_small, aes(x = RESULT_DATE, y = eGFR_u25, group = Patient_ID, color = HTN_group)) +
    geom_line(alpha = 0.15) +
    geom_smooth(aes(group = HTN_group), method = "loess", se = TRUE) +
    labs(
      title = "eGFR trajectories (subset) by hypertension phenotype",
      x = "Date",
      y = "eGFR (CKiD U25, height-based creatinine equation)",
      color = "Group"
    ) +
    theme_minimal()

  save_plot(p1, "Figure1_AdultSlopeDensity.png", dirs, width = 7, height = 4.5)
  save_plot(p2, "Figure2_AdultSlopeBoxplot.png", dirs, width = 7, height = 4.5)
  save_plot(p3, "Figure3_eGFR_Trajectories_Subset.png", dirs, width = 7.5, height = 5)
}

message("Saved kidney/BP analysis outputs to: ", dirs$results_dir)
