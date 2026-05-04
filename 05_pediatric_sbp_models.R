script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("dplyr", "readr", "ggplot2", "splines"))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(splines))

dirs <- ensure_pipeline_dirs(script_dir)

analysis_cohort <- load_derived_rds("analysis_cohort", dirs)
bp_unique_day_analysis <- load_derived_rds("bp_unique_day_analysis", dirs)
traj_long <- load_derived_rds("traj_long", dirs)
table5_patientlevel <- load_derived_rds("table5_patientlevel", dirs)

peds_bp <- bp_unique_day_analysis

if (!"SEX_CD" %in% names(peds_bp)) {
  peds_bp <- peds_bp |>
    left_join(analysis_cohort |> select(Patient_ID, SEX_CD), by = "Patient_ID")
}

peds_bp <- peds_bp |>
  filter(
    is.finite(age_at_bp),
    age_at_bp >= 0,
    age_at_bp < 18,
    is.finite(SBP),
    SBP >= 60,
    SBP <= 220
  ) |>
  filter(!is.na(SEX_CD)) |>
  transmute(
    Patient_ID,
    bp_date,
    SEX_CD,
    age = age_at_bp,
    SBP,
    DBP
  )

fit_sbp <- stats::lm(SBP ~ splines::ns(age, df = 4) + SEX_CD, data = peds_bp)

peds_patient_score <- peds_bp |>
  mutate(SBP_resid = SBP - stats::predict(fit_sbp, newdata = peds_bp)) |>
  group_by(Patient_ID) |>
  summarise(
    n_peds_bp_days = n(),
    Mean_SBP_peds = mean(SBP, na.rm = TRUE),
    Mean_DBP_peds = mean(DBP, na.rm = TRUE),
    SBP_resid_median = median(SBP_resid, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(n_peds_bp_days >= 3)

q33 <- stats::quantile(peds_patient_score$SBP_resid_median, probs = 1 / 3, na.rm = TRUE)
q67 <- stats::quantile(peds_patient_score$SBP_resid_median, probs = 2 / 3, na.rm = TRUE)

peds_patient_score <- peds_patient_score |>
  mutate(
    SBP_tertile = case_when(
      SBP_resid_median <= q33 ~ "Low SBP",
      SBP_resid_median <= q67 ~ "Mid SBP",
      TRUE ~ "High SBP"
    ),
    SBP_tertile = factor(SBP_tertile, levels = c("Low SBP", "Mid SBP", "High SBP"))
  )

traj6 <- traj_long |>
  left_join(
    peds_patient_score |> select(Patient_ID, n_peds_bp_days, Mean_SBP_peds, Mean_DBP_peds, SBP_resid_median, SBP_tertile),
    by = "Patient_ID"
  ) |>
  left_join(table5_patientlevel |> select(Patient_ID, slope), by = "Patient_ID") |>
  left_join(analysis_cohort |> select(Patient_ID, AgeAtFirstT1), by = "Patient_ID") |>
  mutate(
    years_since_dx = as.numeric(years_since_dx),
    eGFR_u25 = as.numeric(eGFR_u25),
    eGFR_u25_cap = pmin(eGFR_u25, 200),
    SBP_tertile = factor(SBP_tertile, levels = c("Low SBP", "Mid SBP", "High SBP"))
  ) |>
  filter(
    !is.na(SBP_tertile),
    is.finite(years_since_dx),
    years_since_dx >= 0,
    is.finite(eGFR_u25),
    eGFR_u25 > 0
  )

table6_bp_by_trajectory <- traj6 |>
  distinct(Patient_ID, AgeAtFirstT1, slope, Mean_SBP_peds, Mean_DBP_peds, n_peds_bp_days, SBP_tertile) |>
  mutate(Trajectory = trajectory_category(slope)) |>
  group_by(Trajectory) |>
  summarise(
    N = n(),
    AgeDx_median = median(AgeAtFirstT1, na.rm = TRUE),
    AgeDx_Q1 = quantile(AgeAtFirstT1, 0.25, na.rm = TRUE),
    AgeDx_Q3 = quantile(AgeAtFirstT1, 0.75, na.rm = TRUE),
    SBP_median = median(Mean_SBP_peds, na.rm = TRUE),
    SBP_Q1 = quantile(Mean_SBP_peds, 0.25, na.rm = TRUE),
    SBP_Q3 = quantile(Mean_SBP_peds, 0.75, na.rm = TRUE),
    DBP_median = median(Mean_DBP_peds, na.rm = TRUE),
    DBP_Q1 = quantile(Mean_DBP_peds, 0.25, na.rm = TRUE),
    DBP_Q3 = quantile(Mean_DBP_peds, 0.75, na.rm = TRUE),
    BPdays_median = median(n_peds_bp_days, na.rm = TRUE),
    BPdays_Q1 = quantile(n_peds_bp_days, 0.25, na.rm = TRUE),
    BPdays_Q3 = quantile(n_peds_bp_days, 0.75, na.rm = TRUE),
    slope_median = median(slope, na.rm = TRUE),
    slope_Q1 = quantile(slope, 0.25, na.rm = TRUE),
    slope_Q3 = quantile(slope, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

write_table_csv(table6_bp_by_trajectory, "Table6_BP_by_Trajectory.csv", dirs)

p6a_simple <- ggplot(traj6, aes(x = years_since_dx, y = eGFR_u25_cap, group = Patient_ID, color = SBP_tertile)) +
  geom_line(alpha = 0.15) +
  labs(
    title = "Figure 6A. eGFR trajectories by pediatric SBP tertile",
    x = "Years since T1D diagnosis",
    y = "eGFR (mL/min/1.73m2)",
    color = "Pediatric SBP tertile"
  ) +
  theme_minimal()

traj_bp_summary <- traj6 |>
  mutate(year_bin = floor(years_since_dx)) |>
  group_by(SBP_tertile, year_bin) |>
  summarise(
    eGFR_median = median(eGFR_u25_cap, na.rm = TRUE),
    eGFR_Q1 = quantile(eGFR_u25_cap, 0.25, na.rm = TRUE),
    eGFR_Q3 = quantile(eGFR_u25_cap, 0.75, na.rm = TRUE),
    n_patients = n_distinct(Patient_ID),
    .groups = "drop"
  ) |>
  filter(n_patients >= 10)

p6b_simple <- ggplot(traj_bp_summary, aes(x = year_bin, y = eGFR_median, color = SBP_tertile)) +
  geom_ribbon(aes(ymin = eGFR_Q1, ymax = eGFR_Q3, fill = SBP_tertile), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1.0) +
  labs(
    title = "Figure 6B. Median eGFR trajectory by pediatric SBP tertile (IQR)",
    x = "Years since diagnosis (binned)",
    y = "Median eGFR (IQR)",
    color = "Pediatric SBP tertile",
    fill = "Pediatric SBP tertile"
  ) +
  theme_minimal()

p6c_simple <- traj6 |>
  distinct(Patient_ID, slope, SBP_tertile) |>
  mutate(Trajectory = trajectory_category(slope)) |>
  ggplot(aes(x = Trajectory, fill = SBP_tertile)) +
  geom_bar(position = "fill") +
  labs(
    title = "Figure 6C. Kidney trajectory class by pediatric SBP tertile",
    x = "Trajectory category",
    y = "Proportion of patients",
    fill = "Pediatric SBP tertile"
  ) +
  theme_minimal()

save_plot(p6a_simple, "Figure6A_BP_spaghetti.png", dirs, width = 7.5, height = 5)
save_plot(p6b_simple, "Figure6B_BP_median.png", dirs, width = 7.5, height = 5)
save_plot(p6c_simple, "Figure6C_BP_trajectory.png", dirs, width = 7.5, height = 5)

if (
  requireNamespace("lme4", quietly = TRUE) &&
  requireNamespace("lmerTest", quietly = TRUE) &&
  requireNamespace("broom.mixed", quietly = TRUE)
) {
  traj6_model <- traj6 |>
    mutate(t_center = years_since_dx - mean(years_since_dx, na.rm = TRUE))

  m1 <- lmerTest::lmer(
    eGFR_u25 ~ t_center * SBP_tertile + (1 + t_center | Patient_ID),
    data = traj6_model,
    REML = TRUE
  )

  table6_mixed <- broom.mixed::tidy(m1, effects = "fixed", conf.int = TRUE) |>
    mutate(
      estimate = round(estimate, 3),
      std.error = round(std.error, 3),
      statistic = round(statistic, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3),
      p.value = signif(p.value, 3)
    )

  m0 <- lmerTest::lmer(
    eGFR_u25 ~ t_center + SBP_tertile + (1 + t_center | Patient_ID),
    data = traj6_model,
    REML = FALSE
  )

  m1_ml <- lmerTest::lmer(
    eGFR_u25 ~ t_center * SBP_tertile + (1 + t_center | Patient_ID),
    data = traj6_model,
    REML = FALSE
  )

  interaction_lrt <- as.data.frame(stats::anova(m0, m1_ml))

  table6_tertile_counts <- traj6_model |>
    group_by(Patient_ID, SBP_tertile) |>
    summarise(
      n_egfr_points = n(),
      followup_years = max(years_since_dx, na.rm = TRUE),
      .groups = "drop"
    ) |>
    group_by(SBP_tertile) |>
    summarise(
      N_patients = n_distinct(Patient_ID),
      egfr_points_median = median(n_egfr_points, na.rm = TRUE),
      egfr_points_Q1 = quantile(n_egfr_points, 0.25, na.rm = TRUE),
      egfr_points_Q3 = quantile(n_egfr_points, 0.75, na.rm = TRUE),
      followup_median = median(followup_years, na.rm = TRUE),
      followup_Q1 = quantile(followup_years, 0.25, na.rm = TRUE),
      followup_Q3 = quantile(followup_years, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  table6_trajectory_summary <- traj6_model |>
    mutate(year_bin = floor(years_since_dx)) |>
    group_by(SBP_tertile, year_bin) |>
    summarise(
      eGFR_median = median(eGFR_u25, na.rm = TRUE),
      eGFR_Q1 = quantile(eGFR_u25, 0.25, na.rm = TRUE),
      eGFR_Q3 = quantile(eGFR_u25, 0.75, na.rm = TRUE),
      n_points = n(),
      n_patients = n_distinct(Patient_ID),
      .groups = "drop"
    )

  write_table_csv(table6_tertile_counts, "Table6_TertileCounts.csv", dirs)
  write_table_csv(table6_mixed, "Table6_MixedEffects.csv", dirs)
  write_table_csv(interaction_lrt, "Table6_Interaction_LRT.csv", dirs)
  write_table_csv(table6_trajectory_summary, "Table6_TertileTrajectorySummary.csv", dirs)

  set.seed(1)
  sampled_ids <- traj6_model |>
    distinct(Patient_ID, SBP_tertile) |>
    group_by(SBP_tertile) |>
    slice_sample(prop = 1) |>
    slice_head(n = 60) |>
    ungroup() |>
    pull(Patient_ID)

  traj6_sample <- traj6_model |> filter(Patient_ID %in% sampled_ids)

  p6a_model <- ggplot(traj6_sample, aes(x = years_since_dx, y = eGFR_u25_cap, group = Patient_ID, linetype = SBP_tertile)) +
    geom_line(alpha = 0.25) +
    coord_cartesian(ylim = c(0, 200)) +
    labs(
      title = "Figure 6A. eGFR trajectories (sampled) by pediatric SBP tertile",
      x = "Years since T1D diagnosis",
      y = "eGFR (CKiD U25, capped at 200)",
      linetype = "Pediatric SBP tertile"
    ) +
    theme_minimal()

  p6b_model <- traj6_model |>
    mutate(year_bin = floor(years_since_dx)) |>
    group_by(SBP_tertile, year_bin) |>
    summarise(
      eGFR_median = median(eGFR_u25_cap, na.rm = TRUE),
      eGFR_Q1 = quantile(eGFR_u25_cap, 0.25, na.rm = TRUE),
      eGFR_Q3 = quantile(eGFR_u25_cap, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    ggplot(aes(x = year_bin, y = eGFR_median)) +
    geom_ribbon(aes(ymin = eGFR_Q1, ymax = eGFR_Q3), alpha = 0.2) +
    geom_line(linewidth = 1.0) +
    facet_wrap(~ SBP_tertile, nrow = 1) +
    coord_cartesian(ylim = c(0, 200)) +
    labs(
      title = "Figure 6B. Median eGFR trajectory by pediatric SBP tertile",
      x = "Years since diagnosis (binned)",
      y = "Median eGFR (IQR, capped)"
    ) +
    theme_minimal()

  predicted_df <- expand.grid(
    years_since_dx = seq(0, min(4, floor(max(traj6_model$years_since_dx, na.rm = TRUE))), by = 0.1),
    SBP_tertile = factor(c("Low SBP", "Mid SBP", "High SBP"), levels = c("Low SBP", "Mid SBP", "High SBP"))
  ) |>
    mutate(t_center = years_since_dx - mean(traj6_model$years_since_dx, na.rm = TRUE))

  predicted_df$predicted <- stats::predict(m1, newdata = predicted_df, re.form = NA)

  p6c_model <- ggplot(predicted_df, aes(x = years_since_dx, y = predicted, linetype = SBP_tertile)) +
    geom_line(linewidth = 1) +
    labs(
      title = "Figure 6C. Mixed-effects predicted eGFR trajectories by pediatric SBP tertile",
      x = "Years since T1D diagnosis",
      y = "Predicted eGFR (mL/min/1.73m2)",
      linetype = "Pediatric SBP tertile"
    ) +
    theme_minimal()

  save_plot(p6a_model, "Figure6A_spaghetti_byTertile.png", dirs, width = 8.2, height = 5.2)
  save_plot(p6b_model, "Figure6B_medianIQR_byTertile.png", dirs, width = 11, height = 4.2)
  save_plot(p6c_model, "Figure6C_predicted_byTertile.png", dirs, width = 7.5, height = 5)
} else {
  message("Skipping mixed-effects outputs because lme4/lmerTest/broom.mixed are not all installed.")
}

save_derived_rds(peds_patient_score, "peds_patient_score", dirs)
save_derived_rds(traj6, "traj6_with_tertiles", dirs)

message("Saved pediatric SBP analysis outputs to: ", dirs$results_dir)
