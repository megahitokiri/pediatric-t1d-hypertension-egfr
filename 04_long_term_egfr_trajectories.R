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

creat_labs_u25 <- load_derived_rds("creat_labs_u25", dirs) |>
  mutate(
    RESULT_DATE = as.Date(RESULT_DATE),
    eGFR_u25 = as.numeric(eGFR_u25)
  ) |>
  filter(!is.na(RESULT_DATE), is.finite(eGFR_u25), eGFR_u25 > 0)

analysis_cohort <- load_derived_rds("analysis_cohort", dirs) |>
  mutate(first_t1_date = as.Date(first_t1_date))

egfr_cap <- 150

creat_labs_u25 <- creat_labs_u25 |>
  mutate(eGFR_u25_cap = pmin(eGFR_u25, egfr_cap))

long_term_ids <- creat_labs_u25 |>
  group_by(Patient_ID) |>
  summarise(
    n_creat = n(),
    first_date = min(RESULT_DATE),
    last_date = max(RESULT_DATE),
    followup_years = as.numeric(difftime(last_date, first_date, units = "days")) / 365.25,
    .groups = "drop"
  ) |>
  inner_join(
    analysis_cohort |> select(Patient_ID, first_t1_date, AgeAtFirstT1),
    by = "Patient_ID"
  ) |>
  filter(n_creat >= 5, followup_years >= 5, !is.na(AgeAtFirstT1), AgeAtFirstT1 < 15)

traj_long <- creat_labs_u25 |>
  inner_join(long_term_ids |> select(Patient_ID), by = "Patient_ID")

if (!"first_t1_date" %in% names(traj_long)) {
  traj_long <- traj_long |>
    left_join(
      analysis_cohort |> select(Patient_ID, first_t1_date),
      by = "Patient_ID"
    )
}

if (!"first_t1_date" %in% names(traj_long)) {
  traj_long <- traj_long |>
    mutate(first_t1_date = dplyr::coalesce(
      if ("first_t1_date.x" %in% names(traj_long)) as.Date(first_t1_date.x) else as.Date(NA),
      if ("first_t1_date.y" %in% names(traj_long)) as.Date(first_t1_date.y) else as.Date(NA)
    ))
}

traj_long <- traj_long |>
  mutate(
    first_t1_date = as.Date(first_t1_date),
    years_since_dx = as.numeric(difftime(RESULT_DATE, first_t1_date, units = "days")) / 365.25
  ) |>
  filter(!is.na(first_t1_date), is.finite(years_since_dx), years_since_dx >= 0) |>
  arrange(Patient_ID, RESULT_DATE)

table5_patientlevel <- traj_long |>
  group_by(Patient_ID) |>
  summarise(
    n_creat = n(),
    years_followed = max(years_since_dx, na.rm = TRUE),
    first_eGFR = first(eGFR_u25_cap),
    last_eGFR = last(eGFR_u25_cap),
    min_eGFR = min(eGFR_u25_cap, na.rm = TRUE),
    slope = as.numeric(stats::coef(stats::lm(eGFR_u25_cap ~ years_since_dx))[2]),
    .groups = "drop"
  )

table5_longterm <- table5_patientlevel |>
  summarise(
    N = n(),
    followup_median = median(years_followed, na.rm = TRUE),
    followup_Q1 = quantile(years_followed, 0.25, na.rm = TRUE),
    followup_Q3 = quantile(years_followed, 0.75, na.rm = TRUE),
    n_creat_median = median(n_creat, na.rm = TRUE),
    n_creat_Q1 = quantile(n_creat, 0.25, na.rm = TRUE),
    n_creat_Q3 = quantile(n_creat, 0.75, na.rm = TRUE),
    slope_median = median(slope, na.rm = TRUE),
    slope_Q1 = quantile(slope, 0.25, na.rm = TRUE),
    slope_Q3 = quantile(slope, 0.75, na.rm = TRUE),
    pct_slope_lt3 = mean(slope < -3, na.rm = TRUE) * 100,
    pct_slope_lt5 = mean(slope < -5, na.rm = TRUE) * 100
  )

traj_summary <- traj_long |>
  mutate(year_bin = floor(years_since_dx)) |>
  group_by(year_bin) |>
  summarise(
    eGFR_median = median(eGFR_u25_cap, na.rm = TRUE),
    eGFR_Q1 = quantile(eGFR_u25_cap, 0.25, na.rm = TRUE),
    eGFR_Q3 = quantile(eGFR_u25_cap, 0.75, na.rm = TRUE),
    n_points = n(),
    n_patients = n_distinct(Patient_ID),
    .groups = "drop"
  )

traj_categories <- table5_patientlevel |>
  mutate(
    Trajectory = case_when(
      slope >= -1 ~ "Stable (>= -1)",
      slope >= -3 ~ "Slow decline (-1 to -3)",
      TRUE ~ "Rapid decline (< -3)"
    )
  )

p5a <- ggplot(traj_long, aes(x = years_since_dx, y = eGFR_u25_cap, group = Patient_ID)) +
  geom_line(alpha = 0.12, linewidth = 0.35) +
  coord_cartesian(ylim = c(30, egfr_cap)) +
  labs(
    title = "Figure 5A. Individual eGFR trajectories from T1D diagnosis (long-term cohort)",
    subtitle = paste0("eGFR estimated using CKiD U25; values capped at ", egfr_cap, " mL/min/1.73m2"),
    x = "Years since T1D diagnosis",
    y = "eGFR (mL/min/1.73m2)"
  ) +
  theme_minimal()

p5b <- ggplot(traj_summary, aes(x = year_bin, y = eGFR_median)) +
  geom_ribbon(aes(ymin = eGFR_Q1, ymax = eGFR_Q3), alpha = 0.25) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 90, linetype = "dashed", alpha = 0.5) +
  coord_cartesian(ylim = c(30, egfr_cap)) +
  labs(
    title = "Figure 5B. Median eGFR trajectory after T1D diagnosis (IQR)",
    subtitle = paste0("Values capped at ", egfr_cap, " mL/min/1.73m2"),
    x = "Years since diagnosis (binned)",
    y = "Median eGFR (IQR)"
  ) +
  theme_minimal()

p5c <- ggplot(traj_categories, aes(x = Trajectory)) +
  geom_bar() +
  labs(
    title = "Figure 5C. Kidney function trajectory categories (>=5 years follow-up)",
    subtitle = paste0("Slopes estimated using capped eGFR (cap=", egfr_cap, ")"),
    x = "",
    y = "Number of patients"
  ) +
  theme_minimal()

write_table_csv(table5_longterm, "Table5_LongTerm.csv", dirs)
write_table_csv(table5_patientlevel, "Table5_LongTerm_patientlevel.csv", dirs)

save_derived_rds(traj_long, "traj_long", dirs)
save_derived_rds(table5_patientlevel, "table5_patientlevel", dirs)

save_plot(p5a, "Figure5A_spaghetti.png", dirs, width = 7.5, height = 5)
save_plot(p5b, "Figure5B_median_IQR.png", dirs, width = 7.5, height = 5)
save_plot(p5c, "Figure5C_trajectory_categories.png", dirs, width = 7.5, height = 4.5)

message("Saved long-term trajectory outputs to: ", dirs$results_dir)
