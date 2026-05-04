script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("dplyr", "readr"))
suppressPackageStartupMessages(library(dplyr))

dirs <- ensure_pipeline_dirs(script_dir)

analysis_cohort <- load_derived_rds("analysis_cohort", dirs)
bp_enriched <- load_derived_rds("bp_enriched", dirs)
bp_unique_day_analysis <- load_derived_rds("bp_unique_day_analysis", dirs)
adult_htn_patient <- load_derived_rds("adult_htn_patient", dirs)

table1_core <- tibble::tribble(
  ~Characteristic, ~Value,
  "N", as.character(nrow(analysis_cohort)),
  "Age at T1D diagnosis, mean +/- SD", fmt_mean_sd(analysis_cohort$AgeAtFirstT1),
  "Age at T1D diagnosis, median [Q1, Q3]", fmt_med_q(analysis_cohort$AgeAtFirstT1),
  "Age at last follow-up, mean +/- SD", fmt_mean_sd(analysis_cohort$AgeAtLastVisit),
  "Age at last follow-up, median [Q1, Q3]", fmt_med_q(analysis_cohort$AgeAtLastVisit),
  "Follow-up years, mean +/- SD", fmt_mean_sd(analysis_cohort$YearsFromFirstT1ToLast),
  "Follow-up years, median [Q1, Q3]", fmt_med_q(analysis_cohort$YearsFromFirstT1ToLast),
  "Visit count, mean +/- SD", fmt_mean_sd(analysis_cohort$VisitCount),
  "Visit count, median [Q1, Q3]", fmt_med_q(analysis_cohort$VisitCount),
  "Unique BP days, mean +/- SD", fmt_mean_sd(analysis_cohort$BP_days),
  "Unique BP days, median [Q1, Q3]", fmt_med_q(analysis_cohort$BP_days)
) |>
  bind_rows(
    analysis_cohort |>
      count(SEX_CD, name = "n") |>
      transmute(
        Characteristic = paste0("Sex: ", SEX_CD),
        Value = paste0(n, " (", round(100 * n / nrow(analysis_cohort), 1), "%)")
      )
  ) |>
  bind_rows(
    analysis_cohort |>
      count(FinalType, name = "n") |>
      transmute(
        Characteristic = paste0("Diabetes type: ", FinalType),
        Value = paste0(n, " (", round(100 * n / nrow(analysis_cohort), 1), "%)")
      )
  ) |>
  bind_rows(
    analysis_cohort |>
      count(AgeFirstT1_bin5, name = "n") |>
      transmute(
        Characteristic = paste0("Age at diagnosis bin: ", AgeFirstT1_bin5),
        Value = paste0(n, " (", round(100 * n / nrow(analysis_cohort), 1), "%)")
      )
  )

table2_pediatric_bp <- bp_unique_day_analysis |>
  filter(age_at_bp < 18) |>
  mutate(
    age_band_5yr = cut(
      age_at_bp,
      breaks = c(0, 5, 10, 15, 18),
      right = FALSE,
      labels = c("0-4", "5-9", "10-14", "15-17")
    )
  ) |>
  group_by(SEX_CD, age_band_5yr) |>
  summarise(
    n_bp_days = n(),
    n_patients = n_distinct(Patient_ID),
    SBP_median = median(SBP, na.rm = TRUE),
    SBP_IQR = IQR(SBP, na.rm = TRUE),
    DBP_median = median(DBP, na.rm = TRUE),
    DBP_IQR = IQR(DBP, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(SEX_CD, age_band_5yr)

adult_bp_summary <- if (nrow(adult_htn_patient) > 0) {
  tibble::tribble(
    ~Characteristic, ~Value,
    "N with adult BP", as.character(nrow(adult_htn_patient)),
    "Adult BP days, median [Q1, Q3]", fmt_med_q(adult_htn_patient$n_BP_days_adult),
    "Adult HTN >=2 days, %", as.character(round(100 * mean(adult_htn_patient$HTN_2plus_days, na.rm = TRUE), 1)),
    "Adult stage 2 HTN >=2 days, %", as.character(round(100 * mean(adult_htn_patient$Stage2_2plus_days, na.rm = TRUE), 1))
  )
} else {
  tibble::tribble(~Characteristic, ~Value, "Adult BP phenotype", "No adult BP records available")
}

table1_full <- bind_rows(table1_core, adult_bp_summary)

write_table_csv(table1_full, "Table1_T1D_T1DImputed_U26.csv", dirs)
write_table_csv(table2_pediatric_bp, "Table2_Pediatric_BP.csv", dirs)

save_derived_rds(table1_full, "table1_full", dirs)
save_derived_rds(table2_pediatric_bp, "table2_pediatric_bp", dirs)

message("Saved Table 1 and Table 2 outputs to: ", dirs$tables_dir)
