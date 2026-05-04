if (!requireNamespace("readr", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE)) {
  stop("Packages `readr` and `dplyr` are required.", call. = FALSE)
}

library(dplyr)

get_env_or_default <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

`%||%` <- function(x, y) {
  if (!is.null(x) && length(x) > 0 && nzchar(x[1])) x else y
}

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- sub("^--file=", "", file_arg[1])
    script_path <- gsub("~\\+~", " ", script_path)
    return(dirname(normalizePath(script_path, mustWork = FALSE)))
  }
  getwd()
}

script_dir <- get_script_dir()
source(file.path(script_dir, "R", "utils.R"))

output_dir <- get_env_or_default("HTA_OUTPUT_DIR", file.path(script_dir, "output"))
egfr_max_height_gap_days <- as.numeric(Sys.getenv("HTA_EGFR_MAX_HEIGHT_GAP_DAYS", unset = "365"))

cohort <- readr::read_csv(file.path(output_dir, "cohort.csv"), show_col_types = FALSE) %>%
  mutate(
    BIRTH_DT = as.Date(BIRTH_DT),
    first_t1_date = as.Date(first_t1_date),
    last_visit_date = as.Date(last_visit_date)
  )

bp <- readr::read_csv(file.path(output_dir, "bp.csv"), show_col_types = FALSE) %>%
  mutate(
    OBSERVATION_DATE = as.Date(OBSERVATION_DATE),
    Systolic = as.numeric(Systolic),
    Diastolic = as.numeric(Diastolic)
  ) %>%
  inner_join(
    cohort %>% select(Patient_ID, SEX_CD, BIRTH_DT, first_t1_date),
    by = "Patient_ID"
  ) %>%
  mutate(
    age_at_bp = as.numeric(difftime(OBSERVATION_DATE, BIRTH_DT, units = "days")) / 365.25,
    years_since_T1D = as.numeric(difftime(OBSERVATION_DATE, first_t1_date, units = "days")) / 365.25
  ) %>%
  filter(!is.na(age_at_bp), !is.na(years_since_T1D), years_since_T1D >= 0)

bp_unique_day <- bp %>%
  group_by(Patient_ID, bp_date = as.Date(OBSERVATION_DATE)) %>%
  summarise(
    SBP = mean(Systolic, na.rm = TRUE),
    DBP = mean(Diastolic, na.rm = TRUE),
    age_at_bp = mean(age_at_bp, na.rm = TRUE),
    years_since_T1D = mean(years_since_T1D, na.rm = TRUE),
    .groups = "drop"
  )

adult_htn_patient <- bp_unique_day %>%
  filter(age_at_bp >= 18) %>%
  mutate(
    HTN_day = SBP >= 130 | DBP >= 80,
    Stage2_day = SBP >= 140 | DBP >= 90
  ) %>%
  group_by(Patient_ID) %>%
  summarise(
    n_BP_days_adult = n(),
    n_HTN_days_adult = sum(HTN_day, na.rm = TRUE),
    n_Stage2_days_adult = sum(Stage2_day, na.rm = TRUE),
    HTN_2plus_days = n_HTN_days_adult >= 2,
    Stage2_2plus_days = n_Stage2_days_adult >= 2,
    mean_SBP_adult = mean(SBP, na.rm = TRUE),
    mean_DBP_adult = mean(DBP, na.rm = TRUE),
    .groups = "drop"
  )

labs <- readr::read_csv(file.path(output_dir, "labs.csv"), show_col_types = FALSE) %>%
  mutate(
    ResultDateTime = as.POSIXct(ResultDateTime, tz = "UTC"),
    RESULT_DATE = as.Date(ResultDateTime),
    Scr = suppressWarnings(as.numeric(OBS_VALUE))
  ) %>%
  filter(!is.na(RESULT_DATE), is.finite(Scr), Scr > 0) %>%
  inner_join(
    cohort %>% select(Patient_ID, SEX_CD, BIRTH_DT, first_t1_date),
    by = "Patient_ID"
  ) %>%
  mutate(
    age_at_lab = as.numeric(difftime(RESULT_DATE, BIRTH_DT, units = "days")) / 365.25,
    years_since_dx = as.numeric(difftime(RESULT_DATE, first_t1_date, units = "days")) / 365.25
  ) %>%
  filter(age_at_lab >= 1, age_at_lab <= 25.99, years_since_dx >= 0)

height_weight_bmi <- readr::read_csv(file.path(output_dir, "height_weight_bmi.csv"), show_col_types = FALSE) %>%
  mutate(
    OBSERVATION_DATE = as.Date(OBSERVATION_DATE),
    HEIGHT_CM = as.numeric(HEIGHT_CM),
    WEIGHT_KG = as.numeric(WEIGHT_KG),
    BMI = as.numeric(BMI)
  ) %>%
  clean_anthropometrics()

labs <- match_nearest_anthro_to_index(
  labs,
  "RESULT_DATE",
  height_weight_bmi,
  prefix = "egfr_height"
) %>%
  filter(
    !is.na(HEIGHT_CM),
    is.finite(days_to_egfr_height),
    days_to_egfr_height <= egfr_max_height_gap_days
  ) %>%
  mutate(eGFR_u25 = ckid_u25_egfr(Scr, age_at_lab, SEX_CD, HEIGHT_CM)) %>%
  filter(is.finite(eGFR_u25), eGFR_u25 > 0)

egfr_patient_summary <- labs %>%
  arrange(Patient_ID, RESULT_DATE) %>%
  group_by(Patient_ID) %>%
  summarise(
    n_creat = n(),
    first_eGFR = first(eGFR_u25),
    last_eGFR = last(eGFR_u25),
    min_eGFR = min(eGFR_u25, na.rm = TRUE),
    any_eGFR_lt90 = any(eGFR_u25 < 90, na.rm = TRUE),
    any_eGFR_lt60 = any(eGFR_u25 < 60, na.rm = TRUE),
    .groups = "drop"
  )

egfr_slopes <- labs %>%
  arrange(Patient_ID, RESULT_DATE) %>%
  group_by(Patient_ID) %>%
  filter(n() >= 3) %>%
  summarise(
    n_creat_used = n(),
    span_years = as.numeric(difftime(max(RESULT_DATE), min(RESULT_DATE), units = "days")) / 365.25,
    egfr_slope_per_year = as.numeric(stats::coef(stats::lm(eGFR_u25 ~ as.numeric(RESULT_DATE)))[2]) * 365.25,
    .groups = "drop"
  ) %>%
  filter(span_years >= 1)

patient_level <- cohort %>%
  left_join(
    bp_unique_day %>%
      group_by(Patient_ID) %>%
      summarise(
        total_BP_days = n(),
        pediatric_BP_days = sum(age_at_bp < 18, na.rm = TRUE),
        adult_BP_days = sum(age_at_bp >= 18, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Patient_ID"
  ) %>%
  left_join(adult_htn_patient, by = "Patient_ID") %>%
  left_join(egfr_patient_summary, by = "Patient_ID") %>%
  left_join(egfr_slopes, by = "Patient_ID")

cohort_summary <- patient_level %>%
  summarise(
    N = n(),
    median_age_dx = median(AgeAtFirstT1, na.rm = TRUE),
    median_age_last = median(AgeAtLastVisit, na.rm = TRUE),
    median_followup_years = median(YearsFromFirstT1ToLast, na.rm = TRUE),
    median_total_BP_days = median(total_BP_days, na.rm = TRUE),
    median_pediatric_BP_days = median(pediatric_BP_days, na.rm = TRUE),
    median_adult_BP_days = median(adult_BP_days, na.rm = TRUE),
    pct_with_adult_htn = mean(HTN_2plus_days, na.rm = TRUE) * 100,
    pct_with_egfr_lt90 = mean(any_eGFR_lt90, na.rm = TRUE) * 100,
    pct_with_egfr_lt60 = mean(any_eGFR_lt60, na.rm = TRUE) * 100
  )

adult_htn_summary <- patient_level %>%
  mutate(
    HTN_group = case_when(
      is.na(HTN_2plus_days) ~ "No adult BP phenotype",
      HTN_2plus_days ~ "Adult HTN",
      TRUE ~ "No adult HTN"
    )
  ) %>%
  group_by(HTN_group) %>%
  summarise(
    N = n(),
    median_followup_years = median(YearsFromFirstT1ToLast, na.rm = TRUE),
    median_last_eGFR = median(last_eGFR, na.rm = TRUE),
    median_egfr_slope = median(egfr_slope_per_year, na.rm = TRUE),
    pct_any_eGFR_lt90 = mean(any_eGFR_lt90, na.rm = TRUE) * 100,
    pct_any_eGFR_lt60 = mean(any_eGFR_lt60, na.rm = TRUE) * 100,
    .groups = "drop"
  )

egfr_slope_summary <- patient_level %>%
  filter(!is.na(egfr_slope_per_year)) %>%
  mutate(
    slope_group = case_when(
      egfr_slope_per_year >= -1 ~ "Stable",
      egfr_slope_per_year >= -3 ~ "Slow decline",
      TRUE ~ "Rapid decline"
    )
  ) %>%
  group_by(slope_group) %>%
  summarise(
    N = n(),
    median_age_dx = median(AgeAtFirstT1, na.rm = TRUE),
    median_followup_years = median(YearsFromFirstT1ToLast, na.rm = TRUE),
    pct_adult_htn = mean(HTN_2plus_days, na.rm = TRUE) * 100,
    median_last_eGFR = median(last_eGFR, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(cohort_summary, file.path(output_dir, "pattern_scan_cohort_summary.csv"))
readr::write_csv(adult_htn_summary, file.path(output_dir, "pattern_scan_adult_htn_summary.csv"))
readr::write_csv(egfr_slope_summary, file.path(output_dir, "pattern_scan_egfr_slope_summary.csv"))
readr::write_csv(patient_level, file.path(output_dir, "pattern_scan_patient_level.csv"))

message("Pattern scan complete. Files written to: ", output_dir)
