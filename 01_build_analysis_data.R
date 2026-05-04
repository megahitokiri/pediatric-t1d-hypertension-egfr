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
egfr_max_height_gap_days <- as.numeric(Sys.getenv("HTA_EGFR_MAX_HEIGHT_GAP_DAYS", unset = "365"))

cohort <- load_extract_csv("cohort.csv", dirs) |>
  mutate(
    BIRTH_DT = as.Date(BIRTH_DT),
    first_t1_date = as.Date(first_t1_date),
    last_visit_date = as.Date(last_visit_date),
    AgeFirstT1_bin5 = cut(
      AgeAtFirstT1,
      breaks = c(0, 5, 10, 15, 20, 25),
      right = FALSE,
      labels = c("0-4", "5-9", "10-14", "15-19", "20-24")
    )
  )

bp <- load_extract_csv("bp.csv", dirs) |>
  mutate(
    OBSERVATION_DATE = as.Date(OBSERVATION_DATE),
    Systolic = as.numeric(Systolic),
    Diastolic = as.numeric(Diastolic)
  ) |>
  inner_join(
    cohort |> select(Patient_ID, SEX_CD, BIRTH_DT, first_t1_date),
    by = "Patient_ID"
  ) |>
  mutate(
    age_at_bp = as.numeric(difftime(OBSERVATION_DATE, BIRTH_DT, units = "days")) / 365.25,
    years_since_T1D = as.numeric(difftime(OBSERVATION_DATE, first_t1_date, units = "days")) / 365.25
  ) |>
  filter(
    !is.na(OBSERVATION_DATE),
    is.finite(Systolic),
    is.finite(Diastolic),
    is.finite(age_at_bp),
    is.finite(years_since_T1D),
    years_since_T1D >= 0
  ) |>
  mutate(
    age_band = case_when(
      age_at_bp >= 0 & age_at_bp < 3 ~ "0-2",
      age_at_bp >= 3 & age_at_bp < 6 ~ "3-5",
      age_at_bp >= 6 & age_at_bp < 9 ~ "6-8",
      age_at_bp >= 9 & age_at_bp < 12 ~ "9-11",
      age_at_bp >= 12 & age_at_bp < 15 ~ "12-14",
      age_at_bp >= 15 & age_at_bp < 18 ~ "15-17",
      TRUE ~ NA_character_
    )
  )

bp_unique_day <- bp |>
  mutate(bp_date = as.Date(OBSERVATION_DATE)) |>
  group_by(Patient_ID, bp_date) |>
  summarise(
    SBP = mean(Systolic, na.rm = TRUE),
    DBP = mean(Diastolic, na.rm = TRUE),
    age_at_bp = mean(age_at_bp, na.rm = TRUE),
    years_since_T1D = mean(years_since_T1D, na.rm = TRUE),
    SEX_CD = first(SEX_CD),
    .groups = "drop"
  )

bp_days_per_patient <- bp_unique_day |>
  count(Patient_ID, name = "BP_days")

analysis_cohort <- cohort |>
  inner_join(bp_days_per_patient, by = "Patient_ID") |>
  filter(BP_days >= 3)

bp_unique_day_analysis <- bp_unique_day |>
  inner_join(
    analysis_cohort |> select(Patient_ID, first_t1_date, AgeAtFirstT1, AgeAtLastVisit),
    by = "Patient_ID"
  )

adult_htn_patient <- bp_unique_day_analysis |>
  filter(age_at_bp >= 18) |>
  mutate(
    HTN_cat_adult = adult_htn_category(SBP, DBP),
    HTN_day = HTN_cat_adult %in% c("Stage1_HTN", "Stage2_HTN"),
    Stage2_day = HTN_cat_adult == "Stage2_HTN"
  ) |>
  group_by(Patient_ID) |>
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

height_weight_bmi <- load_extract_csv("height_weight_bmi.csv", dirs) |>
  mutate(
    OBSERVATION_DATE = as.Date(OBSERVATION_DATE),
    HEIGHT_CM = as.numeric(HEIGHT_CM),
    WEIGHT_KG = as.numeric(WEIGHT_KG),
    BMI = as.numeric(BMI)
  ) |>
  clean_anthropometrics() |>
  filter(!is.na(OBSERVATION_DATE)) |>
  inner_join(analysis_cohort |> select(Patient_ID), by = "Patient_ID")

bp_with_nearest_anthro <- match_nearest_anthro(bp_unique_day_analysis, height_weight_bmi)

labs_raw <- load_extract_csv("labs.csv", dirs)
labs_datetime_col <- if ("ResultDateTime" %in% names(labs_raw)) "ResultDateTime" else "RESULT_DTM"

labs <- labs_raw |>
  mutate(
    ResultDateTime = as.POSIXct(.data[[labs_datetime_col]], tz = "UTC"),
    RESULT_DATE = as.Date(ResultDateTime),
    Scr = suppressWarnings(as.numeric(OBS_VALUE))
  ) |>
  filter(
    ITEM_CODE %in% c("2160-0", "20025"),
    !is.na(RESULT_DATE),
    is.finite(Scr),
    Scr > 0
  ) |>
  inner_join(
    analysis_cohort |> select(Patient_ID, SEX_CD, BIRTH_DT, first_t1_date),
    by = "Patient_ID"
  ) |>
  mutate(
    age_at_lab = as.numeric(difftime(RESULT_DATE, BIRTH_DT, units = "days")) / 365.25,
    years_since_dx = as.numeric(difftime(RESULT_DATE, first_t1_date, units = "days")) / 365.25
  ) |>
  filter(age_at_lab >= 1, age_at_lab <= 25.99, years_since_dx >= 0)

labs_with_nearest_anthro <- match_nearest_anthro_to_index(
  labs,
  "RESULT_DATE",
  height_weight_bmi,
  prefix = "egfr_height"
)

egfr_height_qc <- tibble::tibble(
  Metric = c(
    "Creatinine labs eligible by age/time since T1D",
    "Creatinine labs with any valid matched height",
    paste0("Creatinine labs with valid height within ", egfr_max_height_gap_days, " days"),
    "Patients with eligible creatinine labs",
    "Patients with any valid matched height",
    paste0("Patients with valid height within ", egfr_max_height_gap_days, " days")
  ),
  Value = c(
    nrow(labs),
    sum(!is.na(labs_with_nearest_anthro$HEIGHT_CM)),
    sum(!is.na(labs_with_nearest_anthro$HEIGHT_CM) & labs_with_nearest_anthro$days_to_egfr_height <= egfr_max_height_gap_days),
    dplyr::n_distinct(labs$Patient_ID),
    dplyr::n_distinct(labs_with_nearest_anthro$Patient_ID[!is.na(labs_with_nearest_anthro$HEIGHT_CM)]),
    dplyr::n_distinct(labs_with_nearest_anthro$Patient_ID[
      !is.na(labs_with_nearest_anthro$HEIGHT_CM) &
        labs_with_nearest_anthro$days_to_egfr_height <= egfr_max_height_gap_days
    ])
  )
)

creat_labs_u25 <- labs_with_nearest_anthro |>
  filter(
    !is.na(HEIGHT_CM),
    is.finite(days_to_egfr_height),
    days_to_egfr_height <= egfr_max_height_gap_days
  ) |>
  mutate(eGFR_u25 = ckid_u25_egfr(Scr, age_at_lab, SEX_CD, HEIGHT_CM)) |>
  filter(is.finite(eGFR_u25), eGFR_u25 > 0)

egfr_patient_summary <- creat_labs_u25 |>
  arrange(Patient_ID, RESULT_DATE) |>
  group_by(Patient_ID) |>
  summarise(
    n_creat = n(),
    first_eGFR = first(eGFR_u25),
    last_eGFR = last(eGFR_u25),
    min_eGFR = min(eGFR_u25, na.rm = TRUE),
    any_eGFR_lt90 = any(eGFR_u25 < 90, na.rm = TRUE),
    any_eGFR_lt60 = any(eGFR_u25 < 60, na.rm = TRUE),
    .groups = "drop"
  )

egfr_slopes <- creat_labs_u25 |>
  arrange(Patient_ID, RESULT_DATE) |>
  group_by(Patient_ID) |>
  filter(n() >= 3) |>
  summarise(
    n_creat_used = n(),
    span_years = as.numeric(difftime(max(RESULT_DATE), min(RESULT_DATE), units = "days")) / 365.25,
    egfr_slope_per_year = as.numeric(stats::coef(stats::lm(eGFR_u25 ~ as.numeric(RESULT_DATE)))[2]) * 365.25,
    .groups = "drop"
  ) |>
  filter(span_years >= 1)

kidney_bp_analysis <- analysis_cohort |>
  left_join(egfr_patient_summary, by = "Patient_ID") |>
  left_join(egfr_slopes, by = "Patient_ID") |>
  left_join(adult_htn_patient, by = "Patient_ID")

save_derived_rds(cohort, "cohort_extract", dirs)
save_derived_rds(bp, "bp_enriched", dirs)
save_derived_rds(bp_unique_day, "bp_unique_day", dirs)
save_derived_rds(bp_unique_day_analysis, "bp_unique_day_analysis", dirs)
save_derived_rds(adult_htn_patient, "adult_htn_patient", dirs)
save_derived_rds(height_weight_bmi, "height_weight_bmi", dirs)
save_derived_rds(bp_with_nearest_anthro, "bp_with_nearest_anthro", dirs)
save_derived_rds(analysis_cohort, "analysis_cohort", dirs)
save_derived_rds(labs_with_nearest_anthro, "labs_with_nearest_anthro", dirs)
save_derived_rds(creat_labs_u25, "creat_labs_u25", dirs)
save_derived_rds(egfr_patient_summary, "egfr_patient_summary", dirs)
save_derived_rds(egfr_slopes, "egfr_slopes", dirs)
save_derived_rds(kidney_bp_analysis, "kidney_bp_analysis", dirs)

write_derived_csv(bp_with_nearest_anthro, "bp_with_nearest_anthro.csv", dirs)
write_derived_csv(labs_with_nearest_anthro, "creatinine_labs_with_nearest_height.csv", dirs)
write_derived_csv(creat_labs_u25, "EGFR_u25_labs_long.csv", dirs)
write_derived_csv(egfr_patient_summary, "EGFR_u25_patient_summary.csv", dirs)
write_table_csv(egfr_height_qc, "Table8_eGFR_HeightMatching_QC.csv", dirs)

message("Built analysis datasets in: ", dirs$derived_dir)
