script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("readr", "dplyr", "ggplot2", "broom"))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

dirs <- ensure_pipeline_dirs(script_dir)

python_bin <- Sys.which("python3")
if (!nzchar(python_bin)) {
  stop("python3 is required for height-informed pediatric BP classification.", call. = FALSE)
}

bp_input <- file.path(dirs$derived_dir, "bp_with_nearest_anthro.csv")
bp_module <- file.path(dirname(script_dir), "BP_Percentiles", "BP_percentiles.py")
visit_output <- file.path(dirs$tables_dir, "Table7_PediatricHTN_visit_level.csv")
patient_output <- file.path(dirs$tables_dir, "Table7_PediatricHTN_patient_level.csv")
helper_script <- file.path(script_dir, "python", "classify_pediatric_htn_from_percentiles.py")

if (!file.exists(bp_input)) {
  stop("Missing derived input file: ", bp_input, call. = FALSE)
}
if (!file.exists(bp_module)) {
  stop("Missing BP percentile module: ", bp_module, call. = FALSE)
}

python_args <- shQuote(c(
  helper_script,
  "--input", bp_input,
  "--bp-module", bp_module,
  "--visit-output", visit_output,
  "--patient-output", patient_output,
  "--max-gap-days", "365",
  "--pediatric-age-cutoff", "13"
))

python_out <- system2(
  python_bin,
  python_args,
  stdout = TRUE,
  stderr = TRUE
)

python_status <- attr(python_out, "status")
if (!is.null(python_status) && python_status != 0) {
  stop(
    "Pediatric percentile classifier failed:\n",
    paste(python_out, collapse = "\n"),
    call. = FALSE
  )
}

peds_visit <- readr::read_csv(visit_output, show_col_types = FALSE)
peds_patient <- readr::read_csv(patient_output, show_col_types = FALSE)

table2_height_informed_peds_bp <- peds_visit |>
  mutate(
    age_at_bp = as.numeric(age_at_bp),
    valid_for_percentile_classification = as.integer(valid_for_percentile_classification),
    category_clean = dplyr::na_if(bp_category, ""),
    age_band_method = case_when(
      age_at_bp >= 1 & age_at_bp < 5 ~ "1-4 y",
      age_at_bp >= 5 & age_at_bp < 10 ~ "5-9 y",
      age_at_bp >= 10 & age_at_bp < 13 ~ "10-12 y",
      age_at_bp >= 13 & age_at_bp < 18 ~ "13-17 y",
      TRUE ~ NA_character_
    ),
    method_label = case_when(
      classification_method == "percentile_age_sex_height" ~ "Age/sex/height percentile",
      classification_method == "fixed_threshold_age_13plus" ~ "Adolescent fixed threshold",
      TRUE ~ NA_character_
    )
  ) |>
  filter(
    valid_for_percentile_classification == 1,
    !is.na(category_clean),
    !is.na(age_band_method)
  ) |>
  group_by(SEX_CD, age_band_method, method_label) |>
  summarise(
    n_patients = n_distinct(Patient_ID),
    n_bp_days = n(),
    normal_days = sum(category_clean == "Normal"),
    elevated_days = sum(category_clean == "Elevated"),
    stage1_days = sum(category_clean == "Stage1_HTN"),
    stage2_days = sum(category_clean == "Stage2_HTN"),
    htn_days = sum(category_clean %in% c("Stage1_HTN", "Stage2_HTN")),
    pct_elevated_or_htn = 100 * sum(category_clean %in% c("Elevated", "Stage1_HTN", "Stage2_HTN")) / n(),
    pct_htn = 100 * htn_days / n(),
    pct_stage2 = 100 * stage2_days / n(),
    .groups = "drop"
  ) |>
  arrange(SEX_CD, factor(age_band_method, levels = c("1-4 y", "5-9 y", "10-12 y", "13-17 y")))

analysis_cohort <- load_derived_rds("analysis_cohort", dirs)
kidney_bp_analysis <- load_derived_rds("kidney_bp_analysis", dirs)

model_df <- analysis_cohort |>
  left_join(
    kidney_bp_analysis |>
      select(Patient_ID, egfr_slope_per_year, HTN_2plus_days, n_BP_days_adult),
    by = "Patient_ID"
  ) |>
  left_join(peds_patient, by = "Patient_ID") |>
  mutate(
    median_peds_bmi = as.numeric(median_peds_bmi),
    median_peds_bmi_clean = ifelse(
      is.finite(median_peds_bmi) & median_peds_bmi >= 8 & median_peds_bmi <= 80,
      median_peds_bmi,
      NA_real_
    ),
    peds_htn_burden = as.numeric(peds_htn_burden),
    n_valid_peds_bp_days = as.numeric(n_valid_peds_bp_days),
    rapid_decline = ifelse(!is.na(egfr_slope_per_year), egfr_slope_per_year < -3, NA),
    sex_factor = factor(SEX_CD)
  )

table7_summary <- model_df |>
  summarise(
    N = n(),
    N_with_height_informed_peds_bp_classification = sum(!is.na(n_valid_peds_bp_days) & n_valid_peds_bp_days > 0),
    N_peds_htn_3plus_days = sum(peds_htn_3plus_days %in% 1, na.rm = TRUE),
    pct_peds_htn_3plus_days = mean(peds_htn_3plus_days %in% 1, na.rm = TRUE) * 100,
    N_peds_stage2_3plus_days = sum(peds_stage2_3plus_days %in% 1, na.rm = TRUE),
    pct_peds_stage2_3plus_days = mean(peds_stage2_3plus_days %in% 1, na.rm = TRUE) * 100,
    median_peds_bp_days = median(n_valid_peds_bp_days, na.rm = TRUE),
    N_with_plausible_peds_bmi = sum(!is.na(median_peds_bmi_clean)),
    median_peds_bmi = median(median_peds_bmi_clean, na.rm = TRUE)
  )

adult_htn_model_df <- model_df |>
  filter(
    !is.na(HTN_2plus_days),
    !is.na(peds_htn_3plus_days),
    !is.na(n_valid_peds_bp_days),
    n_valid_peds_bp_days >= 3,
    !is.na(AgeAtFirstT1),
    !is.na(YearsFromFirstT1ToLast)
  )

adult_htn_model_tbl <- tibble::tibble()
if (nrow(adult_htn_model_df) >= 25 && length(unique(adult_htn_model_df$HTN_2plus_days)) == 2) {
  adult_htn_model <- glm(
    HTN_2plus_days ~ peds_htn_3plus_days + median_peds_bmi_clean + AgeAtFirstT1 + sex_factor + YearsFromFirstT1ToLast + n_valid_peds_bp_days,
    data = adult_htn_model_df,
    family = binomial()
  )

  adult_htn_model_tbl <- broom::tidy(adult_htn_model, conf.int = TRUE, exponentiate = TRUE) |>
    mutate(model = "Adult HTN odds ratio from height-informed pediatric HTN phenotype")
}

slope_model_df <- model_df |>
  filter(
    !is.na(egfr_slope_per_year),
    !is.na(peds_htn_burden),
    !is.na(n_valid_peds_bp_days),
    n_valid_peds_bp_days >= 3,
    !is.na(AgeAtFirstT1),
    !is.na(YearsFromFirstT1ToLast)
  )

slope_model_tbl <- tibble::tibble()
rapid_model_tbl <- tibble::tibble()
if (nrow(slope_model_df) >= 25) {
  slope_model <- lm(
    egfr_slope_per_year ~ peds_htn_burden + median_peds_bmi_clean + AgeAtFirstT1 + sex_factor + YearsFromFirstT1ToLast + n_valid_peds_bp_days,
    data = slope_model_df
  )

  slope_model_tbl <- broom::tidy(slope_model, conf.int = TRUE) |>
    mutate(model = "Linear model for eGFR slope")

  rapid_model_df <- slope_model_df |>
    filter(!is.na(rapid_decline))

  if (nrow(rapid_model_df) >= 25 && length(unique(rapid_model_df$rapid_decline)) == 2) {
    rapid_model <- glm(
      rapid_decline ~ peds_htn_3plus_days + median_peds_bmi_clean + AgeAtFirstT1 + sex_factor + YearsFromFirstT1ToLast + n_valid_peds_bp_days,
      data = rapid_model_df,
      family = binomial()
    )

    rapid_model_tbl <- broom::tidy(rapid_model, conf.int = TRUE, exponentiate = TRUE) |>
      mutate(model = "Odds ratio for rapid eGFR decline (< -3)")
  }
}

table7_by_kidney <- model_df |>
  mutate(
    peds_htn_group = case_when(
      is.na(peds_htn_3plus_days) ~ "No classifiable pediatric BP phenotype",
      peds_htn_3plus_days == 1 ~ "Pediatric HTN on >=3 days",
      TRUE ~ "No pediatric HTN"
    )
  ) |>
  group_by(peds_htn_group) |>
  summarise(
    N = n(),
    N_with_slope = sum(!is.na(egfr_slope_per_year)),
    N_with_adult_htn_phenotype = sum(!is.na(HTN_2plus_days)),
    pct_adult_htn = mean(HTN_2plus_days %in% TRUE, na.rm = TRUE) * 100,
    median_egfr_slope = median(egfr_slope_per_year, na.rm = TRUE),
    median_followup_years = median(YearsFromFirstT1ToLast, na.rm = TRUE),
    .groups = "drop"
  )

p7a <- model_df |>
  filter(!is.na(peds_htn_3plus_days), !is.na(HTN_2plus_days)) |>
  mutate(
    peds_htn_group = ifelse(peds_htn_3plus_days == 1, "Pediatric HTN >=3 days", "No pediatric HTN"),
    adult_htn_label = ifelse(HTN_2plus_days, "Adult HTN", "No adult HTN")
  ) |>
  ggplot(aes(x = peds_htn_group, fill = adult_htn_label)) +
  geom_bar(position = "fill") +
  labs(
    title = "Adult HTN by height-informed pediatric BP phenotype",
    x = "",
    y = "Proportion",
    fill = "Adult BP phenotype"
  ) +
  theme_minimal()

p7b <- model_df |>
  filter(!is.na(peds_htn_3plus_days), !is.na(egfr_slope_per_year)) |>
  mutate(
    peds_htn_group = ifelse(peds_htn_3plus_days == 1, "Pediatric HTN >=3 days", "No pediatric HTN")
  ) |>
  ggplot(aes(x = peds_htn_group, y = egfr_slope_per_year)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  geom_hline(yintercept = -3, linetype = "dashed") +
  labs(
    title = "eGFR slope by height-informed pediatric BP phenotype",
    x = "",
    y = "eGFR slope (mL/min/1.73m2/year)"
  ) +
  theme_minimal()

write_table_csv(table7_summary, "Table7_PediatricHTN_summary.csv", dirs)
write_table_csv(table7_by_kidney, "Table7_PediatricHTN_byKidneyOutcome.csv", dirs)
write_table_csv(table2_height_informed_peds_bp, "Table2_HeightInformed_Pediatric_BP_Classification.csv", dirs)
if (nrow(adult_htn_model_tbl) > 0) {
  write_table_csv(adult_htn_model_tbl, "Table7_PediatricHTN_AdultHTN_Model.csv", dirs)
}
if (nrow(slope_model_tbl) > 0) {
  write_table_csv(slope_model_tbl, "Table7_PediatricHTN_eGFRSlope_Model.csv", dirs)
}
if (nrow(rapid_model_tbl) > 0) {
  write_table_csv(rapid_model_tbl, "Table7_PediatricHTN_RapidDecline_Model.csv", dirs)
}

save_plot(p7a, "Figure7A_AdultHTN_by_HeightInformedPedsHTN.png", dirs, width = 7.5, height = 5)
save_plot(p7b, "Figure7B_eGFRSlope_by_HeightInformedPedsHTN.png", dirs, width = 7.5, height = 5)

save_derived_rds(model_df, "pediatric_percentile_model_df", dirs)
save_derived_rds(model_df, "height_informed_pediatric_bp_model_df", dirs)

message("Saved height-informed pediatric BP models to: ", dirs$results_dir)
