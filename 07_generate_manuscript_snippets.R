script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("readr", "dplyr"))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

dirs <- ensure_pipeline_dirs(script_dir)

read_optional_table <- function(filename) {
  path <- file.path(dirs$tables_dir, filename)
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

summary_tbl <- read_optional_table("Table7_PediatricHTN_summary.csv")
outcome_tbl <- read_optional_table("Table7_PediatricHTN_byKidneyOutcome.csv")
adult_model_tbl <- read_optional_table("Table7_PediatricHTN_AdultHTN_Model.csv")
slope_model_tbl <- read_optional_table("Table7_PediatricHTN_eGFRSlope_Model.csv")
rapid_model_tbl <- read_optional_table("Table7_PediatricHTN_RapidDecline_Model.csv")

if (is.null(summary_tbl) || is.null(outcome_tbl)) {
  stop("Table7 summary outputs are missing. Run 06_pediatric_percentile_models.R first.", call. = FALSE)
}

fmt_num <- function(x, digits = 1) {
  if (length(x) == 0 || is.na(x)) return("NA")
  format(round(as.numeric(x), digits), trim = TRUE, nsmall = digits)
}

extract_effect <- function(tbl, term_name) {
  if (is.null(tbl) || !("term" %in% names(tbl))) return(NULL)
  out <- tbl |> filter(term == term_name)
  if (nrow(out) == 0) return(NULL)
  out[1, , drop = FALSE]
}

n_total <- summary_tbl$N[1]
n_classified <- if ("N_with_height_informed_peds_bp_classification" %in% names(summary_tbl)) {
  summary_tbl$N_with_height_informed_peds_bp_classification[1]
} else {
  summary_tbl$N_with_peds_percentile_classification[1]
}
n_peds_htn <- summary_tbl$N_peds_htn_3plus_days[1]
pct_peds_htn <- summary_tbl$pct_peds_htn_3plus_days[1]
n_peds_stage2 <- summary_tbl$N_peds_stage2_3plus_days[1]
pct_peds_stage2 <- summary_tbl$pct_peds_stage2_3plus_days[1]
median_bp_days <- summary_tbl$median_peds_bp_days[1]
median_bmi <- summary_tbl$median_peds_bmi[1]

methods_lines <- c(
  "Pediatric blood pressure was reclassified using a guideline-aligned phenotype based on age, sex, and height.",
  "For visits before age 13 years, blood pressure was categorized using age-, sex-, and height-specific thresholds; for visits from 13 to <18 years, fixed adolescent thresholds were applied.",
  "Height and BMI were assigned from the nearest anthropometric measurement to each blood pressure day, allowing a maximum gap of 365 days.",
  "At the patient level, pediatric hypertension was defined as stage 1 or stage 2 hypertension on at least 3 pediatric blood pressure days.",
  "Kidney function was estimated with the CKiD U25 creatinine equation, using nearest valid height matched to each creatinine measurement within 365 days.",
  "Exploratory associations were then examined with adult hypertension phenotype, eGFR slope, and rapid kidney function decline."
)

results_lines <- c(
  paste0(
    "Among ", n_total, " patients in the analytic cohort, ",
    n_classified, " had at least one pediatric visit eligible for height-informed BP classification."
  ),
  paste0(
    n_peds_htn, " patients (", fmt_num(pct_peds_htn), "%) met the pediatric hypertension phenotype of stage 1/2 hypertension on at least 3 pediatric blood pressure days."
  ),
  paste0(
    n_peds_stage2, " patients (", fmt_num(pct_peds_stage2), "%) met the stage 2 pediatric hypertension phenotype on at least 3 days."
  ),
  paste0(
    "The median number of classified pediatric blood pressure days was ", fmt_num(median_bp_days), ", and the median pediatric BMI at classified visits was ", fmt_num(median_bmi), "."
  )
)

outcome_lines <- apply(outcome_tbl, 1, function(row) {
  paste0(
    "In the group '", row[["peds_htn_group"]], "', N=",
    row[["N"]], ", median eGFR slope=",
    fmt_num(row[["median_egfr_slope"]], 2),
    " mL/min/1.73m2/year and adult hypertension prevalence=",
    fmt_num(row[["pct_adult_htn"]], 1), "%."
  )
})

adult_effect <- extract_effect(adult_model_tbl, "peds_htn_3plus_days")
slope_effect <- extract_effect(slope_model_tbl, "peds_htn_burden")
rapid_effect <- extract_effect(rapid_model_tbl, "peds_htn_3plus_days")

model_lines <- c()
if (!is.null(adult_effect)) {
  model_lines <- c(
    model_lines,
    paste0(
      "In the exploratory adult hypertension model, the pediatric hypertension phenotype was associated with an odds ratio of ",
      fmt_num(adult_effect$estimate, 2),
      " (95% CI ",
      fmt_num(adult_effect$conf.low, 2), " to ",
      fmt_num(adult_effect$conf.high, 2),
      "; p=", fmt_num(adult_effect$p.value, 3), ")."
    )
  )
}
if (!is.null(slope_effect)) {
  model_lines <- c(
    model_lines,
    paste0(
      "In the linear kidney model, pediatric hypertension burden was associated with a beta coefficient of ",
      fmt_num(slope_effect$estimate, 3),
      " (95% CI ",
      fmt_num(slope_effect$conf.low, 3), " to ",
      fmt_num(slope_effect$conf.high, 3),
      "; p=", fmt_num(slope_effect$p.value, 3), ")."
    )
  )
}
if (!is.null(rapid_effect)) {
  model_lines <- c(
    model_lines,
    paste0(
      "In the rapid decline model, the pediatric hypertension phenotype showed an odds ratio of ",
      fmt_num(rapid_effect$estimate, 2),
      " (95% CI ",
      fmt_num(rapid_effect$conf.low, 2), " to ",
      fmt_num(rapid_effect$conf.high, 2),
      "; p=", fmt_num(rapid_effect$p.value, 3), ")."
    )
  )
}
if (length(model_lines) == 0) {
  model_lines <- c(
    "Adjusted exploratory models were generated only when sample size and outcome variation were sufficient; if model tables are absent, the current run should be interpreted as descriptive."
  )
}

interpretation_lines <- c(
  "This height-informed phenotype is methodologically stronger than a simple pediatric SBP tertile because it uses the pediatric clinical framework directly.",
  "For the manuscript, this analysis works best as a secondary or refinement analysis that complements the main BP-trajectory findings and strengthens the clinical relevance of the pediatric blood pressure signal.",
  "The most important sensitivity analyses are the height-to-BP matching window and the threshold used to define the patient-level pediatric hypertension phenotype."
)

manuscript_text <- c(
  "# Table 7 / Figure 7 Draft",
  "",
  "## Methods",
  methods_lines,
  "",
  "## Results",
  results_lines,
  "",
  outcome_lines,
  "",
  model_lines,
  "",
  "## Interpretation",
  interpretation_lines
)

abstract_text <- c(
  "Using a height-informed pediatric blood pressure phenotype that incorporated age, sex, and nearest height, we identified the subset of children who met pediatric hypertension criteria on at least 3 visits.",
  paste0(
    "Of ", n_classified, " children with classifiable pediatric blood pressure data, ",
    n_peds_htn, " (", fmt_num(pct_peds_htn), "%) met the pediatric hypertension phenotype and ",
    n_peds_stage2, " (", fmt_num(pct_peds_stage2), "%) met a stage 2 phenotype."
  ),
  "This analysis adds a clinically grounded pediatric hypertension definition that can be used as a refinement or sensitivity analysis alongside the trajectory-based blood pressure measures."
)

writeLines(manuscript_text, file.path(dirs$manuscript_dir, "Table7_Figure7_draft.md"))
writeLines(abstract_text, file.path(dirs$manuscript_dir, "Table7_abstract_snippet.txt"))

message("Saved manuscript snippets to: ", dirs$manuscript_dir)
