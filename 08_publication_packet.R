script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})
source(file.path(script_dir, "R", "utils.R"))

require_packages(c("readr", "dplyr", "knitr", "rmarkdown"))
suppressPackageStartupMessages(library(dplyr))

dirs <- ensure_pipeline_dirs(script_dir)

read_table <- function(filename, subdir = dirs$tables_dir) {
  path <- file.path(subdir, filename)
  if (!file.exists(path)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE)
}

fmt_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", ifelse(x < 0.001, "<0.001", sprintf("%.3f", x)))
}

fmt_num <- function(x, digits = 1) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(is.na(x), "", format(round(x, digits), nsmall = digits, trim = TRUE))
}

fmt_med_iqr <- function(med, iqr, digits = 1) {
  paste0(fmt_num(med, digits), " (", fmt_num(iqr, digits), ")")
}

kable_md <- function(df, caption = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return("_Not available in this pipeline run._")
  }
  paste(capture.output(knitr::kable(df, format = "pipe", caption = caption)), collapse = "\n")
}

rename_if_present <- function(df, mapping) {
  if (is.null(df)) return(NULL)
  for (old in names(mapping)) {
    if (old %in% names(df)) names(df)[names(df) == old] <- mapping[[old]]
  }
  df
}

cohort_tbl <- read_table("Table1_T1D_T1DImputed_U26.csv")
peds_bp_tbl <- read_table("Table2_Pediatric_BP.csv")
peds_bp_class_tbl <- read_table("Table2_HeightInformed_Pediatric_BP_Classification.csv")
kidney_tbl <- read_table("Table3_adultBPonly.csv")
slope_age_tbl <- read_table("Table4A_slopes_by_age.csv")
slope_htn_tbl <- read_table("Table4B_adult_slopes_by_HTN.csv")
tests_tbl <- read_table("Table4C_tests.csv")
long_term_tbl <- read_table("Table5_LongTerm.csv")
peds_htn_summary <- read_table("Table7_PediatricHTN_summary.csv")
peds_htn_outcomes <- read_table("Table7_PediatricHTN_byKidneyOutcome.csv")
height_qc_tbl <- read_table("Table8_eGFR_HeightMatching_QC.csv")

adult_model_tbl <- read_table("Table7_PediatricHTN_AdultHTN_Model.csv")
slope_model_tbl <- read_table("Table7_PediatricHTN_eGFRSlope_Model.csv")
rapid_model_tbl <- read_table("Table7_PediatricHTN_RapidDecline_Model.csv")

addam_dir <- file.path(dirs$results_dir, "addam_genetics", "tables")
addam_qc <- read_table("Table_ADDAM_Cohort_QC.csv", addam_dir)
addam_ml <- read_table("Table_ADDAM_ML_LOOCV_Performance.csv", addam_dir)
addam_incremental <- read_table("Table_ADDAM_Incremental_LOOCV.csv", addam_dir)
addam_prs <- read_table("Table_ADDAM_PRS_Univariate.csv", addam_dir)

cohort_tbl_pub <- cohort_tbl |>
  rename_if_present(c(Characteristic = "Characteristic", Value = "Value"))

peds_bp_tbl_pub <- peds_bp_tbl |>
  transmute(
    Group = paste(SEX_CD, age_band_5yr, "y"),
    Patients = n_patients,
    `BP days` = n_bp_days,
    `SBP med (IQR)` = fmt_med_iqr(SBP_median, SBP_IQR),
    `DBP med (IQR)` = fmt_med_iqr(DBP_median, DBP_IQR)
  )

peds_bp_class_tbl_pub <- peds_bp_class_tbl |>
  transmute(
    Group = paste(SEX_CD, age_band_method),
    Method = method_label,
    Patients = n_patients,
    `BP days` = n_bp_days,
    `Elevated/HTN` = paste0(elevated_days + htn_days, " (", fmt_num(pct_elevated_or_htn, 1), "%)"),
    `HTN` = paste0(htn_days, " (", fmt_num(pct_htn, 1), "%)"),
    `Stage 2` = paste0(stage2_days, " (", fmt_num(pct_stage2, 1), "%)")
  )

kidney_tbl_pub <- kidney_tbl |>
  transmute(
    `Adult BP phenotype` = HTN_group2,
    Patients = N,
    `Creat/slope n` = paste0(N_with_creat, " / ", N_with_slope),
    `Last eGFR` = last_eGFR,
    `Slope/year` = slope,
    `eGFR <90/<60, %` = paste0(fmt_num(pct_any_eGFR_lt90), " / ", fmt_num(pct_any_eGFR_lt60))
  )

slope_age_tbl_pub <- slope_age_tbl |>
  transmute(
    `Age group` = age_group,
    Patients = N,
    `Slope/year` = slope,
    `Decline >3/>5, %` = paste0(fmt_num(pct_decline_gt3), " / ", fmt_num(pct_decline_gt5))
  )

slope_htn_tbl_pub <- slope_htn_tbl |>
  transmute(
    `Adult BP phenotype` = HTN_group_simple,
    Patients = N,
    `Slope/year` = slope,
    `Decline >3/>5, %` = paste0(fmt_num(pct_decline_gt3), " / ", fmt_num(pct_decline_gt5))
  )

tests_tbl_pub <- tests_tbl |>
  mutate(p_value = fmt_p(p_value)) |>
  rename_if_present(c(Test = "Comparison", p_value = "p value"))

long_term_tbl_pub <- if (!is.null(long_term_tbl) && nrow(long_term_tbl) > 0) {
  row <- long_term_tbl[1, ]
  tibble::tibble(
    Metric = c(
      "Patients with long-term kidney trajectory",
      "Follow-up, median [Q1, Q3], years",
      "Creatinine labs, median [Q1, Q3]",
      "eGFR slope, median [Q1, Q3]",
      "Slope < -3 per year, %",
      "Slope < -5 per year, %"
    ),
    Value = c(
      row$N,
      paste0(fmt_num(row$followup_median, 1), " [", fmt_num(row$followup_Q1, 1), ", ", fmt_num(row$followup_Q3, 1), "]"),
      paste0(fmt_num(row$n_creat_median, 1), " [", fmt_num(row$n_creat_Q1, 1), ", ", fmt_num(row$n_creat_Q3, 1), "]"),
      paste0(fmt_num(row$slope_median, 2), " [", fmt_num(row$slope_Q1, 2), ", ", fmt_num(row$slope_Q3, 2), "]"),
      fmt_num(row$pct_slope_lt3, 1),
      fmt_num(row$pct_slope_lt5, 1)
    )
  )
} else NULL

peds_htn_summary_pub <- if (!is.null(peds_htn_summary) && nrow(peds_htn_summary) > 0) {
  classified_col <- if ("N_with_height_informed_peds_bp_classification" %in% names(peds_htn_summary)) {
    "N_with_height_informed_peds_bp_classification"
  } else {
    "N_with_peds_percentile_classification"
  }
  tibble::tibble(
    Metric = c(
      "Analytic cohort",
      "With classifiable pediatric BP",
      "Pediatric HTN on >=3 days",
      "Pediatric stage 2 HTN on >=3 days",
      "Median classifiable pediatric BP days",
      "With plausible pediatric BMI",
      "Median pediatric BMI"
    ),
    Value = c(
      peds_htn_summary$N[1],
      peds_htn_summary[[classified_col]][1],
      paste0(peds_htn_summary$N_peds_htn_3plus_days[1], " (", fmt_num(peds_htn_summary$pct_peds_htn_3plus_days[1], 1), "%)"),
      paste0(peds_htn_summary$N_peds_stage2_3plus_days[1], " (", fmt_num(peds_htn_summary$pct_peds_stage2_3plus_days[1], 1), "%)"),
      peds_htn_summary$median_peds_bp_days[1],
      peds_htn_summary$N_with_plausible_peds_bmi[1],
      peds_htn_summary$median_peds_bmi[1]
    )
  )
} else NULL

peds_htn_outcomes_pub <- peds_htn_outcomes |>
  transmute(
    `Pediatric BP phenotype` = peds_htn_group,
    Patients = N,
    `eGFR slope` = paste0(N_with_slope, "; ", fmt_num(median_egfr_slope, 2), "/yr"),
    `Adult HTN` = paste0(N_with_adult_htn_phenotype, "; ", fmt_num(pct_adult_htn, 1), "%"),
    `Follow-up` = paste0(fmt_num(median_followup_years, 1), " y")
  )

height_qc_tbl_pub <- height_qc_tbl |>
  rename_if_present(c(Metric = "eGFR height-matching QC metric", Value = "Value"))

model_effects <- list(
  "Adult HTN" = adult_model_tbl,
  "eGFR slope" = slope_model_tbl,
  "Rapid eGFR decline" = rapid_model_tbl
) |>
  lapply(function(tbl) {
    if (is.null(tbl) || nrow(tbl) == 0) return(NULL)
    tbl |>
      filter(term %in% c("peds_htn_3plus_days", "peds_htn_burden", "median_peds_bmi_clean")) |>
      transmute(
        Outcome = dplyr::case_when(
          grepl("Adult HTN", first(model)) ~ "Adult HTN",
          grepl("eGFR slope", first(model)) ~ "eGFR slope",
          grepl("rapid", first(model), ignore.case = TRUE) ~ "Rapid decline",
          TRUE ~ first(model)
        ),
        Term = dplyr::recode(
          term,
          peds_htn_3plus_days = "Pediatric HTN >=3 days",
          peds_htn_burden = "Pediatric HTN burden",
          median_peds_bmi_clean = "Pediatric BMI"
        ),
        Estimate = fmt_num(estimate, 3),
        `95% CI` = paste0(fmt_num(conf.low, 3), " to ", fmt_num(conf.high, 3)),
        `p value` = fmt_p(p.value)
      )
  }) |>
  dplyr::bind_rows()

addam_qc_pub <- if (!is.null(addam_qc) && nrow(addam_qc) > 0) {
  tibble::tibble(
    Metric = c(
      "Participants",
      "Elevated BP cases",
      "Controls",
      "Feature columns",
      "Autoantibody features",
      "Ancestry features",
      "SNP dosage columns",
      "PRS available",
      "BMI/BMI-z available"
    ),
    Value = c(
      addam_qc$n_participants[1],
      paste0(addam_qc$elevated_bp_n[1], " (", fmt_num(addam_qc$elevated_bp_percent[1], 1), "%)"),
      addam_qc$control_n[1],
      addam_qc$feature_columns_n[1],
      addam_qc$autoantibody_features_n[1],
      addam_qc$ancestry_features_n[1],
      addam_qc$snp_dosage_columns_n[1],
      addam_qc$prs_available_n[1],
      paste0(addam_qc$bmi_available_n[1], " / ", addam_qc$bmiz_available_n[1])
    )
  )
} else NULL

addam_ml_pub <- addam_ml |>
  transmute(
    Model = dplyr::recode(
      model,
      `Logistic regression` = "Logit",
      `LASSO logistic regression` = "LASSO",
      `Random forest` = "Random forest",
      `Linear SVM` = "Linear SVM",
      .default = model
    ),
    AUC = fmt_num(auc, 3),
    Accuracy = fmt_num(accuracy, 3),
    `Sens/spec` = paste0(fmt_num(sensitivity, 3), " / ", fmt_num(specificity, 3)),
    Brier = fmt_num(brier_score, 3)
  )

addam_incremental_pub <- addam_incremental |>
  transmute(
    `Model set` = dplyr::case_when(
      model_set == "Clinical + autoantibodies" ~ "Clinical + Abs",
      model_set == "Clinical + autoantibodies + ancestry" ~ "Clinical + Abs + ancestry",
      model_set == "Clinical + autoantibodies + ancestry + PRS" ~ "Clinical + Abs + ancestry + PRS",
      model_set == "Clinical + autoantibodies + PRS" ~ "Clinical + Abs + PRS",
      model_set == "Clinical core" ~ "Clinical core",
      TRUE ~ model_set
    ),
    Predictors = predictors_n,
    AUC = fmt_num(auc, 3),
    `Delta AUC` = fmt_num(delta_auc_vs_clinical, 3),
    Brier = fmt_num(brier_score, 3)
  )

addam_prs_pub <- addam_prs |>
  filter(feature %in% c("GRS", "PRS", "IA2", "X96GAD", "ZnT8", "age_BP_per_year", "Cluster_A1c", "bmi", "bmiz")) |>
  transmute(
    Feature = dplyr::recode(
      feature,
      age_BP_per_year = "Age at BP",
      IA2 = "IA-2",
      X96GAD = "GAD",
      ZnT8 = "ZnT8",
      bmiz = "BMI z",
      bmi = "BMI",
      Cluster_A1c = "A1c",
      GRS = "GRS",
      .default = feature
    ),
    N = n,
    OR = fmt_num(odds_ratio, 2),
    `95% CI` = paste0(fmt_num(ci_low, 2), " to ", fmt_num(ci_high, 2)),
    `p value` = fmt_p(p_value)
  )

pub_date <- format(Sys.Date(), "%B %d, %Y")

publication_md <- c(
  "---",
  "title: \"Publication Table Packet\"",
  "subtitle: \"Elevated Blood Pressure and Kidney Function in Pediatric Type 1 Diabetes\"",
  "output:",
  "  word_document: default",
  "  html_document: default",
  "---",
  "",
  paste0("_Generated from the reproducible pipeline on ", pub_date, "._"),
  "",
  "# Table 1. Analytic Cohort Characteristics",
  kable_md(cohort_tbl_pub),
  "",
  "# Table 2. Height-Informed Pediatric Blood Pressure Classification",
  kable_md(peds_bp_class_tbl_pub),
  "",
  "# Supplemental Table S1. Raw Pediatric Blood Pressure Distribution",
  kable_md(peds_bp_tbl_pub),
  "",
  "# Table 3. Kidney Outcomes by Adult Hypertension Phenotype",
  kable_md(kidney_tbl_pub),
  "",
  "# Table 4A. eGFR Slopes by Age Group",
  kable_md(slope_age_tbl_pub),
  "",
  "# Table 4B. Adult eGFR Slopes by Adult Hypertension Phenotype",
  kable_md(slope_htn_tbl_pub),
  "",
  "# Table 4C. Statistical Comparisons",
  kable_md(tests_tbl_pub),
  "",
  "# Table 5. Longitudinal Kidney Function Trajectory Classes",
  kable_md(long_term_tbl_pub),
  "",
  "# Table 6. Height-Informed Pediatric BP Phenotype",
  kable_md(peds_htn_summary_pub),
  "",
  "# Table 7. Pediatric BP Phenotype by Kidney and Adult BP Outcomes",
  kable_md(peds_htn_outcomes_pub),
  "",
  "# Table 8. CKiD U25 Height Matching QC",
  kable_md(height_qc_tbl_pub),
  "",
  "# Table 9. Exploratory Adjusted Model Effects",
  kable_md(model_effects),
  "",
  "# Supplemental ADDAM Biobank / Genetics Tables",
  "These tables are included when the ADDAM genetics pipeline has been run locally.",
  "",
  "## ADDAM cohort QC",
  kable_md(addam_qc_pub),
  "",
  "## ADDAM ML LOOCV performance",
  kable_md(addam_ml_pub),
  "",
  "## ADDAM incremental model performance",
  kable_md(addam_incremental_pub),
  "",
  "## ADDAM PRS univariate model",
  kable_md(addam_prs_pub)
)

figure_plan_md <- c(
  "# Figure Strategy for the Paper",
  "",
  "## Recommended Main Figures",
  "",
  "1. Study design and cohorts: Utah population-based EHR cohort on the left, ADDAM biobank external/deep phenotype cohort on the right, with clear separation between discovery, phenotype validation, and genetics/PRS analyses.",
  "2. Pediatric BP phenotype over time: replace crude SBP tertile visuals with height-informed pediatric BP categories by years since T1D diagnosis; show classifiable visit density so sparse late follow-up is transparent.",
  "3. Kidney function trajectories: show CKiD U25 height-based eGFR trajectories and median/IQR by years since T1D diagnosis, with values capped only for display and uncapped values used in models unless explicitly stated.",
  "4. Clinical bridge figure: Sankey or stacked bars linking pediatric HTN phenotype to adult HTN phenotype and kidney trajectory category.",
  "5. ADDAM panel: model performance plus feature importance, clearly labeled as deep phenotype/genetics rather than the same pipeline as UDDB.",
  "",
  "## Figures to Move to Supplement",
  "",
  "1. SBP tertile spaghetti plots should become sensitivity/supplementary figures, because the height-informed pediatric BP phenotype is more clinically interpretable.",
  "2. Adult-only HTN density/box plots can support descriptive kidney findings but should not lead the manuscript unless the pediatric phenotype signal remains weak after height-based eGFR recalculation.",
  "3. Any mixed-effects plot with unsupported late follow-up years should be restricted to years with adequate patient counts and annotated with the supported horizon.",
  "",
  "## Visual Rules",
  "",
  "1. Every main figure should show denominator support either directly or in the caption.",
  "2. Avoid AI-generated infographic text inside panels; build final figures from real plots and add clean vector labels in PowerPoint/Illustrator if needed.",
  "3. Use one terminology family throughout: height-informed pediatric BP phenotype, adult HTN phenotype, CKiD U25 height-based eGFR, and ADDAM deep phenotype/genetics."
)

brief_lines <- c(
  "# Manuscript Data Packet for ChatGPT / Drafting",
  "",
  paste0("Generated: ", pub_date),
  "",
  "## Manuscript Direction",
  "",
  "The strongest framing is not simply that blood pressure is elevated in pediatric T1D. The stronger story is that a population-based pediatric T1D EHR cohort shows heterogeneous BP and kidney trajectories, and that replacing crude BP summaries with height-informed pediatric BP phenotyping plus height-based CKiD U25 eGFR creates a more clinically defensible risk phenotype. The ADDAM biobank/genetics analysis should be positioned as complementary deep phenotyping rather than as the same database pipeline.",
  "",
  "## Methods Language to Reuse",
  "",
  "Blood pressure days were collapsed to one mean SBP/DBP value per patient-day. For pediatric visits, BP was classified using a height-informed guideline phenotype: age-, sex-, and height-specific thresholds before age 13 years and fixed adolescent thresholds from age 13 to <18 years. Nearest valid anthropometrics were matched to BP days within a 365-day window.",
  "",
  "Kidney function was estimated using the CKiD U25 creatinine equation, eGFR = kappa x height(m) / serum creatinine, with kappa defined by age and sex. For the primary kidney trajectory analysis, each creatinine value required a nearest valid height measurement within 365 days.",
  "",
  "## Key Tables",
  "",
  "Use the publication table packet as the source of truth. Do not manually copy older CSV outputs from the Figures folder unless the table has been regenerated by this pipeline.",
  "",
  "## Results Anchors",
  "",
  kable_md(cohort_tbl_pub),
  "",
  kable_md(peds_htn_summary_pub),
  "",
  kable_md(peds_htn_outcomes_pub),
  "",
  kable_md(height_qc_tbl_pub),
  "",
  "## Suggested Abstract Interpretation",
  "",
  "Elevated BP was common in pediatric T1D, but the most manuscript-ready claim should emphasize clinically grounded phenotyping rather than raw tertiles. The current pipeline supports reporting how many children had classifiable pediatric BP data, how many met a sustained pediatric HTN phenotype, and whether that phenotype tracked with adult HTN and kidney slope outcomes after height-based CKiD U25 recalculation.",
  "",
  "## Guardrails for ChatGPT Drafting",
  "",
  "1. Do not invent p values, AUCs, or sample sizes not present in the tables.",
  "2. Distinguish UDDB EHR analyses from ADDAM genetics/PRS analyses.",
  "3. Avoid saying pediatric BP causes kidney decline; use associated with, tracked with, or identified a subgroup.",
  "4. State that SBP tertiles are exploratory/sensitivity if they are included.",
  "5. State how eGFR was calculated and how height was matched."
)

publication_rmd <- file.path(dirs$publication_dir, "Publication_Tables.Rmd")
brief_rmd <- file.path(dirs$publication_dir, "Manuscript_Data_Packet.Rmd")
figure_plan_path <- file.path(dirs$publication_dir, "Figure_Strategy.md")

writeLines(publication_md, publication_rmd, useBytes = TRUE)
writeLines(brief_lines, brief_rmd, useBytes = TRUE)
writeLines(figure_plan_md, figure_plan_path, useBytes = TRUE)

render_one <- function(input, output_format, output_file) {
  rmarkdown::render(
    input = input,
    output_format = output_format,
    output_file = output_file,
    output_dir = dirs$publication_dir,
    quiet = TRUE,
    envir = new.env(parent = globalenv())
  )
}

render_one(publication_rmd, "word_document", "HTA_T1D_Publication_Tables.docx")
render_one(publication_rmd, "html_document", "HTA_T1D_Publication_Tables.html")
render_one(brief_rmd, "word_document", "HTA_T1D_Manuscript_Data_Packet.docx")
render_one(brief_rmd, "html_document", "HTA_T1D_Manuscript_Data_Packet.html")

file.copy(brief_rmd, file.path(dirs$publication_dir, "HTA_T1D_Manuscript_Data_Packet.md"), overwrite = TRUE)

message("Saved publication packet to: ", dirs$publication_dir)
