script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})

steps <- c(
  "01_build_analysis_data.R",
  "02_table1_and_bp_summary.R",
  "03_kidney_bp_analysis.R",
  "04_long_term_egfr_trajectories.R",
  "05_pediatric_sbp_models.R",
  "06_pediatric_percentile_models.R",
  "07_generate_manuscript_snippets.R",
  "08_publication_packet.R",
  "09_curate_figures.R"
)

for (step in steps) {
  message("\n==== Running ", step, " ====")
  source(file.path(script_dir, step), local = new.env(parent = globalenv()))
}

message("\nPipeline complete.")
