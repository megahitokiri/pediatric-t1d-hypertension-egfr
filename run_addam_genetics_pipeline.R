script_dir <- local({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    script_path <- gsub("~\\+~", " ", sub("^--file=", "", file_arg[1]))
    dirname(normalizePath(script_path, mustWork = FALSE))
  } else getwd()
})

steps <- file.path(
  "addam_genetics",
  c(
    "01_addam_qc_and_prs_summary.R",
    "02_addam_model_summary.R",
    "02b_addam_incremental_models.R",
    "03_addam_manuscript_snippets.R"
  )
)

for (step in steps) {
  message("\n==== Running ", step, " ====")
  source(file.path(script_dir, step), local = new.env(parent = globalenv()))
}

message("\nADDAM genetics/PRS pipeline complete.")
