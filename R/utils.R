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

pipeline_dirs <- function(script_dir = get_script_dir()) {
  output_dir <- Sys.getenv("HTA_OUTPUT_DIR", unset = "")
  if (!nzchar(output_dir)) {
    output_dir <- file.path(script_dir, "output")
  }

  list(
    script_dir = script_dir,
    output_dir = output_dir,
    derived_dir = file.path(script_dir, "derived"),
    results_dir = file.path(script_dir, "results"),
    tables_dir = file.path(script_dir, "results", "tables"),
    figures_dir = file.path(script_dir, "results", "figures"),
    manuscript_dir = file.path(script_dir, "results", "manuscript"),
    publication_dir = file.path(script_dir, "results", "publication")
  )
}

ensure_pipeline_dirs <- function(script_dir = get_script_dir()) {
  dirs <- pipeline_dirs(script_dir)
  dir.create(dirs$output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$derived_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$results_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$tables_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$figures_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$manuscript_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirs$publication_dir, recursive = TRUE, showWarnings = FALSE)
  dirs
}

require_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

load_extract_csv <- function(filename, dirs) {
  path <- file.path(dirs$output_dir, filename)
  if (!file.exists(path)) {
    stop("Missing extract file: ", path, call. = FALSE)
  }
  readr::read_csv(path, show_col_types = FALSE)
}

save_derived_rds <- function(object, name, dirs) {
  saveRDS(object, file.path(dirs$derived_dir, paste0(name, ".rds")))
}

load_derived_rds <- function(name, dirs) {
  path <- file.path(dirs$derived_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing derived file: ", path, call. = FALSE)
  }
  readRDS(path)
}

write_table_csv <- function(df, filename, dirs) {
  readr::write_csv(df, file.path(dirs$tables_dir, filename))
}

write_derived_csv <- function(df, filename, dirs) {
  readr::write_csv(df, file.path(dirs$derived_dir, filename))
}

save_plot <- function(plot, filename, dirs, width, height, dpi = 300) {
  ggplot2::ggsave(
    filename = file.path(dirs$figures_dir, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi
  )
}

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  sprintf(paste0("%.", digits, "f +/- %.", digits, "f"), mean(x), stats::sd(x))
}

fmt_med_q <- function(x, digits = 2) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_character_)
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  paste0(
    round(q[2], digits), " [",
    round(q[1], digits), ", ",
    round(q[3], digits), "]"
  )
}

adult_htn_category <- function(sbp, dbp) {
  dplyr::case_when(
    sbp < 120 & dbp < 80 ~ "Normal",
    sbp >= 120 & sbp < 130 & dbp < 80 ~ "Elevated",
    (sbp >= 130 & sbp < 140) | (dbp >= 80 & dbp < 90) ~ "Stage1_HTN",
    sbp >= 140 | dbp >= 90 ~ "Stage2_HTN",
    TRUE ~ NA_character_
  )
}

trajectory_category <- function(slope) {
  dplyr::case_when(
    slope >= -1 ~ "Stable",
    slope >= -3 ~ "Slow decline",
    TRUE ~ "Rapid decline"
  )
}

ckid_u25_k_creatinine <- function(age, sex) {
  age <- as.numeric(age)
  sex <- toupper(substr(as.character(sex), 1, 1))
  kappa <- rep(NA_real_, length(age))

  female <- sex == "F"
  male <- sex == "M"

  kappa[female & age >= 1 & age < 12] <- 36.1 * 1.008 ^ (age[female & age >= 1 & age < 12] - 12)
  kappa[male & age >= 1 & age < 12] <- 39.0 * 1.008 ^ (age[male & age >= 1 & age < 12] - 12)
  kappa[female & age >= 12 & age < 18] <- 36.1 * 1.023 ^ (age[female & age >= 12 & age < 18] - 12)
  kappa[male & age >= 12 & age < 18] <- 39.0 * 1.045 ^ (age[male & age >= 12 & age < 18] - 12)
  kappa[female & age >= 18 & age < 26] <- 41.4
  kappa[male & age >= 18 & age < 26] <- 50.8

  kappa
}

ckid_u25_egfr <- function(scr, age, sex, height_cm) {
  scr <- as.numeric(scr)
  age <- as.numeric(age)
  height_m <- as.numeric(height_cm) / 100
  kappa <- ckid_u25_k_creatinine(age, sex)

  egfr <- kappa * (height_m / scr)
  egfr[!is.finite(scr) | !is.finite(age) | !is.finite(height_m) | !is.finite(kappa)] <- NA_real_
  egfr[scr <= 0 | height_m <= 0] <- NA_real_
  egfr
}

clean_anthropometrics <- function(anthropometrics) {
  anthropometrics |>
    dplyr::mutate(
      HEIGHT_CM = dplyr::if_else(HEIGHT_CM >= 45 & HEIGHT_CM <= 230, HEIGHT_CM, NA_real_),
      WEIGHT_KG = dplyr::if_else(WEIGHT_KG >= 2 & WEIGHT_KG <= 250, WEIGHT_KG, NA_real_),
      BMI = dplyr::if_else(BMI >= 8 & BMI <= 80, BMI, NA_real_)
    )
}

match_nearest_anthro_to_index <- function(index_data, index_date_col, anthropometrics, prefix = "anthro") {
  index_data[[index_date_col]] <- as.Date(index_data[[index_date_col]])
  anthropometrics <- anthropometrics |>
    dplyr::mutate(OBSERVATION_DATE = as.Date(OBSERVATION_DATE))

  index_split <- split(index_data, index_data$Patient_ID)
  anthro_split <- split(anthropometrics, anthropometrics$Patient_ID)

  matched <- lapply(names(index_split), function(id) {
    patient_index <- index_split[[id]]
    patient_anthro <- anthro_split[[id]]

    date_col <- paste0(prefix, "_date")
    gap_col <- paste0("days_to_", prefix)

    if (is.null(patient_anthro) || nrow(patient_anthro) == 0) {
      patient_index[[date_col]] <- as.Date(NA)
      patient_index[[gap_col]] <- NA_real_
      patient_index$HEIGHT_CM <- NA_real_
      patient_index$WEIGHT_KG <- NA_real_
      patient_index$BMI <- NA_real_
      return(patient_index)
    }

    idx <- vapply(patient_index[[index_date_col]], function(one_date) {
      which.min(abs(as.numeric(difftime(patient_anthro$OBSERVATION_DATE, one_date, units = "days"))))
    }, integer(1))

    nearest <- patient_anthro[idx, c("OBSERVATION_DATE", "HEIGHT_CM", "WEIGHT_KG", "BMI")]
    names(nearest)[1] <- date_col

    patient_index[[date_col]] <- nearest[[date_col]]
    patient_index[[gap_col]] <- abs(as.numeric(difftime(nearest[[date_col]], patient_index[[index_date_col]], units = "days")))
    patient_index$HEIGHT_CM <- nearest$HEIGHT_CM
    patient_index$WEIGHT_KG <- nearest$WEIGHT_KG
    patient_index$BMI <- nearest$BMI
    patient_index
  })

  dplyr::bind_rows(matched)
}

match_nearest_anthro <- function(bp_unique_day, anthropometrics) {
  match_nearest_anthro_to_index(bp_unique_day, "bp_date", anthropometrics, prefix = "anthro")
}
