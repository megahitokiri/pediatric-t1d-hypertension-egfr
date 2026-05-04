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

  source_files <- vapply(
    sys.frames(),
    function(frame) if (!is.null(frame$ofile)) frame$ofile else NA_character_,
    character(1)
  )
  source_files <- source_files[!is.na(source_files)]
  if (length(source_files) > 0) {
    return(dirname(normalizePath(tail(source_files, 1), mustWork = FALSE)))
  }

  getwd()
}

stop_if_missing <- function(values) {
  missing <- names(values)[vapply(values, function(x) is.null(x) || !nzchar(x), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required environment variables: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

script_dir <- get_script_dir()

read_text_file <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

render_sql <- function(template, replacements) {
  out <- template
  for (key in names(replacements)) {
    out <- gsub(sprintf("\\{\\{%s\\}\\}", key), replacements[[key]], out, fixed = FALSE)
  }
  out
}

choose_backend <- function(connection_args) {
  requested <- tolower(connection_args$backend %||% "auto")
  has_odbc <- requireNamespace("DBI", quietly = TRUE) && requireNamespace("odbc", quietly = TRUE)
  has_rodbc <- requireNamespace("RODBC", quietly = TRUE)

  if (requested %in% c("rodbc", "odbc32", "windows")) {
    if (!has_rodbc) stop("HTA_DB_BACKEND=RODBC was requested, but package `RODBC` is not installed.", call. = FALSE)
    return("RODBC")
  }

  if (requested %in% c("odbc", "dbi")) {
    if (!has_odbc) stop("HTA_DB_BACKEND=odbc was requested, but packages `DBI` and `odbc` are not installed.", call. = FALSE)
    return("odbc")
  }

  if (!requested %in% c("auto", "")) {
    stop("Unsupported HTA_DB_BACKEND: ", connection_args$backend, ". Use auto, RODBC, or odbc.", call. = FALSE)
  }

  if (nzchar(connection_args$dsn %||% "") && has_rodbc) return("RODBC")
  if (.Platform$OS.type == "windows" && has_rodbc) return("RODBC")
  if (has_odbc) return("odbc")
  if (has_rodbc) return("RODBC")

  stop("No supported SQL Server connector found. Install `RODBC` or `DBI` + `odbc` in the VPN-enabled R environment.", call. = FALSE)
}

validate_connection_args <- function(connection_args) {
  if (nzchar(connection_args$dsn %||% "")) {
    stop_if_missing(connection_args[c("dsn")])
    return(invisible(TRUE))
  }

  stop_if_missing(connection_args[c("driver", "server", "database")])
  invisible(TRUE)
}

build_rodbc_connection_string <- function(connection_args) {
  parts <- character()

  if (nzchar(connection_args$dsn %||% "")) {
    parts <- c(parts, paste0("DSN=", connection_args$dsn))
  } else {
    parts <- c(
      parts,
      paste0("Driver={", connection_args$driver, "}"),
      paste0("Server=", connection_args$server)
    )
  }

  if (nzchar(connection_args$database %||% "")) {
    parts <- c(parts, paste0("Database=", connection_args$database))
  }

  if (nzchar(connection_args$uid %||% "")) {
    parts <- c(parts, paste0("UID=", connection_args$uid), paste0("PWD=", connection_args$pwd %||% ""))
  } else {
    parts <- c(parts, paste0("Trusted_Connection=", connection_args$trusted_connection))
  }

  if (nzchar(connection_args$encrypt %||% "")) {
    parts <- c(parts, paste0("Encrypt=", connection_args$encrypt))
  }

  if (nzchar(connection_args$trust_server_certificate %||% "")) {
    parts <- c(parts, paste0("TrustServerCertificate=", connection_args$trust_server_certificate))
  }

  paste0(paste(parts, collapse = ";"), ";")
}

run_query <- function(sql, connection_args) {
  backend <- choose_backend(connection_args)

  if (identical(backend, "odbc")) {
    odbc_args <- if (nzchar(connection_args$dsn %||% "")) {
      list(
        dsn = connection_args$dsn,
        Database = connection_args$database %||% NULL,
        UID = connection_args$uid %||% NULL,
        PWD = connection_args$pwd %||% NULL
      )
    } else {
      list(
        Driver = connection_args$driver,
        Server = connection_args$server,
        Database = connection_args$database,
        Trusted_Connection = connection_args$trusted_connection,
        Encrypt = connection_args$encrypt,
        TrustServerCertificate = connection_args$trust_server_certificate,
        UID = connection_args$uid %||% NULL,
        PWD = connection_args$pwd %||% NULL
      )
    }
    odbc_args <- odbc_args[!vapply(odbc_args, is.null, logical(1))]

    con <- do.call(DBI::dbConnect, c(list(drv = odbc::odbc()), odbc_args))
    on.exit(DBI::dbDisconnect(con), add = TRUE)
    return(DBI::dbGetQuery(con, sql))
  }

  conn_string <- build_rodbc_connection_string(connection_args)
  channel <- RODBC::odbcDriverConnect(conn_string)
  on.exit(RODBC::odbcClose(channel), add = TRUE)
  result <- RODBC::sqlQuery(channel, sql, stringsAsFactors = FALSE)
  if (!is.data.frame(result)) {
    stop("SQL query failed: ", paste(result, collapse = "\n"), call. = FALSE)
  }
  result
}

write_metadata <- function(path, values) {
  lines <- paste(names(values), unlist(values), sep = "=")
  writeLines(lines, con = path, useBytes = TRUE)
}

write_extract_csv <- function(data, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(data, path)
  } else {
    utils::write.csv(data, path, row.names = FALSE, na = "")
  }
}

connection_args <- list(
  backend = get_env_or_default("HTA_DB_BACKEND", if (.Platform$OS.type == "windows") "RODBC" else "auto"),
  dsn = get_env_or_default("HTA_DB_DSN"),
  driver = get_env_or_default("HTA_DB_DRIVER", "ODBC Driver 17 for SQL Server"),
  server = get_env_or_default("HTA_DB_SERVER", "pe-mssql"),
  database = get_env_or_default("HTA_DB_DATABASE", "IRB_00096551"),
  uid = get_env_or_default("HTA_DB_UID"),
  pwd = get_env_or_default("HTA_DB_PWD"),
  trusted_connection = get_env_or_default("HTA_DB_TRUSTED_CONNECTION", "Yes"),
  encrypt = get_env_or_default("HTA_DB_ENCRYPT", "Yes"),
  trust_server_certificate = get_env_or_default("HTA_DB_TRUST_SERVER_CERTIFICATE", "Yes")
)

validate_connection_args(connection_args)
message("Using SQL backend: ", choose_backend(connection_args))
if (nzchar(connection_args$dsn %||% "")) {
  message("Using ODBC DSN: ", connection_args$dsn)
}

output_dir <- get_env_or_default("HTA_OUTPUT_DIR", file.path(script_dir, "output"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

params <- list(
  MAX_FIRST_T1_AGE = get_env_or_default("HTA_MAX_FIRST_T1_AGE", "25"),
  MAX_LAST_VISIT_AGE = get_env_or_default("HTA_MAX_LAST_VISIT_AGE", "26"),
  MIN_FOLLOWUP_YEARS = get_env_or_default("HTA_MIN_FOLLOWUP_YEARS", "1"),
  MIN_VISITS = get_env_or_default("HTA_MIN_VISITS", "3")
)

sql_dir <- file.path(script_dir, "sql")
cohort_ctes_template <- read_text_file(file.path(sql_dir, "00_cohort_ctes.sql"))
cohort_template <- read_text_file(file.path(sql_dir, "01_cohort.sql"))
bp_template <- read_text_file(file.path(sql_dir, "02_bp.sql"))
labs_template <- read_text_file(file.path(sql_dir, "03_labs.sql"))
height_template <- read_text_file(file.path(sql_dir, "04_height_weight_bmi.sql"))

cohort_ctes_sql <- render_sql(cohort_ctes_template, params)
cohort_sql <- render_sql(cohort_template, c(params, COHORT_CTES = cohort_ctes_sql))
bp_sql <- render_sql(bp_template, c(params, COHORT_CTES = cohort_ctes_sql))
labs_sql <- render_sql(labs_template, c(params, COHORT_CTES = cohort_ctes_sql))
height_sql <- render_sql(height_template, c(params, COHORT_CTES = cohort_ctes_sql))

message("Running cohort extraction")
cohort <- run_query(cohort_sql, connection_args)
write_extract_csv(cohort, file.path(output_dir, "cohort.csv"))
saveRDS(cohort, file.path(output_dir, "cohort.rds"))
message("Wrote cohort extract: ", nrow(cohort), " rows")

message("Running BP extraction")
bp <- run_query(bp_sql, connection_args)
write_extract_csv(bp, file.path(output_dir, "bp.csv"))
saveRDS(bp, file.path(output_dir, "bp.rds"))
message("Wrote BP extract: ", nrow(bp), " rows")

message("Running lab extraction")
labs <- run_query(labs_sql, connection_args)
write_extract_csv(labs, file.path(output_dir, "labs.csv"))
saveRDS(labs, file.path(output_dir, "labs.rds"))
message("Wrote lab extract: ", nrow(labs), " rows")

message("Running height/weight/BMI extraction")
height_weight_bmi <- run_query(height_sql, connection_args)
write_extract_csv(height_weight_bmi, file.path(output_dir, "height_weight_bmi.csv"))
saveRDS(height_weight_bmi, file.path(output_dir, "height_weight_bmi.rds"))
message("Wrote height/weight/BMI extract: ", nrow(height_weight_bmi), " rows")

write_metadata(
  file.path(output_dir, "extract_metadata.txt"),
  c(
    list(
      extracted_at = format(Sys.time(), tz = Sys.timezone(), usetz = TRUE),
      backend = choose_backend(connection_args),
      dsn = connection_args$dsn %||% "",
      server = connection_args$server,
      database = connection_args$database
    ),
    params
  )
)

message("Extraction complete. Files written to: ", output_dir)
