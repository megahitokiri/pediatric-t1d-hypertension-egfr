# Windows/VPN setup for the secure SQL Server extract.
#
# Run this once in the same RStudio session before sourcing run_extract.R.
# This file does not store passwords. It assumes Windows/AD integrated auth.

get_setup_dir <- function() {
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

pipeline_dir <- get_setup_dir()

# Option A, recommended if you have a Windows ODBC DSN configured:
# In Windows, open "ODBC Data Sources (64-bit)" and create a System/User DSN
# pointing to server pe-mssql and database IRB_00096551. Then set the DSN name:
# Sys.setenv(HTA_DB_DSN = "IRB_00096551")

# Option B, no DSN needed: connect through RODBC using driver/server/database.
# These defaults match the connection string that works in the secure VM:
# Driver={ODBC Driver 17 for SQL Server}; Server=pe-mssql; Database=IRB_00096551.
# If this driver name is not installed, run:
# unique(RODBC::odbcDataSources())
# or check Windows "ODBC Data Sources (64-bit)" for the exact SQL Server driver.
Sys.setenv(
  HTA_DB_BACKEND = "RODBC",
  HTA_DB_DRIVER = "ODBC Driver 17 for SQL Server",
  HTA_DB_SERVER = "pe-mssql",
  HTA_DB_DATABASE = "IRB_00096551",
  HTA_DB_TRUSTED_CONNECTION = "Yes",
  HTA_DB_ENCRYPT = "Yes",
  HTA_DB_TRUST_SERVER_CERTIFICATE = "Yes",
  HTA_OUTPUT_DIR = file.path(pipeline_dir, "output")
)

# Optional study parameters. Leave these as-is unless you want to change the cohort.
Sys.setenv(
  HTA_MAX_FIRST_T1_AGE = "25",
  HTA_MAX_LAST_VISIT_AGE = "26",
  HTA_MIN_FOLLOWUP_YEARS = "1",
  HTA_MIN_VISITS = "3"
)

message("RODBC environment variables set. Now source run_extract.R in this same R session.")
