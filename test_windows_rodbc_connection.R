# Quick Windows/VPN connection test using the same RODBC pattern that works in SSMS/RStudio.
# Run this before the full extraction if you want to confirm connectivity.

if (!requireNamespace("RODBC", quietly = TRUE)) {
  stop("Package `RODBC` is required. Install it with install.packages('RODBC').", call. = FALSE)
}

driver <- Sys.getenv("HTA_DB_DRIVER", unset = "ODBC Driver 17 for SQL Server")
server <- Sys.getenv("HTA_DB_SERVER", unset = "pe-mssql")
database <- Sys.getenv("HTA_DB_DATABASE", unset = "IRB_00096551")
encrypt <- Sys.getenv("HTA_DB_ENCRYPT", unset = "Yes")
trust_server_certificate <- Sys.getenv("HTA_DB_TRUST_SERVER_CERTIFICATE", unset = "Yes")

conn_string <- paste0(
  "Driver={", driver, "};",
  "Server=", server, ";",
  "Database=", database, ";",
  "Trusted_Connection=Yes;",
  "Encrypt=", encrypt, ";",
  "TrustServerCertificate=", trust_server_certificate, ";"
)

channel <- RODBC::odbcDriverConnect(conn_string)
on.exit(RODBC::odbcClose(channel), add = TRUE)

patients_demo <- RODBC::sqlQuery(
  channel,
  "
  SELECT TOP 10 Patient_ID, SEX_CD, BIRTH_DT
  FROM dbo.Patients_IHC_UUHSC;
  ",
  stringsAsFactors = FALSE
)

if (!is.data.frame(patients_demo)) {
  stop("Connection test failed: ", paste(patients_demo, collapse = "\n"), call. = FALSE)
}

print(patients_demo)
message("RODBC connection test succeeded.")
