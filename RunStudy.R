if (!(cdmName(cdm) %in% c("GOLD", "AURUM", "UiO", "SIDIAP", "CORIVA"))) {
  stop('databaseName should be one of the following: "GOLD", "AURUM", "UiO", "SIDIAP", "CORIVA"')
}

outputFolder <- here(results)
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder)
}

# foler with .RData files (this files will be used during the process and can
# be deleted at the end)
tempData <- "temp_data"
psFolder <- here(results, tempData)
if (!dir.exists(psFolder)) {
  dir.create(psFolder)
}

# create logger ----
log_file <- here(results, "log.txt")
if (file.exists(log_file)) {
  unlink(log_file)
}
logger <- create.logger()
logfile(logger) <- log_file
level(logger) <- "INFO"
info(logger, "CREATE LOGGER")

# STEP 0 Load study parameters and functions ----
info(logger, "STEP 0 INITIAL SETTINGS")

info(logger, "LOAD STUDY PARAMETERS")

minimumFrequencyPropensityScores <- 0.005
sampleSizePropensityScores <- 200000

cohortStem <- tolower(cohortStem)

covidCohortName <- paste0(cohortStem, "_covid")
vaccinatedCohortName <- paste0(cohortStem, "_vaccinated")
symptomsCohortName <- paste0(cohortStem, "_symptoms")
generalConditionsCohortName <- paste0(cohortStem, "_general_conditions")
databaseSpecificCohortName <- paste0(cohortStem, "_database_specific")
longcovidCohortName <- paste0(cohortStem, "_longcovid")
ncoCohortName <- paste0(cohortStem, "_negative_control_outcomes")
indexCohortName <- paste0(cohortStem, "_index")

allCohortNames <- c(
  covidCohortName, vaccinatedCohortName, symptomsCohortName,
  generalConditionsCohortName, databaseSpecificCohortName, 
  longcovidCohortName, ncoCohortName, indexCohortName
)

info(logger, "LOAD STUDY FUNCTIONS")
source(here("functions.R"))

snapshotDb <- dplyr::as_tibble(
  do.call(cbind.data.frame,CDMConnector::snapshot(cdm))
)
write_csv(snapshotDb, here(results, paste0("cdm_snapshot_", cdmName(cdm), ".csv")))

# STEP 1 Instantiate cohorts ----
if (isTRUE(instantiateCohorts)) {
  info(logger, "STEP 1 INSTANTIATE COHORS")

  ## 1.1 Instantiate json cohorts ----
  info(logger, "INSTANTIATE JSON COHORTS")
  source(here("1_InstantiateCohorts", "instantiate_json_cohorts.R"))

  ## 1.2 Instantiate long covid cohorts ----
  info(logger, "INSTANTIATE LONGCOVID COHORTS")
  source(here("1_InstantiateCohorts", "instantiate_longcovid_cohort.R"))

  ## 1.4 Instantiate nco cohorts ----
  info(logger, "INSTANTIATE NCO COHORTS")
  source(here("1_InstantiateCohorts", "instantiate_nco_cohort.R"))

  ## 1.5 Instantiate index cohorts ----
  info(logger, "INSTANTIATE INDEX COHORTS")
  source(here("1_InstantiateCohorts", paste0("instantiate_index_cohort_", cdmName(cdm), ".R")))
  
  ## 1.6 Export all cohort details in a single file
  info(logger, "EXPORT ALL COHORT DETAILS")
  allCohortDetails <- collectCohortSummary(cdm, allCohortNames) %>%
    mutate(cohort_table_name = gsub(paste0(cohortStem, "_"), "", cohort_table_name))
  write_csv(
    allCohortDetails,
    file = here(results, paste0("cohort_details_", cdmName(cdm), ".csv"))
  )
} else {
  info(logger, "STEP 1 READ COHORS")
  cdm <- cdmFromCon(
    con = db,
    cdmSchema = cdmDatabaseSchema,
    writeSchema = resultsDatabaseSchema,
    cohortTables = allCohortNames,
    cdmName = cdmName(cdm)
  )
}

# STEP 3 Propensity scores ----
if (runPropensityScores) {
  info(logger, "STEP 3 PROPENSITY SCORES")
  ## 3.1 Feature Extraction ----
  info(logger, "STEP 3 EXTRACTING FEATURES")
  source(here("3_PropensityScores", "extract_features.R"))
  ## 3.2 Define comparisons ----
  info(logger, "STEP 3 DEFINE COMPARISONS")
  source(here("3_PropensityScores", "define_comparisons.R"))
  ## 3.3 Compute propensity scores ----
  info(logger, "STEP 3 COMPUTE WEIGHTS")
  source(here("3_PropensityScores", "compute_propensity_scores.R"))
  ## 3.4 Standardized mean differences ----
  info(logger, "STEP 3 EVALUATE CONFOUNDING")
  source(here("3_PropensityScores", "evaluate_observed_confounding.R"))
  ## 3.5 Table one ----
  info(logger, "STEP 3 GET CHARACTERISTICS TABLE")
  source(here("3_PropensityScores", "get_characteristics_table.R"))
} else {
  comparisons <- read_csv(here(results, paste0("comparisons_", cdmName(cdm), ".csv")), show_col_types = FALSE)
}

# STEP 4 Outcome model ----
if (runOutcomeModel){
  source(here("4_OutcomeModel", "outcome_model.R"))
} else {
  ## skipping analysis
  info(logger, "STEP 4 SKIPPING ANALYSIS")
}

file_ext <- function(x) {
  pos <- regexpr("\\.([[:alnum:]]+)$", x)
  ifelse(pos > -1L, substring(x, pos + 1L), "")
}

filesToZip <- list.files(outputFolder, full.names = TRUE)
filesToZip <- filesToZip[file_ext(filesToZip) %in% c("csv", "txt")]
zip::zipr(
  zipfile = file.path(here(paste0("Results_", cdmName(cdm), ".zip"))),
  files = filesToZip
)

