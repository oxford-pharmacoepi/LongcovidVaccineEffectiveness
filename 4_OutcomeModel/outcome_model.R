# create survival object
anyVaccineId <- cohortSet(cdm[[vaccinatedCohortName]]) %>%
  filter(cohort_name == "any_vaccine") %>%
  pull("cohort_definition_id")
anyCovidId <- cohortSet(cdm[[covidCohortName]]) %>%
  filter(cohort_name == "covid19") %>%
  pull("cohort_definition_id")
postAcuteCovidId <- cohortSet(cdm[[symptomsCohortName]]) %>%
  filter(cohort_name == "post_acute_covid19") %>%
  pull("cohort_definition_id")

unbalancedCovariates <- read_csv(here(results, paste0("asmd_", cdmName(cdm), ".csv")), show_col_types = FALSE) %>%
  filter(asmd_adjusted > 0.1) %>%
  filter(variable != "gp")
if (nrow(unbalancedCovariates) > 0) {
  unbalancedCovariates <- unbalancedCovariates %>%
    rowwise() %>%
    mutate(variable = if_else(
      grepl(")_window:", variable), 
      getCovariateId(variable),
      variable
    ))
}

longcovidCohortSet <- cohortSet(cdm[[longcovidCohortName]]) %>% collect()
ncoCohortSet <- cohortSet(cdm[[ncoCohortName]]) %>% collect()

events <- cdm[[indexCohortName]] %>%
  select("subject_id", "cohort_start_date") %>%
  distinct() %>%
  compute() %>%
  addEvent(cdm, "observation_period", name = "leave_db", eventDate = "observation_period_end_date") %>%
  addEvent(cdm, vaccinatedCohortName, anyVaccineId, window = c(1, NA), name = "next_vaccine") %>%
  addEvent(cdm, covidCohortName, anyCovidId, window = c(1, NA), name = "next_covid") %>%
  addEvent(cdm, "death", name = "death", eventDate = "death_date") %>%
  addMultipleEvent(cdm, longcovidCohortName, longcovidCohortSet$cohort_definition_id, c(1, NA), tolower(longcovidCohortSet$cohort_name)) %>%
  addMultipleEvent(cdm, ncoCohortName, ncoCohortSet$cohort_definition_id, c(1, NA), tolower(ncoCohortSet$cohort_name)) %>%
  addMultipleEvent(cdm, symptomsCohortName, postAcuteCovidId, c(1, NA), "next_post_acute_covid19")

# compute estimates
comparisonIds <-  comparisons$comparison_id[comparisons$skip == 0]
comparisonIds <- comparisonIds[lapply(comparisonIds, function(x) {
  unbalancedCovariates %>%
    filter(.data$comparison_id == .env$x) %>%
    pull("variable") %>%
    length() <= 10
}) %>% 
  unlist()]
outcomeNames <- c(
  ncoCohortSet$cohort_name,
  longcovidCohortSet$cohort_name, 
  "next_post_acute_covid19",
  "next_covid"
) %>%
  tolower()
censoringMethods = c("leave", "leave+vaccine")

for (comparison_id in comparisonIds) {
  
  exposure_cohort_id <- comparisons %>%
    filter(.data$comparison_id == .env$comparison_id) %>%
    pull("exposure_cohort_id")
  
  comparator_cohort_id <- comparisons %>%
    filter(.data$comparison_id == .env$comparison_id) %>%
    pull("comparator_cohort_id")
  
  unvaccinated <- grepl(
    "unvaccinated",
    comparisons$comparator_name[comparisons$comparison_id == comparison_id]
  )
  
  collectedCohort <-  cdm[[indexCohortName]] %>%
    filter(cohort_definition_id %in% !!c(exposure_cohort_id, comparator_cohort_id)) %>%
    mutate(group = if_else(
      cohort_definition_id == !!exposure_cohort_id,
      "exposure",
      "comparator"
    )) %>% 
    select("group", "subject_id", "cohort_start_date") %>%
    inner_join(events, by = c("subject_id", "cohort_start_date")) %>%
    collect() %>%
    mutate(comparison_id = .env$comparison_id) %>%
    mutate(leave_db = if_else(
      .env$unvaccinated & .data$group == "comparator", 
      pmin(.data$leave_db, .data$next_vaccine, na.rm = T), 
      .data$leave_db  
    ))
  
  unbalancedCovariatesComparison <- unbalancedCovariates %>%
    filter(.data$comparison_id == .env$comparison_id) %>%
    pull("variable")
  
  # load ps data
  load(paste0(psFolder, "/ps_model_", comparison_id, ".RData"))
  rm(countsFeatures, glmResult, selectedLassoFeatures)
  collectedCohort <- collectedCohort %>%
    inner_join(ps, by = c("subject_id", "group")) %>%
    mutate(weight = if_else(group == "exposure", 1-ps, ps))
  
  # add covariates that we are going to adjust
  nFeatureUnbalanced <- unbalancedCovariatesComparison[substr(unbalancedCovariatesComparison, 1, 1) != "f"]
  if (length(nFeatureUnbalanced) > 0) {
    load(here(results, tempData, "cohort.RData"))
    collectedCohort <- collectedCohort %>%
      left_join(
        cohort %>%
          select(all_of(c("subject_id", "cohort_start_date", nFeatureUnbalanced))),
        by = c("subject_id", "cohort_start_date")
      )
    rm(cohort)
  }
  
  # add features
  featureUnbalanced <- unbalancedCovariatesComparison[substr(unbalancedCovariatesComparison, 1, 1) == "f"]
  if (length(unbalancedCovariatesComparison) > 0) {
    load(here(results, tempData, "features.RData"))
    collectedCohort <- collectedCohort %>%
      left_join(
        features %>%
          filter(feature %in% .env$featureUnbalanced) %>%
          inner_join(
            collectedCohort %>%
              select("subject_id", "cohort_start_date"), 
            by = c("subject_id", "cohort_start_date")
          ) %>%
          mutate(value = 1) %>%
          pivot_wider(names_from = "feature", values_fill = 0),
        by = c("subject_id", "cohort_start_date")
      ) %>%
      mutate(across(all_of(featureUnbalanced), ~ if_else(is.na(.), 0, .)))
    rm(features)
  }
  collectedCohort <- collectedCohort %>% select(-"subject_id")
  save(
    collectedCohort, 
    file = here(results, tempData, paste0("outcome_", comparison_id, ".RData"))
  )
}

log4r::info(logger, "start")
for (comparison_id in comparisonIds) {
  tryCatch({
    load(here(results, tempData, paste0("outcome_", comparison_id, ".RData")))
    
    unbalancedCovariatesComparison <- unbalancedCovariates %>%
      filter(.data$comparison_id == .env$comparison_id) %>%
      pull("variable")
    
    # get the estimates and kaplan plots
    x <- getEstimates(
      collectedCohort, 
      outcomeNames, 
      censoringMethods,
      unbalancedCovariatesComparison
    )
    x$result <- x$result %>% mutate(comparison_id = comparison_id)
    x$survival_plot <- x$survival_plot %>% 
      filter(
        outcome_name %in% c(
          "longcovid_post_acute_covid19_28_365", 
          "longcovid_post_acute_covid19_90_365", 
          "longcovid_any_symptom_28_365", 
          "longcovid_any_symptom_90_365", 
          "next_post_acute_covid19",
          "next_covid"
        )
      ) %>%
      mutate(comparison_id = comparison_id)
    log4r::info(logger, paste0("finished comparison:", comparison_id))
    readr::write_csv(x$result, here::here(results, tempData, paste0("estimates_comparison_", comparison_id, ".csv")))
    readr::write_csv(x$survival_plot, here::here(results, tempData, paste0("survival_comparison_", comparison_id, ".csv"))
    )
  },
  error = function(e) {
    log4r::info(logger, paste0("failed comparison:", comparison_id))
  })
}

survivalData <- lapply(comparisonIds, function(x) {
  read_csv(here(results, tempData, paste0("survival_comparison_", x, ".csv")), show_col_types = FALSE) 
}) %>%
  bind_rows() %>%
  mutate(cdm_name = cdmName(cdm)) %>%
  write_csv(here(results, paste0("survival_plot_", cdmName(cdm), ".csv")))

dataEstimates <- lapply(comparisonIds, function(x) {
  read_csv(here(results, tempData, paste0("estimates_comparison_", x, ".csv")), show_col_types = FALSE)
})
dataEstimates <- bind_rows(dataEstimates)

outcomeNamesGroups <- tibble(
  outcome_name = ncoCohortSet$cohort_name,
  outcome_group = "nco"
) %>% 
  union_all(
    tibble(
      outcome_name = longcovidCohortSet$cohort_name
    ) %>%
      mutate(outcome_group = if_else(grepl("28_365", outcome_name), "longcovid28", "longcovid90"))
  ) %>%
  union_all(
    tibble(
      outcome_name = c("next_post_acute_covid19", "next_covid"), 
      outcome_group = "sensitivity"
    )
  ) %>%
  union_all(
    tibble(
      outcome_name = pacsCohortSet$cohort_name
    ) %>%
      mutate(outcome_group = case_when(
        grepl("leave_0-30", outcome_name) ~ "pacs 0-30",
        grepl("leave_31-90", outcome_name) ~ "pacs 31-90",
        grepl("leave_91-180", outcome_name) ~ "pacs 91-180",
        grepl("leave_181-365", outcome_name) ~ "pacs 181-365",
        grepl("covid+leave_0-30", outcome_name) ~ "pacs 0-30 covid",
        grepl("covid+leave_31-90", outcome_name) ~ "pacs 31-90 covid",
        grepl("covid+leave_91-180", outcome_name) ~ "pacs 91-180 covid",
        grepl("covid+leave_181-365", outcome_name) ~ "pacs 181-365 covid",
        grepl("vax+leave_0-30", outcome_name) ~ "pacs 0-30 vax",
        grepl("vax+leave_31-90", outcome_name) ~ "pacs 31-90 vax",
        grepl("vax+leave_91-180", outcome_name) ~ "pacs 91-180 vax",
        grepl("vax+leave_181-365", outcome_name) ~ "pacs 181-365 vax",
        grepl("covid+vax+leave_0-30", outcome_name) ~ "pacs 0-30 covid+vax",
        grepl("covid+vax+leave_31-90", outcome_name) ~ "pacs 31-90 covid+vax",
        grepl("covid+vax+leave_91-180", outcome_name) ~ "pacs 91-180 covid+vax",
        grepl("covid+vax+leave_181-365", outcome_name) ~ "pacs 181-365 covid+vax"
      ))
  ) %>%
  union_all(
    tibble(
      outcome_name = consequencesCohortSet$cohort_name,
      outcome_group = "consequences"
    )
  ) %>%
  mutate(outcome_name = tolower(outcome_name))

dataEstimates <- dataEstimates %>%
  left_join(outcomeNamesGroups, by = "outcome_name")

dataCalibration <- dataEstimates %>% 
  filter(.data$variable == "groupexposure") %>%
  filter(outcome_group == "nco")

write_csv(
  dataCalibration %>% mutate(cdm_name = cdmName(cdm)),
  here(results, paste0("negative_control_outcomes_", cdmName(cdm), ".csv"))
)

# add calibration
for (comparison_id in comparisonIds) {
  for (censoring_method in censoringMethods) {
    dataToFitModel <- dataCalibration %>%
      filter(.data$censoring_method == .env$censoring_method) %>%
      filter(.data$comparison_id == .env$comparison_id) %>%
      mutate(true_coef = 0) %>%
      filter(!is.na(hr)) 
    calibrationModel <- fitSystematicErrorModel(
      dataToFitModel$coef, 
      dataToFitModel$se_coef, 
      dataToFitModel$true_coef
    )
    dataOutcome <- dataEstimates %>%
      filter(.data$censoring_method == .env$censoring_method) %>%
      filter(.data$comparison_id == .env$comparison_id)
    dataOutcome <- dataOutcome %>%
      bind_cols(calibrateConfidenceInterval(
        dataOutcome$coef, 
        dataOutcome$se_coef, 
        calibrationModel
      )) %>%
      select(-c("coef", "hr", "se_coef", "z", "p", "lower_hr", "upper_hr")) %>%
      rename(
        "coef" = "logRr",
        "se_coef" = "seLogRr"
      ) %>%
      mutate(
        hr = exp(.data$coef),
        z = NA,
        p = NA,
        lower_hr = exp(logLb95Rr),
        upper_hr = exp(logUb95Rr),
        adjustment = "calibrated"
      ) %>%
      select(-c("logLb95Rr", "logUb95Rr"))
    
    dataEstimates <- dataEstimates %>%
      union_all(dataOutcome)
  }
}

write_csv(
  dataEstimates %>% mutate(cdm_name = cdmName(cdm)), 
  here(results, paste0("outcome_estimates_", cdmName(cdm), ".csv")) 
)
