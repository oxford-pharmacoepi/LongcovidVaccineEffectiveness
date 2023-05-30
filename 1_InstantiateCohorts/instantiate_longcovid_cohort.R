attritionLC <- function(x, attrition, reason, ids) {
  newLine <- x %>%
    dplyr::group_by(.data$cohort_definition_id) %>%
    dplyr::mutate(obs = paste0(.data$subject_id, "_", .data$cohort_start_date)) %>%
    dplyr::summarise(
      number_records = dplyr::n_distinct(.data$obs),
      number_subjects = dplyr::n_distinct(.data$subject_id),
      .groups = "drop"
    ) %>%
    dplyr::collect() %>%
    dplyr::right_join(
      dplyr::tibble(cohort_definition_id = ids),
      by = "cohort_definition_id"
    ) %>%
    dplyr::mutate(
      number_records = dplyr::if_else(
        is.na(.data$number_records), 0, as.numeric(.data$number_records)
      ),
      number_subjects = dplyr::if_else(
        is.na(.data$number_subjects), 0, as.numeric(.data$number_subjects)
      )
    ) %>%
    dplyr::mutate(reason = .env$reason)
  if (is.null(attrition)) {
    attrition <- newLine %>%
      dplyr::mutate(
        reason_id = 1,
        excluded_records = NA,
        excluded_subjects = NA,
        step = as.character(NA)
      )
  } else {
    prevId <- max(attrition$reason_id)
    attrition <- attrition %>%
      dplyr::union_all(
        newLine %>%
          mutate(
            reason_id = .env$prevId + 1,
            excluded_records = 
              .env$attrition$number_records[prevId] - .data$number_records,
            excluded_subjects = 
              .env$attrition$number_subjects[prevId] - .data$number_subjects,
            step = as.character(NA)
          )
      )
  }
  return(attrition)
}

generateDenominatorPopulation <- function(cdm,
                                          cohortName,
                                          targetCohortId,
                                          influenzaId,
                                          anyCovidId) {
  target <- cdm[[cohortName]] %>%
    dplyr::filter(.data$cohort_definition_id == .env$targetCohortId) %>%
    dplyr::mutate(cohort_definition_id = 1) %>%
    dplyr::compute()
  
  influenza <- cdm[[cohortName]] %>%
    dplyr::filter(.data$cohort_definition_id == .env$influenzaId) %>%
    dplyr::select("subject_id", "influenza_date" = "cohort_start_date") %>%
    dplyr::compute()
  
  anyCovid <- cdm[[cohortName]] %>%
    dplyr::filter(.data$cohort_definition_id == .env$anyCovidId) %>%
    dplyr::select("subject_id", "any_covid_date" = "cohort_start_date") %>%
    dplyr::compute()
  
  attrition <- attritionLC(target, NULL, "Qualifying initial records", 1)
  
  # No prior target cohort in the previous 42 days
  target <- target %>%
    dplyr::group_by(.data$subject_id) %>%
    dbplyr::window_order(.data$cohort_start_date) %>%
    dplyr::mutate(target_id = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::compute()
  target <- target %>%
    dplyr::left_join(
      target %>%
        dplyr::mutate(target_id = .data$target_id + 1) %>%
        dplyr::select(
          "subject_id",
          "prev_target" = "cohort_start_date",
          "target_id"
        ),
      by = c("subject_id", "target_id")
    ) %>%
    dplyr::filter(
      is.na(.data$prev_target) |
        .data$cohort_start_date - .data$prev_target > 42
    ) %>%
    dplyr::select(-"target_id", -"prev_target") %>%
    dplyr::compute()
  
  attrition <- attritionLC(target, attrition, "No prior target in the previous 42 days", 1)
  
  # 180 days of prior observation
  target <- target %>%
    dplyr::inner_join(
      cdm[["observation_period"]] %>%
        dplyr::select(
          "subject_id" = "person_id", "observation_period_start_date"
        ),
      by = "subject_id"
    ) %>%
    dplyr::mutate(
      prior_history =
        .data$cohort_start_date - .data$observation_period_start_date
    ) %>%
    dplyr::filter(.data$prior_history >= 180) %>%
    dplyr::select(-"observation_period_start_date", -"prior_history") %>%
    dplyr::compute()
  
  attrition <- attritionLC(target, attrition, "180 days of prior history", 1)
  
  # find cohort_end_date
  target <- target %>%
    dplyr::left_join(
      cdm[["observation_period"]] %>%
        dplyr::select("subject_id" = "person_id", "observation_period_end_date"),
      by = "subject_id"
    ) %>%
    addEvent(cdm, cohortName, influenzaId, c(0, NA), "influenza_date") %>%
    addEvent(cdm, cohortName, anyCovidId, c(43, NA), "any_covid_date") %>%
    addEvent(cdm, "death", NULL, c(0, NA), "death_date", "death_date") %>%
    dplyr::mutate(cohort_end_date = as.Date(.data$cohort_start_date + lubridate::days(365))) %>%
    dplyr::mutate(cohort_end_date = pmin(
      .data$cohort_end_date,
      .data$observation_period_end_date,
      .data$influenza_date,
      .data$any_covid_date,
      .data$death_date,
      na.rm = TRUE
    )) %>%
    dplyr::select(
      "cohort_definition_id",
      "subject_id",
      "cohort_start_date",
      "cohort_end_date"
    ) %>%
    dplyr::compute()
  
  attr(target, "attrition") <- attrition %>%
    dplyr::mutate(step = "Generate denominator long-covid cohort")
  attr(target, "CohortSet") <- dplyr::tibble(
    cohort_definition_id = 1,
    cohort_name = "Denominator long-covid cohort"
  )
  
  return(target)
}

generateLongCovidCohorts <- function(cdm,
                                     targetCohortName,
                                     symptomCohortName,
                                     symptomCohortIds,
                                     symptomNames,
                                     windows,
                                     anySymptomIds) {
  target <- cdm[[targetCohortName]]
  attrition0 <- attr(target, "attrition")
  cohort_definition_id <- 1
  
  for (k in 1:length(symptomCohortIds)) {
    symptomCohortId <- symptomCohortIds[k]
    symptom <- cdm[[symptomCohortName]] %>%
      dplyr::filter(.data$cohort_definition_id == .env$symptomCohortId) %>%
      dplyr::select("subject_id", "symptom_date" = "cohort_start_date") %>%
      dplyr::compute()
    
    # add symptoms
    target_symptom <- target %>%
      dplyr::left_join(symptom, by = "subject_id") %>%
      dplyr::compute()
    
    # No symptom prior 180 days
    symptom180 <- target_symptom %>%
      dplyr::filter(.data$symptom_date < .data$cohort_start_date) %>%
      dplyr::filter(.data$symptom_date >= .data$cohort_start_date - lubridate::days(180)) %>%
      dplyr::select("subject_id", "cohort_start_date") %>%
      dplyr::compute()
    target_symptom <- target_symptom %>%
      dplyr::anti_join(symptom180, by = c("subject_id", "cohort_start_date")) %>%
      dplyr::compute()
    
    attrition.k <- attritionLC(
      target_symptom, 
      attrition0 %>% mutate(cohort_definition_id = .env$cohort_definition_id), 
      "No prior symptom in the previous 180 days",
      1
    )
    
    # symptom before the end date
    target_symptom <- target_symptom %>%
      dplyr::filter(.data$symptom_date <= .data$cohort_end_date) %>%
      dplyr::compute()
    
    for (i in 1:length(windows)) {
      
      # see if they have the symptom
      window_time <- windows[i]
      cohort_name <- paste0(targetCohortName, "_", symptomNames[k], "_", window_time, "_365")
      target_symptom_i <- target_symptom %>%
        dplyr::filter(.data$symptom_date >= .data$cohort_start_date + lubridate::days(.env$window_time)) %>%
        dplyr::mutate(cohort_definition_id = .env$cohort_definition_id) %>%
        dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") %>%
        dplyr::distinct() %>%
        dplyr::compute()
      
      attrition.k.i <- attritionLC(
        target_symptom_i, 
        attrition.k, 
        paste0("Symptom appears between ", window_time, " days and 365 days after the index date"),
        cohort_definition_id
      ) %>%
        mutate(cohort_definition_id = .env$cohort_definition_id)
      
      analysis <- dplyr::tibble(
        cohort_definition_id = .env$cohort_definition_id,
        cohort_name = .env$cohort_name,
        window = .env$window_time,
        target_cohort_name = .env$targetCohortName,
        symptom_name = .env$symptomNames[k]
      )
      
      if (cohort_definition_id == 1) {
        cohort <- target_symptom_i
        attrition <- attrition.k.i
        cohortSet <- analysis
      } else {
        cohort <- cohort %>%
          dplyr::union_all(target_symptom_i)
        cohortSet <- rbind(
          cohortSet,
          analysis
        )
        attrition <- rbind(
          attrition,
          attrition.k.i
        )
      }
      
      cohort_definition_id <- cohort_definition_id + 1
    }
  }
  
  cohort <- cohort %>% dplyr::compute()
  
  all_symptoms <- cdm[[symptomCohortName]] %>%
    dplyr::filter(.data$cohort_definition_id %in% .env$anySymptomIds) %>%
    dplyr::select(
      "cohort_definition_id",
      "subject_id",
      "symptom_date" = "cohort_start_date"
    ) %>%
    dplyr::inner_join(
      target %>%
        dplyr::select("subject_id", "cohort_start_date"),
      by = "subject_id"
    ) %>%
    dplyr::filter(.data$symptom_date < .data$cohort_start_date) %>%
    dplyr::filter(.data$symptom_date >= .data$cohort_start_date - lubridate::days(180)) %>%
    dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date") %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data$subject_id, .data$cohort_start_date) %>%
    dplyr::tally() %>%
    dplyr::filter(.data$n == !!length(anySymptomIds)) %>%
    dplyr::compute()
  
  target <- target %>%
    anti_join(all_symptoms, by = c("subject_id", "cohort_start_date")) %>%
    compute()
  
  attrition.all <- attritionLC(
    target, 
    attrition0, 
    "No prior symptom in the previous 180 days",
    1
  )
  
  anySymptomNames <- symptomNames[symptomCohortIds %in% anySymptomIds]
  
  for (i in 1:length(windows)) {
    window_time <- windows[i]
    cohort_name <- paste0(targetCohortName, "_any_symptom_", window_time, "_365")
    
    cohort.i <- cohort %>%
      dplyr::filter(
        .data$cohort_definition_id %in%  
          !!(cohortSet %>%
               dplyr::filter(.data$window == .env$window_time) %>%
               dplyr::filter(.data$symptom_name %in% .env$anySymptomNames) %>%
               dplyr::pull("cohort_definition_id"))
      ) %>%
      dplyr::mutate(cohort_definition_id = .env$cohort_definition_id) %>%
      dplyr::distinct() %>%
      dplyr::compute()
    
    # compile all the ones in the same window
    cohort <- cohort %>%
      dplyr::union_all(cohort.i) %>%
      dplyr::compute()
    
    attrition.all.i <- attritionLC(
      cohort.i, 
      attrition.all, 
      paste0("Symptom appears between ", window_time, " days and 365 days after the index date"),
      cohort_definition_id
    ) %>%
      dplyr::mutate(cohort_definition_id = .env$cohort_definition_id)
    
    analysis <- dplyr::tibble(
      cohort_definition_id = .env$cohort_definition_id,
      cohort_name = .env$cohort_name,
      window = .env$window_time,
      target_cohort_name = .env$targetCohortName,
      symptom_name = "any_symptom"
    )
    
    cohortSet <- rbind(cohortSet, analysis)
    attrition <- rbind(attrition, attrition.all.i) %>%
      mutate(step = if_else(is.na(step), "Symptom specific attrition", step))
    
    cohort_definition_id <- cohort_definition_id + 1
  }
  
  attr(cohort, "cohortSet") <- cohortSet
  attr(cohort, "attrition") <- attrition
  
  return(cohort)
}

anyCovidId <- cohortSet(cdm[[covidCohortName]]) %>% 
  filter(cohort_name == "covid19") %>% 
  pull("cohort_definition_id")
influenzaId <- cohortSet(cdm[[covidCohortName]]) %>% 
  filter(cohort_name == "influenza") %>% 
  pull("cohort_definition_id")
anySymptomIds <- cohortSet(cdm[[symptomsCohortName]]) %>%
  filter(cohort_name != "post_acute_covid19") %>%
  pull("cohort_definition_id")

cdm[["longcovid"]] <- generateDenominatorPopulation(
  cdm = cdm,
  cohortName = covidCohortName,
  targetCohortId = anyCovidId,
  influenzaId = influenzaId,
  anyCovidId = anyCovidId
)

symptomsCohortSet <- cohortSet(cdm[[symptomsCohortName]]) %>%
  collect()

longCovid <- generateLongCovidCohorts(
  cdm = cdm,
  targetCohortName = "longcovid",
  symptomCohortName = symptomsCohortName,
  symptomCohortIds = symptomsCohortSet$cohort_definition_id,
  symptomNames = symptomsCohortSet$cohort_name,
  windows = c(28, 90),
  anySymptomIds = anySymptomIds
)

longcovidAttrition <- attr(longCovid, "attrition")
longcovidCohortSet <- attr(longCovid, "cohortSet")

cdm[[longcovidCohortName]] <- newGeneratedCohortSet(
  cohortRef = computeQuery(longCovid, longcovidCohortName, FALSE, attr(cdm, "write_schema"), TRUE),
  cohortSetRef = insertTable(longcovidCohortSet, cdm, paste0(longcovidCohortName, "_set")),
  cohortAttritionRef = insertTable(longcovidAttrition, cdm, paste0(longcovidCohortName, "_attrition")),
  cohortCountRef = insertTable(getCohortCount(longCovid), cdm, paste0(longcovidCohortName, "_count"))
)
