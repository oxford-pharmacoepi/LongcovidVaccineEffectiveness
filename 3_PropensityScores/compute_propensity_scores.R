covid_test_id <- cohortSet(cdm[[covidCohortName]]) %>%
  filter(cohort_name == "any_covid19_test") %>%
  pull("cohort_definition_id")

cohort <- cdm[[indexCohortName]] %>%
  select("subject_id", "cohort_start_date") %>%
  distinct() %>%
  compute() %>%
  addAge(cdm) %>%
  addAgeGroup(ageGroup = list(
    c(15,19), c(20,24), c(25,29), c(30,34), c(35,39), c(40,44), c(45,49),
    c(50,54), c(55,59), c(60,64), c(65,69), c(70,74), c(75,79), c(80,84),
    c(85,89), c(90,94), c(95,99), c(100,150)
  )) %>%
  addPriorHistory(cdm) %>%
  addSex(cdm)

cohort <- cohort %>%
  addNumberVisit(cdm, c(NA, -366), "number_visit_3")

cdm[["visit_occurrence"]] <- cdm[["visit_occurrence"]] %>%
  dplyr::filter(lubridate::year(.data$visit_start_date) >= 2020) %>%
  dplyr::compute()

cohort <- cohort %>%
  addNumberVisit(cdm, c(-365,  -31), "number_visit_2") %>%
  addNumberVisit(cdm, c( -30,   -1), "number_visit_1") %>%
  addNumberEvent(cdm, covidCohortName, list(cohort_definition_id = covid_test_id), c(  NA, -366), "number_covid_test_3") %>%
  addNumberEvent(cdm, covidCohortName, list(cohort_definition_id = covid_test_id), c(-365,  -31), "number_covid_test_2") %>%
  addNumberEvent(cdm, covidCohortName, list(cohort_definition_id = covid_test_id), c( -30,   -1), "number_covid_test_1") %>%
  addRegionalData(cdm) %>%
  addGPIdentifier(cdm) %>%
  collect()

save(cohort, file = here(results, tempData, "cohort.RData"))

rm(cohort)

comparisons <- comparisons %>%
  mutate(
    skip = NA,
    number_exposures = NA,
    number_comparators = NA
  )

psCoefficients <- NULL

for (k in comparisons$comparison_id) {
  
  load(here(results, tempData, "cohort.RData"))
  
  cohort.k <- cdm[[indexCohortName]] %>%
    filter(cohort_definition_id %in% !!c(comparisons$exposure_cohort_id[comparisons$comparison_id == k], 
                                         comparisons$comparator_cohort_id[comparisons$comparison_id == k])) %>%
    mutate(group = if_else(
      cohort_definition_id == !!comparisons$exposure_cohort_id[comparisons$comparison_id == k],
      "exposure",
      "comparator"
    )) %>% 
    select("group", "subject_id", "cohort_start_date") %>%
    collect() %>%
    inner_join(cohort, by = c("subject_id", "cohort_start_date")) %>%
    mutate(gp = as.character(gp)) %>%
    mutate(group = factor(group, c("comparator", "exposure")))
  
  rm(cohort)
  
  comparisons$number_exposures[comparisons$comparison_id == k] <- sum(cohort.k$group == "exposure")
  comparisons$number_comparators[comparisons$comparison_id == k] <- sum(cohort.k$group == "comparator")
  
  if (sum(cohort.k$group == "exposure") < 1 || 
      sum(cohort.k$group == "comparator") < 1) {
    
    comparisons$skip[comparisons$comparison_id == k] <- 1
    cat(paste0("Skipping comparison ", k, ". A group is empty.\n"))
    
  } else {
    
    comparisons$skip[comparisons$comparison_id == k] <- 0
    cat(paste0("Starting comparison ", k, ".\n"))
    
    load(here(results, tempData, "features.RData"))
    
    features <- features %>%
      inner_join(
        cohort.k %>%
          select("subject_id", "group", "cohort_start_date"),
        by = c("subject_id", "cohort_start_date")
      ) %>%
      select(-"cohort_start_date")
    
    countsFeatures <- features %>%
      group_by(feature) %>%
      summarise(n = n()) %>%
      mutate(freq = n / nrow(cohort.k)) %>%
      filter(freq >= minimumFrequencyPropensityScores) %>%
      select(feature)
    
    if (nrow(cohort.k) > sampleSizePropensityScores) {
      subcohort.k <- cohort.k[sample(1:nrow(cohort.k), sampleSizePropensityScores),]
    } else {
      subcohort.k <- cohort.k
    }
    
    subfeatures <- features %>%
      inner_join(countsFeatures, by = "feature") %>%
      inner_join(
        subcohort.k %>% select("subject_id", "group"),
        by = c("subject_id", "group")
      ) %>%
      mutate(value = 1) %>%
      pivot_wider(names_from = "feature", values_from = "value", values_fill = 0) %>%
      right_join(
        subcohort.k %>% select("subject_id", "group"),
        by = c("subject_id", "group")
      ) %>%
      mutate(across(starts_with("f"), ~ if_else(is.na(.x), 0, .x)))
    
    if (nrow(subfeatures) != nrow(subcohort.k)) {stop("different numer of rows in subfeatures and subcohort.k")}
    
    x <- data.matrix(subfeatures %>% select(-c("subject_id", "group")))
    y <- subfeatures$group
    lambdas <- 10^seq(2, -3, by = -.1)
    lasso_reg <- cv.glmnet(x, y, lambda = lambdas, standardize = TRUE, nfolds = 5, family = "binomial", alpha = 1)
    coef.lasso_reg <- coef(lasso_reg, s = lasso_reg$lambda.1se)
    selectedLassoFeatures <- names(coef.lasso_reg[(coef.lasso_reg[,1]!=0),1])
    selectedLassoFeatures <- selectedLassoFeatures[selectedLassoFeatures != "(Intercept)"]
    
    rm(x, y, lasso_reg)
    
    subfeatures <- subfeatures %>%
      select(all_of(c("subject_id", "group", selectedLassoFeatures))) %>%
      right_join(
        subcohort.k %>% select(-c("cohort_start_date", "region")),
        by = c("subject_id", "group")
      ) %>%
      mutate(age2 = age ** 2)
    
    if (nrow(subfeatures) != nrow(subcohort.k)) {stop("different numer of rows in subfeatures and subcohort.k")}
    
    if (cdmName(cdm) %in% c("CORIVA", "UiO")) {
      subfeatures <- subfeatures %>%
        select(-c("gp", "new_cases_7d_100k", "new_tests_7d_100k", "dose1_fraction", "dose2_fraction"))
    }
    
    glmResult <- glm(group ~ . - subject_id, data = subfeatures, family = binomial(link = "logit"))
    
    if (is.null(selectedLassoFeatures)) {
      features <- cohort.k %>%
        select(-c("cohort_start_date", "region")) %>%
        mutate(age2 = age ** 2)
    } else {
      features <- features %>%
        inner_join(tibble(feature = selectedLassoFeatures), by = "feature") %>%
        mutate(value = 1) %>%
        pivot_wider(names_from = "feature", values_from = "value", values_fill = 0) %>%
        right_join(
          cohort.k %>% select(-c("cohort_start_date", "region")),
          by = c("subject_id", "group")
        ) %>%
        mutate(across(starts_with("f"), ~ if_else(is.na(.x), 0, .x))) %>%
        mutate(age2 = age ** 2)
    }
    
    if (nrow(features) != nrow(cohort.k)) {stop("different numer of rows in subfeatures and subcohort.k")}
    
    features <- features %>%
      mutate(age_group = if_else(age_group %in% subfeatures$age_group, age_group, as.character(NA)))
    
    if (cdmName(cdm) %in% c("CORIVA", "UiO")) {
      features <- features %>%
        select(-c("gp", "new_cases_7d_100k", "new_tests_7d_100k", "dose1_fraction", "dose2_fraction"))
    } else {
      features <- features %>%
        mutate(gp = if_else(gp %in% subfeatures$gp, gp, as.character(NA)))
    }
    
    ps <- features %>%
      select("subject_id", "group") %>%
      bind_cols(
        predict.glm(glmResult, features, "response") %>%
          as_tibble() %>%
          rename("ps" = "value")
      ) %>%
      filter(!is.na(ps))
    
    save(
      ps, 
      selectedLassoFeatures, 
      countsFeatures, 
      glmResult, 
      file = paste0(psFolder, "/ps_model_", k, ".RData")
    )
    
    psCoefficients <- psCoefficients %>%
      union_all(
        coefficients(glmResult) %>% 
          as_tibble(rownames = "covariate") %>%
          mutate(comparison_id = k) %>%
          relocate("comparison_id")
      )
    
    cat(paste0("Comparison ", k, " finished.\n"))
    
  }
  
}

write_csv(
  comparisons %>% mutate(cdm_name = cdmName(cdm)),
  file = here(results, paste0("comparisons_", cdmName(cdm), ".csv"))
)

write_csv(
  psCoefficients %>% mutate(cdm_name = cdmName(cdm)),
  file = here(results, paste0("ps_coefficients_", cdmName(cdm), ".csv"))
)

