cdm[["visit_occurrence"]] <- cdm[["visit_occurrence"]] %>%
  dplyr::filter(lubridate::year(.data$visit_start_date) >= 2020) %>%
  dplyr::compute()

covid_test_id <- cohortSet(cdm[[covidCohortName]]) %>%
  filter(cohort_name == "any_covid19_test") %>%
  pull("cohort_definition_id")

generalConditionsCohortSet <- cohortSet(cdm[[generalConditionsCohortName]]) %>%
  collect()

# add needed covariates
cohort <- cdm[[indexCohortName]] %>%
  select("subject_id", "cohort_start_date") %>%
  distinct() %>%
  compute() %>%
  addAge(cdm) %>%
  addAgeGroup(ageGroup = list(
    c(15, 19), c(20, 24), c(25, 29), c(30, 34), c(35, 39), c(40, 44), c(45, 49),
    c(50, 54), c(55, 59), c(60, 64), c(65, 69), c(70, 74), c(75, 79), c(80, 84),
    c(85, 89), c(90, 94), c(95, 99), c(100, 150)
  )) %>%
  addPriorHistory(cdm) %>%
  addSex(cdm) %>%
  addNumberVisit(cdm, c(-180, -1)) %>%
  addNumberEvent(cdm, covidCohortName, list(cohort_definition_id = covid_test_id), c(-180, -1), "number_pcrs") %>%
  addGPIdentifier(cdm) %>%
  PatientProfiles::addCohortIntersectFlag(cdm, generalConditionsCohortName, window = c(-Inf, 0), nameStyle = "{cohort_name}")

if (cdmName(cdm) %in% c("AURUM", "GOLD")) {
  cohort <- cohort %>%
    inner_join(
      cdm$person %>%
        select("subject_id" = "person_id", "care_site_id") %>%
        inner_join(
          cdm$care_site %>%
            select("care_site_id", "location_id") %>%
            inner_join(
              cdm$location %>%
                select("location_id", "region" = "location_source_value"),
              by = "location_id"
            ),
          by = "care_site_id"
        ) %>%
        select("subject_id", "region"),
      by = "subject_id"
    ) %>%
    mutate(region = tolower(region)) %>%
    mutate(region = if_else(
      region %in% c("south east coast", "south central"),
      "south east",
      region
    ))
} else if (cdmName(cdm) == "SIDIAP") {
  cohort <- cohort %>%
    inner_join(
      cdm$person %>%
        select("subject_id" = "person_id", "location_id") %>%
        inner_join(
          cdm$location %>%
            select("location_id", "region" = "location_source_value"),
          by = "location_id"
        ) %>%
        select(-"location_id"),
      by = "subject_id"
    )
} else {
  cohort <- cohort %>%
    mutate(region = as.character(NA))
}
cohort <- cohort %>%
  collect()

tableCharacteristics <- NULL
for (k in comparisons$comparison_id) {
  if (comparisons$skip[comparisons$comparison_id == k] == 0) {
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
      mutate(gp = as.character(gp))
    tableCharacteristics.k <- cohort.k %>%
      group_by(group, age_group) %>%
      summarise(value = as.character(n()), .groups = "drop") %>%
      mutate(
        variable = paste0("age_group_", age_group),
        estimate = "count"
      ) %>%
      select("group", "variable", "estimate", "value") %>%
      union_all(
        cohort.k %>%
          group_by(group) %>%
          summarise(
            `number-subjects_count` = n_distinct(subject_id),
            `sex-female_count` = sum(sex == "Female"),
            `age_median` = median(age),
            `age_quantile25` = quantile(age, 0.25),
            `age_quantile75` = quantile(age, 0.75),
            `prior-history_median` = median(prior_history),
            `prior-history_quantile25` = quantile(prior_history, 0.25),
            `prior-history_quantile75` = quantile(prior_history, 0.75),
            `number-visits_median` = median(number_visit),
            `number-visits_quantile25` = quantile(number_visit, 0.25),
            `number-visits_quantile75` = quantile(number_visit, 0.75),
            `number-pcrs_median` = median(number_pcrs),
            `number-pcrs_quantile25` = quantile(number_pcrs, 0.25),
            `number-pcrs_quantile75` = quantile(number_pcrs, 0.75),
            .groups = "drop"
          ) %>%
          pivot_longer(!"group", names_to = "variables", values_to = "value") %>%
          separate(col = "variables", sep = "_", into = c("variable", "estimate")) %>%
          mutate(variable = gsub("-", "_", variable)) %>%
          mutate(value = as.character(value))
      ) %>%
      union_all(
        cohort.k %>%
          group_by(group, region) %>%
          summarise(value = as.character(n()), .groups = "drop") %>%
          mutate(variable = paste0("region: ", region)) %>%
          mutate(estimate = "count") %>%
          select(-"region")
      ) %>%
      union_all(
        cohort.k %>%
          group_by(group, gp) %>%
          summarise(count = n(), .groups = "drop") %>%
          group_by(group) %>%
          summarise(
            `gp-distinct_counts` = n_distinct(gp),
            `gp-median_individuals` = median(count),
            `gp-quantile25_individuals` = quantile(count, 0.25),
            `gp-quantile75_individuals` = quantile(count, 0.75),
            .groups = "drop"
          ) %>%
          pivot_longer(!"group", names_to = c("variable", "estimate"), names_sep = "-", values_to = "value") %>%
          mutate(value = as.character(value))
      ) %>%
      union_all(
        cohort.k %>%
          group_by(group) %>%
          summarise(across(all_of(generalConditionsCohortSet$cohort_name), ~ sum(.x))) %>%
          ungroup() %>%
          pivot_longer(!"group", names_to = "variable", values_to = "value") %>%
          mutate(estimate = "count") %>%
          mutate(value = as.character(value))
      ) %>%
      union_all(
        asmdBinary(cohort.k, c("sex", generalConditionsCohortSet$cohort_name)) %>%
          mutate(estimate = "asmd") %>%
          rename("value" = "asmd") %>%
          select(-"asmd_type") %>%
          mutate(group = "comparison") %>%
          mutate(value = as.character(value))
      ) %>%
      union_all(
        asmdCategorical(cohort.k, c("age_group", "region", "gp")) %>%
          mutate(estimate = "asmd") %>%
          rename("value" = "asmd") %>%
          select(-"asmd_type") %>%
          mutate(group = "comparison") %>%
          mutate(value = as.character(value))
      ) %>%
      union_all(
        asmdContinuous(cohort.k, c("prior_history", "age", "number_visit", "number_pcrs")) %>%
          mutate(estimate = "asmd") %>%
          rename("value" = "asmd") %>%
          select(-"asmd_type") %>%
          mutate(group = "comparison") %>%
          mutate(value = as.character(value))
      ) %>%
      mutate(comparison_id = comparisons$comparison_id[comparisons$comparison_id == k])
    
    tableCharacteristics <- tableCharacteristics %>%
      union_all(tableCharacteristics.k)
  }
}

write_csv(
  tableCharacteristics %>% mutate(cdm_name = cdmName(cdm)), 
  here(results, paste0("table_characteristics_crude_", cdmName(cdm),  ".csv"))
)

tableCharacteristicsWeighted <- NULL
comparisonIds <-  comparisons$comparison_id[comparisons$skip == 0]
for (k in comparisonIds) {
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
    mutate(gp = as.character(gp))
  ## add ps
  load(paste0(psFolder, "/ps_model_", k, ".RData"))
  rm(countsFeatures, glmResult, selectedLassoFeatures)
  cohort.k <- cohort.k %>%
    inner_join(ps, by = c("subject_id", "group")) %>%
    mutate(weight = if_else(group == "exposure", 1-ps, ps))
  rm(ps)
  # get characteristics
  tableCharacteristicsWeighted.k <- cohort.k %>%
    group_by(group, age_group) %>%
    summarise(value = as.character(round(sum(.data$weight))), .groups = "drop") %>%
    mutate(
      variable = paste0("age_group_", age_group),
      estimate = "count"
    ) %>%
    select("group", "variable", "estimate", "value") %>%
    union_all(
      cohort.k %>%
        group_by(group) %>%
        summarise(
          `number-subjects_count` = round(sum(.data$weight)),
          `sex-female_count` = round(sum(.data$weight[sex == "Female"])),
          `age_median` = weighted.median(age, weight),
          `age_quantile25` = weighted.quantile(age, weight, 0.25),
          `age_quantile75` = weighted.quantile(age, weight, 0.75),
          `prior-history_median` = weighted.median(prior_history, weight),
          `prior-history_quantile25` = weighted.quantile(prior_history, weight, 0.25),
          `prior-history_quantile75` = weighted.quantile(prior_history, weight, 0.75),
          `number-visits_median` = weighted.median(number_visit, weight),
          `number-visits_quantile25` = weighted.quantile(number_visit, weight, 0.25),
          `number-visits_quantile75` = weighted.quantile(number_visit, weight, 0.75),
          `number-pcrs_median` = weighted.median(number_pcrs, weight),
          `number-pcrs_quantile25` = weighted.quantile(number_pcrs, weight, 0.25),
          `number-pcrs_quantile75` = weighted.quantile(number_pcrs, weight, 0.75),
          .groups = "drop"
        ) %>%
        pivot_longer(!"group", names_to = "variables", values_to = "value") %>%
        separate(col = "variables", sep = "_", into = c("variable", "estimate")) %>%
        mutate(variable = gsub("-", "_", variable)) %>%
        mutate(value = as.character(round(value)))
    ) %>%
    union_all(
      cohort.k %>%
        group_by(group, region) %>%
        summarise(value = as.character(round(sum(.data$weight))), .groups = "drop") %>%
        mutate(variable = paste0("region: ", region)) %>%
        mutate(estimate = "count") %>%
        select(-"region")
    ) %>%
    union_all(
      cohort.k %>%
        group_by(group, gp) %>%
        summarise(count = round(sum(.data$weight)), .groups = "drop") %>%
        group_by(group) %>%
        summarise(
          `gp-distinct_counts` = n_distinct(gp),
          `gp-median_individuals` = median(count),
          `gp-quantile25_individuals` = quantile(count, 0.25),
          `gp-quantile75_individuals` = quantile(count, 0.75),
          .groups = "drop"
        ) %>%
        pivot_longer(!"group", names_to = c("variable", "estimate"), names_sep = "-", values_to = "value") %>%
        mutate(value = as.character(value))
    ) %>%
    union_all(
      cohort.k %>%
        group_by(group) %>%
        summarise(across(all_of(generalConditionsCohortSet$cohort_name), ~ sum(.x * weight))) %>%
        ungroup() %>%
        pivot_longer(!"group", names_to = "variable", values_to = "value") %>%
        mutate(estimate = "count") %>%
        mutate(value = as.character(round(value)))
    ) %>%
    union_all(
      asmdBinary(cohort.k, c("sex", generalConditionsCohortSet$cohort_name), weight = "weight") %>%
        mutate(estimate = "asmd") %>%
        rename("value" = "asmd") %>%
        select(-"asmd_type") %>%
        mutate(group = "comparison") %>%
        mutate(value = as.character(value))
    ) %>%
    union_all(
      asmdCategorical(cohort.k, c("age_group", "region", "gp"), weight = "weight") %>%
        mutate(estimate = "asmd") %>%
        rename("value" = "asmd") %>%
        select(-"asmd_type") %>%
        mutate(group = "comparison") %>%
        mutate(value = as.character(value))
    ) %>%
    union_all(
      asmdContinuous(cohort.k, c("prior_history", "age", "number_visit", "number_pcrs"), weight = "weight") %>%
        mutate(estimate = "asmd") %>%
        rename("value" = "asmd") %>%
        select(-"asmd_type") %>%
        mutate(group = "comparison") %>%
        mutate(value = as.character(value))
    ) %>%
    mutate(comparison_id = k)
  
  tableCharacteristicsWeighted <- tableCharacteristicsWeighted %>%
    union_all(tableCharacteristicsWeighted.k)
}

write_csv(
  tableCharacteristicsWeighted %>% mutate(cdm_name = cdmName(cdm)), 
  here(results, paste0("table_characteristics_weighted_", cdmName(cdm),  ".csv"))
)
