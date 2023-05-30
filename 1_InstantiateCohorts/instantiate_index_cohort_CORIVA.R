## Database Specific arameters ----
info(logger, "LOAD DATABASE SPECIFIC STUDY PARAMETERS")

start.date1 <- as.Date("2021-01-28") # Start date of the study 1 (first vaccination included)
end.date1 <- as.Date("2021-04-20") # End date of the study 1 (last vaccination included)
start.date2 <- as.Date("2021-04-21") # Start date of the study 2 (first vaccination included)
end.date2 <- as.Date("2021-05-01") # End date of the study 2 (last vaccination included)
start.date3 <- as.Date("2021-05-02") # Start date of the study 3 (first vaccination included)
end.date3 <- as.Date("2021-06-12") # End date of the study 3 (last vaccination included)
start.date4 <- as.Date("2021-06-13") # Start date of the study 4 (first vaccination included)
end.date4 <- as.Date("2021-08-01") # End date of the study 4 (last vaccination included)
followback <- 180 # Minimum follow back needed
follow.date <- start.date1 - followback

covid_id <- cohortSet(cdm[[covidCohortName]]) %>%
  filter(cohort_name == "covid19") %>%
  pull("cohort_definition_id")

any_vaccine_id <- cohortSet(cdm[[vaccinatedCohortName]]) %>%
  filter(cohort_name == "any_vaccine") %>%
  pull("cohort_definition_id")

at_risk_id <- cohortSet(cdm[[databaseSpecificCohortName]]) %>%
  filter(grepl("atrisk", cohort_name)) %>%
  pull("cohort_definition_id")

high_risk_id <- cohortSet(cdm[[databaseSpecificCohortName]]) %>%
  filter(grepl("highrisk", cohort_name)) %>%
  pull("cohort_definition_id")


## create useful temp tables ----

# covid cohort
covid <- cdm[[covidCohortName]] %>%
  dplyr::filter(.data$cohort_definition_id == .env$covid_id) %>%
  dplyr::select(
    "person_id" = "subject_id",
    "covid_date" = "cohort_start_date"
  ) %>%
  dplyr::compute()

# vaccinated cohort
vaccinated <- cdm[[vaccinatedCohortName]] %>%
  dplyr::inner_join(
    cohortSet(cdm[[vaccinatedCohortName]]) %>%
      select("cohort_definition_id", "vax_type" = "cohort_name"),
    by = "cohort_definition_id", 
    copy = TRUE
  ) %>%
  dplyr::select(
    "person_id" = "subject_id",
    "vaccination_date" = "cohort_start_date",
    "vax_type"
  ) %>%
  dplyr::compute()
vaccinated <- vaccinated %>%
  dplyr::union(
    vaccinated %>%
      dplyr::filter(.data$vax_type != "any_vaccine") %>%
      dplyr::group_by(.data$person_id, .data$vaccination_date) %>%
      dplyr::tally() %>%
      dplyr::filter(.data$n > 1) %>%
      dplyr::ungroup() %>%
      dplyr::select(-"n") %>%
      dplyr::mutate(vax_type = "exclude")
  ) %>%
  dplyr::compute()

# first vaccine ever
first_vaccine <- vaccinated %>%
  dplyr::filter(.data$vax_type == "any_vaccine") %>%
  dplyr::group_by(.data$person_id) %>%
  dplyr::summarise(
    cohort_start_date = min(.data$vaccination_date, na.rm = TRUE)
  ) %>%
  dplyr::compute()

# double first vaccine
double_first_vaccine <- first_vaccine %>%
  dplyr::inner_join(
    vaccinated %>%
      dplyr::filter(.data$vax_type == "exclude") %>%
      dplyr::select("person_id", "cohort_start_date" = "vaccination_date"),
    by = c("person_id", "cohort_start_date")
  ) %>%
  dplyr::compute()

# high risk
high_risk <- cdm[[databaseSpecificCohortName]] %>%
  dplyr::inner_join(
    tibble::tibble(
      cohort_definition_id = .env$high_risk_id
    ),
    by = "cohort_definition_id",
    copy = TRUE
    ) %>%
  dplyr::select(
    "person_id" = "subject_id", 
    "high_risk_date" = "cohort_start_date"
    ) %>%
  dplyr::distinct() %>%
  dplyr::compute()

# at risk
at_risk <- cdm[[databaseSpecificCohortName]] %>%
  dplyr::inner_join(
    tibble::tibble(
      cohort_definition_id = .env$at_risk_id
    ),
    by = "cohort_definition_id",
    copy = TRUE
  ) %>%
  dplyr::select(
    "person_id" = "subject_id",
    "at_risk_date" = "cohort_start_date"
  ) %>%
  dplyr::distinct() %>%
  dplyr::compute()



## Obtain cohorts of the studies ----

### we start with all population
individuals <- cdm[["person"]] %>%
  dplyr::select("person_id", "year_of_birth", "care_site_id")

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- dplyr::tibble(
  number_observations = number_observations,
  reason = "Initial population",
  excluded = NA
)

## Study 1 ----
# in obs at start.date1
individuals <- individuals %>%
  filterInObservation(cdm, start.date1)

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("Not in observation at ", .env$start.date1),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# >= 19 years at 2021-01-01
individuals <- individuals %>%
  dplyr::mutate(age = 2021 - .data$year_of_birth) %>%
  dplyr::filter(.data$age >= 19) %>%
  dplyr::left_join(
    cdm$observation_period %>%
      select("person_id", "observation_period_start_date"),
    by = "person_id"
  ) %>%
  dplyr::select("person_id", "age", "observation_period_start_date") %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "< 19 in 2021",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# 180 days of observation
individuals <- individuals %>%
  dplyr::filter(.data$observation_period_start_date <= .env$follow.date) %>%
  dplyr::select("person_id", "age") %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "180 days of prior observation",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No covid-19
individuals <- individuals %>%
  dplyr::anti_join(
    covid %>%
      dplyr::filter(.data$covid_date <= .env$start.date1),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "No prior history of covid-19",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No vaccine
individuals <- individuals %>%
  dplyr::anti_join(
    vaccinated %>%
      dplyr::filter(.data$vaccination_date < .env$start.date1),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "No prior history of vaccination",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No double vaccinated before end.date1
individuals <- individuals %>%
  dplyr::anti_join(
    double_first_vaccine %>%
      dplyr::filter(.data$cohort_start_date <= .env$end.date1),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No double vaccinated before ", end.date1),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# Group 1
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(.data$age >= 80) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(group = as.character(NA)) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      is.na(.data$group) & .data$groupN == 1,
      "group1",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

# Group 2
# at risk, 70-80 years
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 70) %>%
      dplyr::inner_join(
        at_risk %>%
          dplyr::filter(.data$at_risk_date <= .env$start.date1) %>%
          dplyr::select("person_id") %>%
          dplyr::distinct(),
        by = "person_id"
      ) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(group = dplyr::if_else(
    !is.na(.data$groupN),
    "group2",
    .data$group
  )) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

# Group 3
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(.data$age >= 70) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(group = as.character(NA)) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      is.na(.data$group) & .data$groupN == 1,
      "group3",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

study1_cohort <- individuals %>%
  dplyr::filter(!is.na(.data$group)) %>%
  dplyr::select("person_id", "group") %>%
  dplyr::compute()

study1_groups <- study1_cohort %>%
  dplyr::group_by(.data$group) %>%
  dplyr::tally() %>%
  dplyr::collect()

study1_cohort <- study1_cohort %>%
  dplyr::select("person_id") %>%
  dplyr::compute()

individuals <- individuals %>%
  dplyr::anti_join(study1_cohort, by = "person_id") %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "Study 1",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

individuals <- individuals %>%
  dplyr::union(
    study1_cohort %>%
      dplyr::anti_join(
        first_vaccine %>%
          dplyr::filter(.data$cohort_start_date >= .env$start.date1) %>%
          dplyr::filter(.data$cohort_start_date <= .env$end.date1),
        by = "person_id"
      ) %>%
      dplyr::mutate(group = "unvaccinated study 1"),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "Unvaccinated study 1",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

## Study 2 ----
# in obs at start.date2
individuals <- individuals %>%
  filterInObservation(cdm, start.date2)

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("In observation at ", start.date2),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No covid-19 before start.date2
individuals <- individuals %>%
  dplyr::anti_join(
    covid %>%
      dplyr::filter(.data$covid_date <= .env$start.date2),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No covid-19 before ", start.date2),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No vaccine before start.date2
individuals <- individuals %>%
  dplyr::anti_join(
    vaccinated %>%
      dplyr::filter(.data$vaccination_date < .env$start.date2),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No vaccination before ", start.date2),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No double vaccinated before end.date2
individuals <- individuals %>%
  dplyr::anti_join(
    double_first_vaccine %>%
      dplyr::filter(.data$cohort_start_date <= .env$end.date2),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No double vaccinated before ", end.date2),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# Group 4
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 65) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group4",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

# Group 5
# high risk
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::inner_join(
        high_risk %>%
          dplyr::filter(.data$high_risk_date <= .env$start.date2) %>%
          dplyr::select("person_id") %>%
          dplyr::distinct(),
        by = "person_id"
      ) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group5",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()


study2_cohort <- individuals %>%
  dplyr::filter(!is.na(.data$group)) %>%
  dplyr::select("person_id", "group") %>%
  dplyr::compute()

study2_groups <- study2_cohort %>%
  dplyr::group_by(.data$group) %>%
  dplyr::tally() %>%
  dplyr::collect()

study2_cohort <- study2_cohort %>%
  dplyr::select("person_id") %>%
  dplyr::compute()

individuals <- individuals %>%
  dplyr::anti_join(study2_cohort, by = "person_id") %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "Study 2",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

individuals <- individuals %>%
  dplyr::union(
    study2_cohort %>%
      dplyr::anti_join(
        first_vaccine %>%
          dplyr::filter(.data$cohort_start_date >= .env$start.date2) %>%
          dplyr::filter(.data$cohort_start_date <= .env$end.date2),
        by = "person_id"
      ) %>%
      dplyr::mutate(group = "unvaccinated study 2"),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "Unvaccinated study 2",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

## Study 3 ----
# in obs at start.date3
individuals <- individuals %>%
  filterInObservation(cdm, start.date3)

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("In observation at ", start.date3),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No covid-19 before start.date3
individuals <- individuals %>%
  dplyr::anti_join(
    covid %>%
      dplyr::filter(.data$covid_date <= .env$start.date3),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No covid-19 before ", start.date3),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No vaccine before start.date3
individuals <- individuals %>%
  dplyr::anti_join(
    vaccinated %>%
      dplyr::filter(.data$vaccination_date < .env$start.date3),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No vaccination before ", start.date3),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No double vaccinated before end.date3
individuals <- individuals %>%
  dplyr::anti_join(
    double_first_vaccine %>%
      dplyr::filter(.data$cohort_start_date <= .env$end.date3),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No double vaccinated before ", end.date3),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)


# Group 6
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 60) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group6",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

# Group 7
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 50) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group7",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

# Group 8
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 40) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group8",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

study3_cohort <- individuals %>%
  dplyr::filter(!is.na(.data$group)) %>%
  dplyr::select("person_id", "group") %>%
  dplyr::compute()

study3_groups <- study3_cohort %>%
  dplyr::group_by(.data$group) %>%
  dplyr::tally() %>%
  dplyr::collect()

study3_cohort <- study3_cohort %>%
  dplyr::select("person_id") %>%
  dplyr::compute()

individuals <- individuals %>%
  dplyr::anti_join(study3_cohort, by = "person_id") %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "Study 3",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

individuals <- individuals %>%
  dplyr::union(
    study3_cohort %>%
      dplyr::anti_join(
        first_vaccine %>%
          dplyr::filter(.data$cohort_start_date >= .env$start.date3) %>%
          dplyr::filter(.data$cohort_start_date <= .env$end.date3),
        by = "person_id"
      ) %>%
      dplyr::mutate(group = "unvaccinated study 3"),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = "Unvaccinated study 3",
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

## Study 4 ----
# in obs at start.date4
individuals <- individuals %>%
  filterInObservation(cdm, start.date4)

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("In observation at ", start.date4),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No covid-19 before start.date4
individuals <- individuals %>%
  dplyr::anti_join(
    covid %>%
      dplyr::filter(.data$covid_date <= .env$start.date4),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No covid-19 before ", start.date4),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No vaccine before start.date4
individuals <- individuals %>%
  dplyr::anti_join(
    vaccinated %>%
      dplyr::filter(.data$vaccination_date < .env$start.date4),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No vaccination before ", start.date4),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# No double vaccinated before end.date4
individuals <- individuals %>%
  dplyr::anti_join(
    double_first_vaccine %>%
      dplyr::filter(.data$cohort_start_date <= .env$end.date4),
    by = "person_id"
  ) %>%
  dplyr::compute()

number_observations <- individuals %>%
  dplyr::tally() %>%
  dplyr::pull() %>%
  as.numeric()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = number_observations,
    reason = paste0("No double vaccinated before ", end.date4),
    excluded = last(attrition$number_observations) - .env$number_observations
  )
)

# Group 9
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 30) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group9",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()

# Group 10
individuals <- individuals %>%
  dplyr::left_join(
    individuals %>%
      dplyr::filter(is.na(.data$group)) %>%
      dplyr::filter(.data$age >= 19) %>%
      dplyr::mutate(groupN = 1) %>%
      dplyr::select("person_id", "groupN"),
    by = "person_id"
  ) %>%
  dplyr::mutate(
    group = dplyr::if_else(
      !is.na(.data$groupN),
      "group10",
      .data$group
    )
  ) %>%
  dplyr::select(-"groupN") %>%
  dplyr::compute()


study4_cohort <- individuals %>%
  dplyr::filter(!is.na(.data$group)) %>%
  dplyr::select("person_id", "group") %>%
  dplyr::compute()


study4_groups <- study4_cohort %>%
  dplyr::group_by(.data$group) %>%
  dplyr::tally() %>%
  dplyr::collect()

study4_cohort <- study4_cohort %>%
  dplyr::select("person_id") %>%
  dplyr::compute()

attrition <- rbind(
  attrition,
  dplyr::tibble(
    number_observations = 0,
    reason = "Study 4",
    excluded = .env$number_observations
  )
)

attritionStudies <- attrition %>%
  mutate(reason_id = row_number()) %>%
  mutate(cdm_name = cdmName(cdm))

write_csv(
  attritionStudies,
  file = here::here(results, paste0("study_attrition_", cdmName(cdm) ,".csv"))
)

vaccinatedCohortSet <- cdm[[vaccinatedCohortName]] %>%
  cohortSet() %>%
  collect()

## Study specific attrition ----
create_study_cohort <- function(study_cohort, s, start.date, end.date) {
  any_vaccine <- 10 * s - 9
  astrazeneca <- 10 * s - 8
  pfizer <- 10 * s - 7
  moderna <- 10 * s - 6
  janssen <- 10 * s - 5
  un_any_vaccine <- 10 * s - 4
  un_astrazeneca <- 10 * s - 3
  un_pfizer <- 10 * s - 2
  un_moderna <- 10 * s - 1
  un_janssen <- 10 * s
  
  study_cohort <- study_cohort %>%
    dplyr::left_join(
      first_vaccine %>%
        dplyr::filter(.data$cohort_start_date >= .env$start.date) %>%
        dplyr::filter(.data$cohort_start_date <= .env$end.date),
      by = "person_id"
    ) %>%
    dplyr::left_join(
      vaccinated %>%
        dplyr::select(
          "person_id", "cohort_start_date" = "vaccination_date", "vax_type"
        ),
      by = c("person_id", "cohort_start_date")
    ) %>%
    dplyr::compute()
  
  # unvaccinated vs vaccine specific
  for (vax_type in vaccinatedCohortSet$cohort_name) {
    study_cohort_index_dates <- study_cohort %>%
      dplyr::filter(.data$vax_type == .env$vax_type) %>%
      dplyr::pull("cohort_start_date")
    if (length(study_cohort_index_dates) > 0){
      study_cohort_index_dates <- study_cohort_index_dates %>%
        base::table() %>%
        dplyr::as_tibble() %>%
        dplyr::rename("cohort_start_date" = ".") %>%
        dplyr::mutate(cohort_start_date = as.Date(.data$cohort_start_date)) %>%
        dplyr::mutate(n = .data$n / sum(.data$n)) %>%
        dplyr::mutate(max_rand = cumsum(.data$n)) %>%
        dplyr::mutate(min_rand = .data$max_rand - .data$n) %>%
        dplyr::select("cohort_start_date", "min_rand", "max_rand") %>%
        dplyr::mutate(to_merge = 1)
      unvaccinated_id <- eval(parse(text = paste0("un_", vax_type)))
      vaccinated_id <- eval(parse(text = vax_type))
      study_cohort_uv <- study_cohort %>%
        dplyr::filter(is.na(.data$vax_type)) %>%
        dplyr::mutate(rand = dbplyr::sql("random()")) %>%
        dplyr::mutate(to_merge = 1) %>%
        dplyr::select(-"cohort_start_date") %>%
        dplyr::left_join(study_cohort_index_dates, by = "to_merge", copy = TRUE) %>%
        dplyr::filter(.data$rand >= .data$min_rand & .data$rand < .data$max_rand) %>%
        dplyr::mutate(cohort_definition_id = .env$unvaccinated_id) %>%
        dplyr::select("cohort_definition_id", "person_id", "cohort_start_date") %>%
        dplyr::union(
          study_cohort %>%
            dplyr::filter(.data$vax_type == .env$vax_type) %>%
            dplyr::mutate(cohort_definition_id = .env$vaccinated_id) %>%
            dplyr::select("cohort_definition_id", "person_id", "cohort_start_date")
        ) %>%
        dplyr::mutate(cohort_end_date = .data$cohort_start_date) %>%
        dplyr::rename("subject_id" = "person_id") %>%
        dplyr::compute()
      if (vax_type == vaccinatedCohortSet$cohort_name[1]) {
        study_cohort_all <- study_cohort_uv
      } else {
        study_cohort_all <- study_cohort_all %>%
          dplyr::union_all(study_cohort_uv)
      }
    }
  }
  
  study_cohort_all <- study_cohort_all %>%
    addEvent(cdm, "observation_period", NULL, c(NA, NA), "leave_db", "observation_period_end_date") %>%
    addEvent(cdm, "death", NULL, c(NA, NA), "death", "death_date") %>%
    addEvent(cdm, covidCohortName, covid_id, c(NA, NA), "first_covid")
  
  attrition <- study_cohort_all %>%
    dplyr::group_by(.data$cohort_definition_id) %>%
    dplyr::tally() %>%
    dplyr::collect() %>%
    dplyr::mutate(
      number_observations = as.numeric(.data$n), 
      excluded = as.numeric(NA),
      reason = paste0("Initial qualifying events"),
      step = 1
    ) %>%
    dplyr::select(-"n")
  
  study_cohort_all <- study_cohort_all %>%
    dplyr::filter(.data$leave_db > .data$cohort_start_date) %>%
    dplyr::filter(is.na(.data$death) | .data$death > .data$cohort_start_date) %>%
    dplyr::compute()
  
  attrition <- attrition %>%
    dplyr::union_all(
      study_cohort_all %>%
        dplyr::group_by(.data$cohort_definition_id) %>%
        dplyr::tally() %>%
        dplyr::collect() %>%
        dplyr::right_join(
          attrition %>% 
            dplyr::filter(.data$step == 1) %>%
            dplyr::select(
              "cohort_definition_id", 
              "previous_observations" = "number_observations"
            ),
          by = "cohort_definition_id"
        ) %>%
        dplyr::mutate(
          number_observations = dplyr::if_else(is.na(.data$n), 0, as.numeric(.data$n)), 
          excluded = .data$previous_observations - .data$number_observations,
          reason = paste0("Left before index date"),
          step = 2,
        ) %>%
        dplyr::select(-"n", -"previous_observations")
    )
  
  study_cohort_all <- study_cohort_all %>%
    dplyr::filter(is.na(.data$first_covid) | .data$first_covid > .data$cohort_start_date) %>%
    dplyr::select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") %>%
    dplyr::compute()
  
  attrition <- attrition %>%
    dplyr::union_all(
      study_cohort_all %>%
        dplyr::group_by(.data$cohort_definition_id) %>%
        dplyr::tally() %>%
        dplyr::collect() %>%
        dplyr::right_join(
          attrition %>% 
            dplyr::filter(.data$step == 2) %>%
            dplyr::select(
              "cohort_definition_id", 
              "previous_observations" = "number_observations"
            ),
          by = "cohort_definition_id"
        ) %>%
        dplyr::mutate(
          number_observations = dplyr::if_else(is.na(.data$n), 0, as.numeric(.data$n)), 
          excluded = .data$previous_observations - .data$number_observations,
          reason = paste0("Covid-19 before index date"),
          step = 3,
        ) %>%
        dplyr::select(-"n", -"previous_observations")
    )
  
  attr(study_cohort_all, "attrition") <- attrition
  attr(study_cohort_all, "cohortSet") <- dplyr::tibble(
    cohort_definition_id = c(
      any_vaccine, astrazeneca, pfizer, moderna, janssen, un_any_vaccine, 
      un_astrazeneca, un_pfizer, un_moderna, un_janssen
    ),
    cohort_name = c(
      paste0("study ", s, ": any vaccine vaccinated"), 
      paste0("study ", s, ": astrazeneca vaccinated"), 
      paste0("study ", s, ": pfizer vaccinated"),
      paste0("study ", s, ": moderna vaccinated"), 
      paste0("study ", s, ": janssen vaccinated"), 
      paste0("study ", s, ": any vaccine unvaccinated"),
      paste0("study ", s, ": astrazeneca unvaccinated"), 
      paste0("study ", s, ": pfizer unvaccinated"),
      paste0("study ", s, ": moderna unvaccinated"), 
      paste0("study ", s, ": janssen unvaccinated")
    )
  )
  
  return(study_cohort_all)
}

study1_cohort <- create_study_cohort(study1_cohort, 1, start.date1, end.date1)
study2_cohort <- create_study_cohort(study2_cohort, 2, start.date2, end.date2)
study3_cohort <- create_study_cohort(study3_cohort, 3, start.date3, end.date3)
study4_cohort <- create_study_cohort(study4_cohort, 4, start.date4, end.date4)

indexCohort <- study1_cohort %>%
  dplyr::union_all(study2_cohort) %>%
  dplyr::union_all(study3_cohort) %>%
  dplyr::union_all(study4_cohort) %>%
  dplyr::compute()

indexCohortSet <- attr(study1_cohort, "cohortSet") %>%
  dplyr::union_all(attr(study2_cohort, "cohortSet")) %>%
  dplyr::union_all(attr(study3_cohort, "cohortSet")) %>%
  dplyr::union_all(attr(study4_cohort, "cohortSet"))

indexAttrition <- attr(study1_cohort, "attrition") %>%
  dplyr::union_all(attr(study2_cohort, "attrition")) %>%
  dplyr::union_all(attr(study3_cohort, "attrition")) %>%
  dplyr::union_all(attr(study4_cohort, "attrition")) %>%
  mutate(
    number_records = number_observations,
    number_subjects = number_observations,
    excluded_records = excluded,
    excluded_subjects = excluded,
    reason_id = step
  ) %>%
  select(
    "cohort_definition_id", "number_records", "number_subjects", "reason_id",
    "reason", "excluded_records", "excluded_subjects"
  ) %>%
  arrange(cohort_definition_id, reason_id)

indexCohortCount <- getCohortCount(indexCohort)

cdm[[indexCohortName]] <- newGeneratedCohortSet(
  cohortRef = computeQuery(indexCohort, indexCohortName, FALSE, attr(cdm, "write_schema"), TRUE),
  cohortSetRef = insertTable(indexCohortSet, cdm, paste0(indexCohortName, "_set")),
  cohortAttritionRef = insertTable(indexAttrition, cdm, paste0(indexCohortName, "_attrition")),
  cohortCountRef = insertTable(indexCohortCount, cdm, paste0(indexCohortName, "_count"))
)
