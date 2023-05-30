immunisation <- tbl(db, sql("SELECT * FROM source.immunisation"))
lookup <- tbl(db, sql("SELECT * FROM source.lookup"))
lookuptype <- tbl(db, sql("SELECT * FROM source.lookuptype"))
medical <- tbl(db, sql("SELECT * FROM source.medical"))

subset_immunisation <- immunisation %>% 
  inner_join(
    medical %>%
      filter(read_code %in% c("65F0100", "65F0200", "65F0.00", "9Niq200")) %>%
      select("medcode"),
    by = "medcode"
  ) %>%
  left_join(
    lookup %>%
      filter(lookup_type_id == 39) %>%
      select("immstype" = "code", "vax_type" = "text"),
    by = "immstype"
  ) %>%
  select("subject_id" = "patid", "cohort_start_date" = "eventdate", "vax_type") %>%
  mutate("cohort_end_date" = cohort_start_date) %>%
  compute()

vaccinated <- subset_immunisation %>%
  select(-"vax_type") %>%
  mutate(cohort_definition_id = 1) %>%
  union_all(
    subset_immunisation %>%
      inner_join(
        tibble(
          vax_type = c("COVOXFORD", "COVJANSSEN", "COVMODERNA", "COVPFIZER"),
          cohort_definition_id = c(2, 3, 4, 5)
        ),
        by = "vax_type",
        copy = TRUE
      ) %>%
      select(-"vax_type")
  )  %>%
  distinct() %>%
  select("cohort_definition_id", "subject_id", "cohort_start_date", "cohort_end_date") %>%
  compute()

vaccinatedCohortSet <- tibble(
  cohort_definition_id = c(1, 2, 3, 4, 5),
  cohort_name = c("any_vaccine", "astrazeneca", "janssen", "moderna", "pfizer")
)

vaccinatedCount <- getCohortCount(vaccinated)

vaccinatedAttrition <- vaccinatedCount %>%
  mutate(
    reason_id = 1,
    reason = "Qualifying initial records",
    excluded_records = 0,
    excluded_subjects = 0
  )

cdm[[vaccinatedCohortName]] <- newGeneratedCohortSet(
  cohortRef = computeQuery(vaccinated, vaccinatedCohortName, FALSE, attr(cdm, "write_schema"), TRUE),
  cohortSetRef = insertTable(vaccinatedCohortSet, cdm, paste0(vaccinatedCohortName, "_set")),
  cohortAttritionRef = insertTable(vaccinatedAttrition, cdm, paste0(vaccinatedCohortName, "_attrition")),
  cohortCountRef = insertTable(vaccinatedCount, cdm, paste0(vaccinatedCohortName, "_count"))
)
