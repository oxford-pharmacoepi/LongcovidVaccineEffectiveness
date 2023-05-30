ncoCohortSet <- read.csv(here("Data", "NCO.csv")) %>%
  mutate(cohort_definition_id = row_number()) %>%
  select("cohort_definition_id", "cohort_name" = "OutcomeName", "concept_id" = "ConceptId") %>%
  as_tibble()

ncoCohort <- cdm$condition_occurrence %>%
  inner_join(
    ncoCohortSet %>% 
      select(
        "cohort_definition_id", 
        "condition_concept_id" = "concept_id"
      ),
    by = "condition_concept_id",
    copy = TRUE
  ) %>%
  select(
    "cohort_definition_id",
    "subject_id" = "person_id",
    "cohort_start_date" = "condition_start_date"
  ) %>%
  mutate(cohort_end_date = cohort_start_date) %>%
  compute()

ncoCohortCount <- getCohortCount(ncoCohort)

ncoAttrition <- ncoCohortCount %>%
  mutate(
    reason_id = 1,
    reason = "Qualifying initial records",
    excluded_records = 0,
    excluded_subjects = 0
  )

cdm[[ncoCohortName]] <- newGeneratedCohortSet(
  cohortRef = computeQuery(ncoCohort, ncoCohortName, FALSE, attr(cdm, "write_schema"), TRUE),
  cohortSetRef = insertTable(ncoCohortSet, cdm, paste0(ncoCohortName, "_set")),
  cohortAttritionRef = insertTable(ncoAttrition, cdm, paste0(ncoCohortName, "_attrition")),
  cohortCountRef = insertTable(ncoCohortCount, cdm, paste0(ncoCohortName, "_count"))
)
