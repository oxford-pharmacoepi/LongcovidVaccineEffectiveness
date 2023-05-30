indexCohortSet <- cohortSet(cdm[[indexCohortName]]) %>%
  collect()
eci <- c(1:3, 2)
cci <- c(6:8, 3)
comparisons <- tibble(
  exposure_cohort_id = c(eci, eci+10, eci+20, eci+30),
  comparator_cohort_id = c(cci, cci+10, cci+20, cci+30)
) %>%
  inner_join(
    indexCohortSet %>% 
      select(
        "exposure_cohort_id" = "cohort_definition_id", 
        "exposure_name" = "cohort_name"
      ),
    by = "exposure_cohort_id"
  ) %>%
  inner_join(
    indexCohortSet %>% 
      select(
        "comparator_cohort_id" = "cohort_definition_id", 
        "comparator_name" = "cohort_name"
      ),
    by = "comparator_cohort_id"
  ) %>%
  mutate(study = as.numeric(substr(exposure_name, 7, 7))) %>%
  mutate(comparator_name = substr(comparator_name, 10, nchar(comparator_name))) %>%
  mutate(exposure_name = substr(exposure_name, 10, nchar(exposure_name))) %>%
  mutate(comparison_id = row_number()) %>%
  mutate(comparison_name = paste0("Study_", study, "_", gsub(" ", "_", exposure_name), "_vs_", gsub(" ", "_", comparator_name))) %>%
  select(comparison_id, comparison_name, study, exposure_cohort_id, exposure_name, comparator_cohort_id, comparator_name)
