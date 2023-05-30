allSubjects <- cdm[[indexCohortName]] %>%
  select("person_id" = "subject_id", "cohort_start_date") %>%
  distinct() %>%
  compute()
features <- cdm$condition_occurrence %>%
  inner_join(allSubjects, by = "person_id") %>%
  select(
    "subject_id" = "person_id",
    "cohort_start_date",
    "concept_id" = "condition_concept_id", 
    "date" = "condition_start_date"
  ) %>%
  filter(date < cohort_start_date) %>%
  mutate(dif_time = !!datediff("cohort_start_date", "date")) %>%
  mutate(window = as.integer(if_else(dif_time >= -30, 1, if_else(dif_time >= -365, 2, 3)))) %>%
  mutate(feature = paste0("f", concept_id, "_", window)) %>%
  select("subject_id", "cohort_start_date", "feature") %>%
  distinct() %>%
  union_all(
    cdm$drug_era %>%
      inner_join(allSubjects, by = "person_id") %>%
      select(
        "subject_id" = "person_id", 
        "cohort_start_date",
        "concept_id" = "drug_concept_id", 
        "date" = "drug_era_start_date"
      ) %>%
      mutate(date = as.Date(date)) %>%
      filter(date < cohort_start_date) %>%
      mutate(dif_time = !!datediff("cohort_start_date", "date")) %>%
      filter(dif_time >= -365) %>%
      mutate(window = as.integer(if_else(dif_time >= -30, 1, 2))) %>%
      mutate(feature = paste0("f", concept_id, "_", window)) %>%
      select("subject_id", "cohort_start_date", "feature") %>%
      distinct()
  ) %>%
  collect()
save(features, file = here(results, tempData, "features.RData"))
rm(features)
