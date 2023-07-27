all_asmd <- NULL
ps_all <- NULL

for (k in comparisons$comparison_id) {
  if (comparisons$skip[comparisons$comparison_id == k] == 0) {
    load(paste0(psFolder, "/ps_model_", k, ".RData"))
    rm(glmResult)
    load(here(results, tempData, "features.RData"))
    load(here(results, tempData, "cohort.RData"))
    
    cohort.k <- cdm[[indexCohortName]] %>%
      filter(cohort_definition_id %in% !!c(
        comparisons$exposure_cohort_id[comparisons$comparison_id == k], 
        comparisons$comparator_cohort_id[k]
      )) %>%
      mutate(group = if_else(
        cohort_definition_id == !!comparisons$exposure_cohort_id[comparisons$comparison_id == k],
        "exposure",
        "comparator"
      )) %>% 
      select("group", "subject_id", "cohort_start_date") %>%
      collect() %>%
      inner_join(cohort, by = c("subject_id", "cohort_start_date"))
    
    colnames.cohort <- colnames(cohort)
    rm(cohort)
    
    features <- features %>%
      inner_join(
        cohort.k %>% 
          select("subject_id", "group", "cohort_start_date"), 
        by = c("subject_id", "cohort_start_date")
      ) %>%
      select(-"cohort_start_date") %>%
      inner_join(countsFeatures, by = "feature")
    
    splitEvery <- 500
    binaryVariables <-  split(
      countsFeatures$feature, 
      ceiling(seq_along(countsFeatures$feature)/splitEvery)
    )
    featuresGroup <- binaryVariables %>%
      lapply(function(x) {
        dplyr::tibble(feature = x)
      })
    binaryVariables[[length(binaryVariables)]] <- c(binaryVariables[[length(binaryVariables)]], "sex")
    
    asmd_comparison <- NULL
    for (i in seq_along(featuresGroup)) {
      x <- features %>%
        inner_join(
          featuresGroup[[i]], 
          by = "feature"
        ) %>%
        mutate(value = 1) %>%
        pivot_wider(names_from = "feature", values_from = "value", values_fill = 0) %>%
        right_join(
          cohort.k, 
          by = c("subject_id", "group")
        ) %>%
        mutate(across(starts_with("f"), ~ if_else(is.na(.x), 0, .x))) %>%
        mutate(cohort_start_date = as.numeric(cohort_start_date))
      # unadjusted asmd
      asmd_comparison <- asmd_comparison %>%
        union_all(
          asmdBinary(x, binaryVariables[[i]]) %>%
            mutate(adjustment = "asmd_unadjusted")
        )
      # add weights
      x <- x %>% 
        inner_join(ps, by = c("subject_id", "group")) %>%
        mutate(ps = if_else(group == "exposure", 1-ps, ps)) %>%
        filter(!is.na(ps))
      # adjusted asmd
      asmd_comparison <- asmd_comparison %>%
        union_all(
          asmdBinary(x, binaryVariables[[i]], weight = "ps") %>%
            mutate(adjustment = "asmd_adjusted")
        )
    }
    
    categoricalVariables <- c("gp", "region", "age_group")
    continuousVariables <- colnames.cohort[!(
      colnames.cohort %in% c(
        "subject_id", "gp", "region", "sex", "age_group", "cohort_start_date"
      )
    )]
    
    asmd_comparison <- asmd_comparison %>%
      union_all(
        asmdContinuous(cohort.k, continuousVariables) %>%
          mutate(adjustment = "asmd_unadjusted")
      ) %>%
      union_all(
        asmdCategorical(cohort.k, categoricalVariables) %>%
          mutate(adjustment = "asmd_unadjusted")
      )
    
    cohort.k <- cohort.k %>% 
      inner_join(ps, by = c("subject_id", "group")) %>%
      mutate(ps = if_else(group == "exposure", 1-ps, ps))
    
    asmd_comparison <- asmd_comparison %>%
      union_all(
        asmdContinuous(cohort.k, continuousVariables, weight = "ps") %>%
          mutate(adjustment = "asmd_adjusted")
      ) %>%
      union_all(
        asmdCategorical(cohort.k, categoricalVariables, weight = "ps") %>%
          mutate(adjustment = "asmd_adjusted")
      ) %>%
      pivot_wider(names_from = "adjustment", values_from = "asmd") %>%
      mutate(comparison_id = comparisons$comparison_id[k])
    
    old_names <- cdm$concept %>%
      select("concept_id", "concept_name") %>%
      inner_join(
        asmd_comparison %>% 
          select("variable") %>%
          filter(substr(variable,1,1) == "f") %>%
          separate(col = "variable", sep = "_", into = c("concept_id", "window"), remove = FALSE) %>%
          mutate(concept_id = as.numeric(substr(.data$concept_id, 2, nchar(.data$concept_id)))),
        by = "concept_id",
        copy = TRUE
      ) %>%
      collect() %>%
      mutate(window = if_else(window == 1, "-30;-1", if_else(window == "2", "-365;-31", "-Any;-366"))) %>%
      mutate(new_name = paste0(concept_name, "(", concept_id, ")_window:", window)) %>%
      select(variable, new_name)
    
    asmd_comparison <- asmd_comparison %>%
      left_join(old_names, by = "variable") %>%
      mutate(variable = if_else(is.na(new_name), variable, new_name)) %>%
      select(-"new_name")
    
    all_asmd <- all_asmd %>% union_all(asmd_comparison)
    ps_all <- ps_all %>%
      union_all(
        cohort.k %>%
          select("group", "ps") %>%
          mutate(ps = if_else(group == "exposure", 1-ps, ps)) %>%
          mutate(comparison_id = comparisons$comparison_id[k])
      )
      
  }
}

write_csv(
  all_asmd %>% mutate(cdm_name = cdmName(cdm)),
  file = here(results, paste0("asmd_", cdmName(cdm), ".csv"))
)

ps_all <- ps_all %>%
  mutate(ps = round(ps, 2)) %>%
  group_by(.data$comparison_id, .data$group, .data$ps) %>%
  summarise(n = as.numeric(n()), .groups = "drop")

write_csv(
  ps_all %>% mutate(cdm_name = cdmName(cdm)),
  file = here(results, paste0("ps_distribution_", cdmName(cdm), ".csv"))
)
