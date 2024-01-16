# Functions ----
nice <- function(x, dec = 0) {
  trimws(format(round(as.numeric(x), dec), big.mark = ",", nsmall = dec, scientific = FALSE))
}
nice_perc <- function(x, number, dec = 0) {
  ifelse(x >= 5, 
         paste0(nice(x), " (", nice(100 * x / number, dec), "%)"),
         "< 5")
}
nice_per10mil <- function(x, number, dec = 0) {
  ifelse(x >= 5, 
         paste0(nice(x), " (", nice(10000 * x / number, dec), ")"),
         "< 5")
}
clean_columns <- function(x, names) {
  for (nam in names) {
    x <- x %>%
      mutate(!!nam := gsub("_", " ", .data[[nam]])) %>%
      mutate(!!nam := paste0(toupper(substr(.data[[nam]], 1, 1)), substr(.data[[nam]], 2, nchar(.data[[nam]]))))
  }
  return(x)
}

readFiles <- function(x) {
  dataFiles <- list.files(here("data"), full.names = TRUE)
  dataFiles <- dataFiles[grepl(x, dataFiles)]
  if (x == "lsc_post_acute_covid19") {
    colType <- cols(
      cohort_definition_id = col_integer(),
      cohort_name = col_character(),
      table_name = col_character(),
      window_name = col_character(),
      concept_id = col_integer(),
      concept_name = col_character(),
      concept_count = col_character(),
      denominator_count = col_character(),
      overlap = col_logical(),
      concept_type = col_character(),
      cdm_name = col_character()
    )
  } else {
    colType <- NULL
  }
  result <- lapply(dataFiles, function(f){read_csv(f, show_col_types = FALSE, col_types = colType)}) %>%
    bind_rows()
  return(result)
}
renameComparisonName <- function(x, new_name) {
  x %>%
    inner_join(new_name, by = "comparison_name") %>%
    select(-"comparison_name") %>%
    rename("comparison_name" = "comparison_new_name")
}

table_one_study <- function(table_char, study_id) {
  tbl_char_temp <- table_char %>%
    filter(study == study_id) %>%
    select(-study, -comparison_id, -cdm_name) %>%
    filter(!grepl("region:", .data$variable)) %>%
    filter(group %in% c("exposure", "comparator")) %>%
    pivot_wider(names_from = c("variable", "estimate"), values_from = "value") %>%
    mutate(`Number individuals` = nice(number_subjects_count)) %>%
    mutate(`Age median [Q25-Q75]` = paste0(nice(age_median), " [", nice(age_quantile25), "-", nice(age_quantile75), "]")) %>%
    mutate(`Prior history median [Q25-Q75]` = paste0(nice(prior_history_median), " [", nice(prior_history_quantile25), "-", nice(prior_history_quantile75), "]")) %>%
    mutate(`Number visits median [Q25-Q75]` = paste0(nice(number_visits_median), " [", nice(number_visits_quantile25), "-", nice(number_visits_quantile75), "]")) %>%
    mutate(`Number PCRs median [Q25-Q75]` = paste0(nice(number_pcrs_median), "[", nice(number_pcrs_quantile25), "-", nice(number_pcrs_quantile75), "]")) %>%
    mutate(`Age group` = as.character(NA)) %>%
    mutate(across(starts_with("age_group"), ~ nice_perc(.x, number_subjects_count))) %>%
    mutate(across(all_of(covariates),  ~ nice_perc(.x, number_subjects_count))) %>%
    mutate(`Sex Female (%)` =  nice_perc(sex_female_count, number_subjects_count)) %>%
    mutate(`Region` = as.character(NA)) %>%
    mutate(GP = as.character(NA)) %>%
    mutate("Distinct counts" = gp_distinct_counts) %>%
    mutate("Individuals median [Q25 -Q75]" = paste0(gp_median_individuals, " [", gp_quantile25_individuals, "-", gp_quantile75_individuals, "]")) %>%
    select(
      "group", "Number individuals", "Age median [Q25-Q75]", "Age group", starts_with("age_group"), "Sex Female (%)",
      "Prior history median [Q25-Q75]", "Number visits median [Q25-Q75]",
      "Number PCRs median [Q25-Q75]", "Region", starts_with("region:"), "GP", "Distinct counts", "Individuals median [Q25 -Q75]",
      !!covariates
    ) %>%
    rename_with(~ gsub("age_group_", "", gsub("_count", "", gsub(";", " to ", .))), starts_with("age_group")) %>%
    rename_with(~ gsub("region: ", "", gsub("_count", "", .)), starts_with("region: ")) %>%
    rename_with(~ gsub("_count", "", .), all_of(covariates)) %>%
    mutate(across(.fns = as.character)) %>%
    {if ("100 to 150" %in% names(.)) {
      relocate(., "100 to 150", .after = "95 to 99") 
    } else .} %>%
    pivot_longer(!"group", names_to = "covariate", values_to = "value") %>%
    pivot_wider(names_from = group, values_from = value) %>%
    left_join(
      table_char %>%
        filter(group == "comparison") %>%
        filter(study == study_id) %>%
        select(-study, -comparison_id, -cdm_name) %>%
        mutate(variable = if_else(variable == "sex", "Sex Female (%)", variable)) %>%
        mutate(variable = if_else(variable == "age_group", "Age group", variable)) %>%
        mutate(variable = if_else(variable == "region", "Region", variable)) %>%
        mutate(variable = if_else(variable == "gp", "GP", variable)) %>%
        mutate(variable = if_else(variable == "prior_history", "Prior history median [Q25-Q75]", variable)) %>%
        mutate(variable = if_else(variable == "age", "Age median [Q25-Q75]", variable)) %>%
        mutate(variable = if_else(variable == "number_visit", "Number of GP visits median [Q25-Q75]", variable)) %>%
        mutate(variable = if_else(variable == "number_pcrs", "Number of PCR tests median [Q25-Q75]", variable)) %>%
        select("covariate" = "variable", "asmd" = "value"),
      by = "covariate"
    ) %>%
    filter(!is.na(comparator) | !is.na(asmd)) %>%
    mutate(covariate = gsub("_", " ", .data$covariate)) %>%
    mutate(covariate = stringr::str_to_sentence(.data$covariate)) %>%
    mutate(covariate = if_else(covariate == "Copd", "COPD", covariate)) %>%
    mutate(covariate = if_else(covariate == "Gerd", "GERD", covariate)) %>%
    mutate(covariate = if_else(covariate == "Number individuals", "Number individuals, N", covariate)) %>%
    mutate(covariate = if_else(covariate == "Sex female (%)", "Sex: Female, N(%)", covariate)) %>%
    mutate(covariate = gsub(" median", ", median", .data$covariate)) 
  
  comorbidities_out <- c("Chronic liver disease", "Hiv", "Infertility", "Inflammarory bowel disease")
  
  tbl_char_temp <- tbl_char_temp %>%
    filter(covariate != "Distinct counts") %>%
    filter(covariate != "Individuals, median [q25 -q75]") %>%
    mutate(covariate = gsub("q25-q75", "Q25-Q75", .data$covariate)) %>%
    filter(!grepl(" to ", .data$covariate)) %>%
    filter(!covariate %in% comorbidities_out) %>%
    filter(!grepl("Gp", .data$covariate)) %>%
    filter(!grepl("group", .data$covariate)) %>%
    filter(!grepl("Region", .data$covariate))
  
  ind <- which(tbl_char_temp$covariate == "Anxiety")     
  tbl_char_temp <- tbl_char_temp[1:(ind-1),] %>%
    full_join(tibble(covariate = "Comorbidities")) %>%
    union_all(tbl_char_temp[ind:nrow(tbl_char_temp),]) %>%
    mutate(covariate = if_else(covariate == "Comorbidities", "Comorbidities, N(%)", covariate)) %>%
    mutate(asmd = round(asmd, 3)) 
  
  return(tbl_char_temp)
}

getTableOne <- function(table_char, cve) {
  studies <- unique(table_char$study)
  if (table_char$cdm_name %>% unique() == "SIDIAP" & 
      table_char$exposure_name %>% unique() == "astrazeneca vaccinated") {
    studies <- c(2,3)
  }
  
  table_char <- table_char %>%
    select(-exposure_name, -comparator_name)
  
  count <- 1
  for (ii in studies) {
    if (cve) {
      if (count == 1) {
        if (ii == 1) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 1_BNT162b2" = "exposure",
                   # "Cohort 1_ChAdOx1" = "comparator",
                   "Cohort 1_ASMD" = "asmd") %>%
           { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 1_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 1_ChAdOx1" = "comparator")
            } }
        } else if (ii == 2) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 2_BNT162b2" = "exposure",
                   # "Cohort 2_ChAdOx1" = "comparator",
                   "Cohort 2_ASMD" = "asmd") %>%
            { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 2_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 2_ChAdOx1" = "comparator")
            } }
        } else if (ii == 3) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 3_BNT162b2" = "exposure",
                   # "Cohort 3_ChAdOx1" = "comparator",
                   "Cohort 3_ASMD" = "asmd") %>%
            { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 3_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 3_ChAdOx1" = "comparator")
            } }
        } else if (ii == 4) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 4_BNT162b2" = "exposure",
                   # "Cohort 4_ChAdOx1" = "comparator",
                   "Cohort 4_ASMD" = "asmd") %>%
            { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 4_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 4_ChAdOx1" = "comparator")
            } }
        }
      } else {
        if (ii == 1) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 1_BNT162b2" = "exposure",
                               # "Cohort 1_ChAdOx1" = "comparator",
                               "Cohort 1_ASMD" = "asmd") %>%
            { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 1_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 1_ChAdOx1" = "comparator")
            } },
            by = "covariate")
        } else if (ii == 2) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 2_BNT162b2" = "exposure",
                               # "Cohort 2_ChAdOx1" = "comparator",
                               "Cohort 2_ASMD" = "asmd")  %>%
                        { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
                          rename(., "Cohort 2_mRNA-1273" = "comparator")
                        } else {
                          rename(., "Cohort 2_ChAdOx1" = "comparator")
                        } },
                      by = "covariate")
        } else if (ii == 3) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 3_BNT162b2" = "exposure",
                               # "Cohort 3_ChAdOx1" = "comparator",
                               "Cohort 3_ASMD" = "asmd") %>%
            { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 3_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 3_ChAdOx1" = "comparator")
            } }, 
            by = "covariate")
          
        } else if (ii == 4) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 4_BNT162b2" = "exposure",
                               # "Cohort 4_ChAdOx1" = "comparator",
                               "Cohort 4_ASMD" = "asmd") %>%
            { if (table_char$cdm_name[1] %in% c("CORIVA", "UiO")) {
              rename(., "Cohort 4_mRNA-1273" = "comparator")
            } else {
              rename(., "Cohort 4_ChAdOx1" = "comparator")
            } },
            by = "covariate")
        }
      }
    } else {
      if (count == 1) {
        if (ii == 1) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 1_Vaccinated" = "exposure",
                   "Cohort 1_Unvaccinated" = "comparator",
                   "Cohort 1_ASMD" = "asmd")
        } else if (ii == 2) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 2_Vaccinated" = "exposure",
                   "Cohort 2_Unvaccinated" = "comparator",
                   "Cohort 2_ASMD" = "asmd")
        } else if (ii == 3) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 3_Vaccinated" = "exposure",
                   "Cohort 3_Unvaccinated" = "comparator",
                   "Cohort 3_ASMD" = "asmd")
        } else if (ii == 4) {
          tbl_char_studies <- table_one_study(table_char, ii) %>% 
            rename("Cohort 4_Vaccinated" = "exposure",
                   "Cohort 4_Unvaccinated" = "comparator",
                   "Cohort 4_ASMD" = "asmd")
        }
      } else {
        if (ii == 1) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 1_Vaccinated" = "exposure",
                               "Cohort 1_Unvaccinated" = "comparator",
                               "Cohort 1_ASMD" = "asmd"), by = "covariate")
        } else if (ii == 2) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 2_Vaccinated" = "exposure",
                               "Cohort 2_Unvaccinated" = "comparator",
                               "Cohort 2_ASMD" = "asmd"), by = "covariate")
        } else if (ii == 3) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 3_Vaccinated" = "exposure",
                               "Cohort 3_Unvaccinated" = "comparator",
                               "Cohort 3_ASMD" = "asmd"), by = "covariate")
        } else if (ii == 4) {
          tbl_char_studies <- tbl_char_studies %>%
            left_join(table_one_study(table_char, ii) %>% 
                        rename("Cohort 4_Vaccinated" = "exposure",
                               "Cohort 4_Unvaccinated" = "comparator",
                               "Cohort 4_ASMD" = "asmd"), by = "covariate")
        }
      }
    }
    
    
    count <- count + 1
  }
  
  border <- fp_border_default(width = 0.5, color = "gray")
  
  ind <- which(tbl_char_studies$covariate == "Comorbidities, N(%)")
  
  return(tbl_char_studies %>%
           rename(" " = "covariate") %>%
           flextable() %>%
           separate_header() %>%
           # vline(j = c('Cohort 1_ASMD', 'Cohort 2_ASMD', 'Cohort 3_ASMD'), 
           # border = border, part = "body") %>%
           hline(border = border) %>%
           set_table_properties(width = 1, layout = "autofit") %>%
           bold(part = "header") %>%
           fontsize(size = 9, part = "header") %>%
           fontsize(size = 8, part = "body")  %>%
           padding(padding = 3, part = "all") %>%
           align(align = "center", part = "header") %>%
           align(align = "left", part = "body") %>%
           bg(j = seq(2, ncol(tbl_char_studies), 2), bg = "#ededed") %>%
           bold(i = c(1:7, ind), j = 1))
}