library(here)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(readr)

# functions ----
clean_columns <- function(x, names) {
  for (nam in names) {
    x <- x %>%
      mutate(!!nam := gsub("_", " ", .data[[nam]])) %>%
      mutate(!!nam := paste0(toupper(substr(.data[[nam]], 1, 1)), substr(.data[[nam]], 2, nchar(.data[[nam]]))))
  }
  return(x)
}
nice <- function(x, obscure = FALSE, dec = 0) {
  y <- trimws(format(round(as.numeric(x), dec), big.mark = ",", nsmall = dec, scientific = FALSE))
  if (obscure) {
    y[x<5 & x>0] <- "<5"
  }
  return(y)
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
  } else if (x == "cdm_snapshot") {
    colType <- cols(
      cdm_source_name = col_character(),
      cdm_version = col_character(),
      cdm_holder = col_character(),
      cdm_release_date = col_date(),
      vocabulary_version = col_character(),
      person_cnt = col_integer(),
      observation_period_cnt = col_integer(),
      cdm_schema = col_character(),
      write_schema = col_character(),
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

# prepare data ----
comparisons <- readFiles("comparison")
estimates <- readFiles("outcome_estimates")

comparison_new_name <- comparisons %>%
  select(comparison_name) %>%
  distinct() %>%
  mutate(comparison_new_name = gsub("tudy_", "", comparison_name)) %>%
  mutate(comparison_new_name = gsub("_", " ", comparison_new_name)) %>%
  mutate(comparison_new_name = gsub(" unvaccinated", "", comparison_new_name)) %>%
  mutate(comparison_new_name = gsub(" vaccinated", "", comparison_new_name)) %>%
  mutate(comparison_new_name = gsub(" vaccine", "", comparison_new_name)) %>%
  mutate(comparison_new_name = gsub(" vs ", " ", comparison_new_name)) %>%
  separate(comparison_new_name, c("study", "vax1", "vax2")) %>%
  mutate(vax2 = if_else(vax1 == vax2, "(VE)", paste0("- ", vax2, " (CVE)"))) %>%
  mutate(comparison_new_name = paste(study, vax1, vax2)) %>%
  select("comparison_name", "comparison_new_name")
comparisons <- renameComparisonName(comparisons, comparison_new_name)
# Exclude CVE SIDIAP
comparisons <- comparisons %>%
  filter(!(cdm_name %in% c("SIDIAP", "CORIVA")) | grepl("\\(VE\\)", .data$comparison_name)) %>%
  filter(grepl("any|astrazeneca|pfizer", exposure_name)) %>%
  filter(grepl("any|astrazeneca|pfizer", comparator_name))
estimates <- estimates %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name")) %>%
  filter(outcome_group %in% c("nco", "longcovid28", "longcovid90", "sensitivity")) %>%
  #filter(outcome_name != "next_covid") %>%
  mutate(outcome_name = if_else(outcome_name == "next_post_acute_covid19", "post_acute_covid19", outcome_name))


estimate <- estimates %>%
  filter(
    model == "finegray", 
    outcome_name %in% c("next_covid", "longcovid_any_symptom_90_365"),
    censoring_method == "leave",
    adjustment == "calibrated",
    grepl("any", comparison_name)
  ) %>%
  mutate(
    hr = paste0(round(hr, 2), " [", round(lower_hr, 2), "-", round(upper_hr, 2), "]")
  ) %>%
  mutate(
    outcome_name = if_else(outcome_name == "next_covid", "covid", "longcovid")
  ) %>%
  select(comparison_name, hr, outcome_name, cdm_name) %>%
  pivot_wider(names_from = "outcome_name", values_from = "hr")

write_csv(estimate, here("figures", "longcovid_covid.csv"))


estimate <- estimates %>%
  filter(
    outcome_name %in% c("longcovid_any_symptom_90_365"),
    censoring_method == "leave",
    adjustment == "calibrated",
    grepl("any", comparison_name)
  ) %>%
  mutate(
    hr = paste0(round(hr, 2), " [", round(lower_hr, 2), "-", round(upper_hr, 2), "]")
  ) %>%
  mutate(
    outcome_name = if_else(outcome_name == "next_covid", "covid", "longcovid")
  ) %>%
  select(comparison_name, hr, outcome_name, cdm_name, model) %>%
  pivot_wider(names_from = "model", values_from = "hr")

write_csv(estimate, here("figures", "fg_vs_cox.csv"))
