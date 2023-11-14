library(dplyr)
library(tidyr)
library(EmpiricalCalibration)
library(readr)
library(survival)
library(ggplot2)
library(ggfortify)
library(here)
library(meta)
library(plotly)
library(survminer)

# functions ----
clean_columns <- function(x, names) {
  for (nam in names) {
    x <- x %>%
      mutate(!!nam := gsub("_", " ", .data[[nam]])) %>%
      mutate(!!nam := paste0(toupper(substr(.data[[nam]], 1, 1)), substr(.data[[nam]], 2, nchar(.data[[nam]]))))
  }
  return(x)
}
nice <- function(x, dec = 0) {
  trimws(format(round(as.numeric(x), dec), big.mark = ",", nsmall = dec, scientific = FALSE))
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
survival_plot <- readFiles("survival_plot")
comparisons <- readFiles("comparison")

# plots ----

outcome <- "longcovid_any_symptom_90_365"
censor <- c("leave", "leave+vaccine")
cids <- c(1, 9, 17, 25)
db <- c("AURUM", "CORIVA", "GOLD", "SIDIAP")

x <- survival_plot %>%
  filter(
    outcome_name %in% outcome & censoring_method %in% censor &
      comparison_id %in% cids & cdm_name %in% db  
  ) %>%
  mutate(status = if_else(event == outcome_name, 1, 0)) %>%
  group_by(outcome_name, censoring_method, cdm_name, comparison_id, group, time, status) %>%
  summarise(weight = sum(count), .groups = "drop")

for (out in outcome) {
  for (cen in censor) {
    for (comp in cids) {
      for (dbName in db) {
        xx <- x %>%
          filter(
            outcome_name == out & censoring_method == cen &
              comparison_id == comp & cdm_name == dbName  
          )
        s <- comparisons %>% 
          filter(
            comparison_id == comp & grepl("any vaccine", exposure_name) & 
              grepl("any vaccine", comparator_name) & cdm_name == dbName
          ) %>% 
          pull("study")
        km <- survfit(Surv(time, status) ~ group, data = xx, weights = weight)
        p <- ggsurvplot(km, fun = "cloglog")
        ggsave(here("figures", "loglog", paste0(dbName, "_Study", s, "_", cen, "_", out, ".png")))
      }
    }
  }
}

