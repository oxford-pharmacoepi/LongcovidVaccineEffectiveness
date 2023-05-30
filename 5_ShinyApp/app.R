library(shiny)
library(shinydashboard)
library(shinythemes)
library(shinyWidgets)
library(shinycssloaders)
library(DT)
library(plotly)
library(here)
library(scales)
library(dplyr)
library(tibble)
library(stringr)
library(tidyr)
library(EmpiricalCalibration)
library(readr)
library(cmprsk)
library(survival)
library(ggplot2)
library(ggfortify)

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
asmd <- readFiles("asmd")
cdm_snapshot <- readFiles("cdm_snapshot")
cohort_details <- readFiles("cohort_details")
comparisons <- readFiles("comparison")
nco <- readFiles("negative_control_outcomes")
estimates <- readFiles("outcome_estimates")
ps_distribution <- readFiles("ps_distribution")
ps_coefficients <- readFiles("ps_coefficients")
study_attrition <- readFiles("study_attrition")
table_characteristics <- readFiles("table_characteristics_crude")
table_characteristics_weighted <- readFiles("table_characteristics_weighted")
survival_plot <- readFiles("survival_plot")

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
cohort_details <- cohort_details %>%
  filter(cohort_table_name %in% c(
    "covid", "vaccinated", "symptoms", "general_conditions", 
    "database_specific", "longcovid", "negative_control_outcomes",
    "index"
  )) %>%
  filter(
    !(cohort_table_name %in% c("index", "vaccinated")) | 
      grepl("any|astrazeneca|pfizer", cohort_name)
  )
cdm_snapshot <- cdm_snapshot %>%
  select(-c("cdm_schema", "write_schema")) %>%
  mutate(person_cnt = nice(person_cnt, T)) %>%
  mutate(observation_period_cnt = nice(observation_period_cnt, T)) %>%
  relocate("cdm_name")
cohort_count <- cohort_details %>%
  group_by(cdm_name, cohort_table_name, cohort_name) %>%
  filter(reason_id == max(reason_id)) %>%
  ungroup() %>%
  select(cdm_name, cohort_table_name, cohort_name, number_records, number_subjects) %>%
  mutate(cohort_table_name = gsub("_cohort", "", cohort_table_name)) %>%
  mutate(cohort_table_name = gsub("longcovid_project_", "", cohort_table_name)) %>%
  clean_columns(c("cohort_table_name", "cohort_name"))
cohort_attrition <- cohort_details %>%
  clean_columns(c("cohort_table_name", "cohort_name")) %>%
  mutate(number_records = nice(number_records, T)) %>%
  mutate(number_subjects = nice(number_subjects, T)) %>%
  mutate(excluded_records = nice(excluded_records, T)) %>%
  mutate(excluded_subjects = nice(excluded_subjects, T))
study_attrition <- study_attrition %>%
  mutate(excluded = if_else(is.na(excluded), 0, excluded)) %>%
  mutate(number_observations = nice(number_observations, T)) %>%
  mutate(excluded = nice(excluded, T))
nco <- nco %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name"))
estimates <- estimates %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name")) %>%
  filter(outcome_group %in% c("nco", "longcovid28", "longcovid90", "sensitivity")) %>%
  filter(outcome_name != "next_covid") %>%
  mutate(outcome_name = if_else(outcome_name == "next_post_acute_covid19", "post_acute_covid19", outcome_name))
covariates <- c("anxiety", "asthma", "chronic_kidney_disease", "chronic_liver_disease", "copd", "dementia", "depressive_disorder", "diabetes", "gerd", "heart_failure", "hiv", "hypertension", "hypothyroidism", "infertility", "inflammarory_bowel_disease", "malignant_neoplastic_disease", "myocardial_infarction", "osteoporosis", "pneumonia", "rheumatoid_arthritis", "stroke", "venous_thromboembolism")
covariates <- paste0(covariates, "_count")
initial_counts <- cohort_count %>%
  filter(cohort_table_name %in% c(
    "Covid", "Vaccinated", "Symptoms", "General conditions", 
    "Database specific", "negative control outcomes")
  ) %>%
  mutate(number_records = nice(number_records, T)) %>%
  mutate(number_subjects = nice(number_subjects, T)) %>%
  pivot_wider(names_from = cdm_name, values_from = c(number_records, number_subjects))
counts_before <- cohort_attrition %>%
  filter(.data$cohort_table_name == "Index") %>%
  filter(reason_id == 1) %>%
  select("cohort_name", "number_subjects", "cdm_name") %>%
  mutate(Study = substr(cohort_name, 7, 7) %>% as.numeric()) %>%
  mutate("Vaccine brand" = gsub(" vaccinated", "", gsub(" unvaccinated", "", substr(cohort_name, 10, nchar(cohort_name))))) %>%
  mutate(exposure = if_else(grepl("unvaccinated", cohort_name), "unvaccinated", "vaccinated")) %>%
  select(-"cohort_name") %>%
  pivot_wider(names_from = exposure, values_from = number_subjects)
counts_after <- cohort_count %>%
  filter(.data$cohort_table_name == "Index") %>%
  select("cohort_name", "number_subjects", "cdm_name") %>%
  mutate(Study = substr(cohort_name, 7, 7) %>% as.numeric()) %>%
  mutate("Vaccine brand" = gsub(" vaccinated", "", gsub(" unvaccinated", "", substr(cohort_name, 10, nchar(cohort_name))))) %>%
  mutate(exposure = if_else(grepl("unvaccinated", cohort_name), "unvaccinated", "vaccinated")) %>%
  select(-"cohort_name") %>%
  pivot_wider(names_from = exposure, values_from = number_subjects)
attrition_index <- cohort_attrition %>%
  filter(cohort_table_name == "Index") %>%
  mutate(Study = substr(cohort_name, 7, 7) %>% as.numeric()) %>%
  mutate("Vaccine brand" = gsub(" vaccinated", "", gsub(" unvaccinated", "", substr(cohort_name, 10, nchar(cohort_name))))) %>%
  mutate(Exposure = if_else(grepl("unvaccinated", cohort_name), "unvaccinated", "vaccinated")) %>%
  select("cdm_name", "Study", "Vaccine brand", "Exposure", "number_subjects", "reason", "excluded_subjects")
outcome_count <- cohort_count %>%
  filter(cohort_table_name %in% c(
    "Longcovid", "Negative control outcomes"
  )) %>%
  select(-"cohort_table_name") %>%
  mutate(outcome_name = tolower(cohort_name)) %>%
  left_join(
    estimates %>%
      select("outcome_name", "outcome_group") %>%
      distinct() %>%
      mutate(outcome_name = gsub("_", " ", outcome_name)) %>%
      mutate(outcome_name = gsub(" days", " ", outcome_name)),
    by = "outcome_name"
  ) %>%
  select(-"outcome_name") %>%
  mutate(number_records = nice(number_records, T)) %>%
  mutate(number_subjects = nice(number_subjects, T)) %>%
  pivot_wider(names_from = cdm_name, values_from = c("number_records", "number_subjects"))
outcome_attrition <- cohort_attrition %>%
  filter(cohort_table_name %in% c("Longcovid", "Pacs", "Negative control outcomes")) %>%
  select(-"cohort_table_name") %>%
  mutate(outcome_name = tolower(cohort_name)) %>%
  left_join(
    estimates %>%
      select("outcome_name", "outcome_group") %>%
      distinct() %>%
      mutate(outcome_name = gsub("_", " ", outcome_name)) %>%
      mutate(outcome_name = gsub(" days", " ", outcome_name)),
    by = "outcome_name"
  ) %>%
  select("cdm_name", "cohort_name", "number_records", "number_subjects", "reason", "excluded_records", "excluded_subjects")
survival_plot <- survival_plot %>%
  inner_join(
    comparisons %>% 
      select(comparison_id, comparison_name, exposure_name, comparator_name, cdm_name),
    by = c("comparison_id", "cdm_name")
  ) %>%
  mutate(
    group = if_else(.data$group == "comparator", comparator_name, exposure_name),
    status = if_else(.data$event %in% c("censor", "death"), 0, 1)
  ) %>%
  group_by(group, time, status, censoring_method, outcome_name, comparison_name, cdm_name) %>%
  summarise(weight = sum(.data$count), .groups = "drop")

# ui shiny ----
ui <- dashboardPage(
  dashboardHeader(title = "Menu"),
  ## menu ----
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        "Databases",
        tabName = "databases"
      ),
      menuItem(
        "Cohort details",
        tabName = "cohort_details",
        menuSubItem("Initial cohort counts", tabName = "initial_counts"),
        menuSubItem("Study attrition", tabName = "study_attrition"),
        menuSubItem("Outcome cohort counts", tabName = "outcome_counts"),
        menuSubItem("Outcome cohort attrition", tabName = "outcome_attition")
      ),
      menuItem("Propensity scores", tabName = "ps_scores"),
      menuItem("Negative control outcomes", tabName = "nco"),
      menuItem("Estimates", tabName = "estimates"),
      menuItem("Kaplan-Meier plots", tabName = "km")
    )
  ),
  ## body ----
  dashboardBody(
    tabItems(
      ### databases ----
      tabItem(
        tabName = "databases", 
        h3("Details of the databases that participated in the study"),
        p("Identifier 'cdm_name' is the one used in the multiple selection panels of the shiny"),
        DTOutput("cdm_snapshot")
      ),
      ### cohort details ----
      #### initial counts ----
      tabItem(
        tabName = "initial_counts", 
        h3("Here you can see the counts for the initial cohorts"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "initial_counts_type",
            label = "Type of count",
            choices = c("number_records", "number_subjects"),
            selected = "number_records",
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "initial_counts_table",
            label = "Cohort group",
            choices = unique(initial_counts$cohort_table_name),
            selected = unique(initial_counts$cohort_table_name),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        DTOutput("initial_counts_out")
      ),
      #### Study attrition ----
      tabItem(
        tabName = "study_attrition",
        h2("Choose a database and see the attrition details and counts"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "study_attrition_db",
            label = "Choose database",
            choices = unique(study_attrition$cdm_name),
            selected = unique(study_attrition$cdm_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        h3("Study attrition"),
        DTOutput("study_attrition"),
        h3("Counts before index date sorting"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "before_study",
            label = "Study",
            choices = unique(counts_before$Study),
            selected = unique(counts_before$Study),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "before_vax",
            label = "Vaccine brand",
            choices = unique(counts_before[["Vaccine brand"]]),
            selected = unique(counts_before[["Vaccine brand"]]),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        DTOutput("counts_before_sorting"),
        h3("Attrition index date"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "attrition_study",
            label = "Study",
            choices = unique(attrition_index$Study),
            selected = unique(attrition_index$Study)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "attrition_brand",
            label = "Vaccine brand",
            choices = unique(attrition_index[["Vaccine brand"]]),
            selected = unique(attrition_index[["Vaccine brand"]])[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "attrition_exposure",
            label = "Exposure",
            choices = unique(attrition_index$Exposure),
            selected = unique(attrition_index$Exposure)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        DTOutput("attrition_index_date"),
        h4("Counts after index date sorting"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "after_study",
            label = "Study",
            choices = unique(counts_after$Study),
            selected = unique(counts_after$Study),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "after_vax",
            label = "Vaccine brand",
            choices = unique(counts_after[["Vaccine brand"]]),
            selected = unique(counts_after[["Vaccine brand"]]),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        DTOutput("counts_after_sorting")
      ),
      #### outcome cohort counts ----
      tabItem(
        tabName = "outcome_counts",
        h3("Please select the outcome"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "outcome_count_group",
            label = "Outcome group",
            choices = unique(outcome_count$outcome_group),
            selected = unique(outcome_count$outcome_group),
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "outcome_counts_type",
            label = "Counts type",
            choices = c("number_records", "number_subjects"),
            selected = "number_records",
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        DTOutput("outcome_counts_out")
      ),
      #### outcome cohort attrition ----
      tabItem(
        tabName = "outcome_attition",
        h3("Please select database and outcome to see attrition"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "outcome_attrition_db",
            label = "Database",
            choices = unique(outcome_attrition$cdm_name),
            selected = unique(outcome_attrition$cdm_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "outcome_attrition_name",
            label = "Outcome name",
            choices = unique(outcome_attrition$cohort_name),
            selected = unique(outcome_attrition$cohort_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        DTOutput("outcome_attrition")
      ),
      ###  ps ----
      tabItem(
        tabName = "ps_scores",
        h3("See propensity scores details"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "ps_cdm_name",
            label = "Database name",
            choices = unique(comparisons$cdm_name),
            selected = unique(comparisons$cdm_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        DTOutput("comparison_details"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "ps_comparison_name",
            label = "Select comparison",
            choices = unique(comparisons$comparison_name),
            selected = unique(comparisons$comparison_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Large scale balancing",
            DTOutput("tbl_lsc") %>% withSpinner()
          ),
          tabPanel(
            "Propensity scores plot",
            plotlyOutput("plot_ps", height = "800px") %>% withSpinner()
          ),
          tabPanel(
            "Propensity scores coefficients",
            checkboxInput(inputId = "gp_considered", label = "GP", value = TRUE),
            DTOutput("coef_ps") %>% withSpinner()
          ),
          tabPanel(
            "Table characteristics",
            DTOutput("tbl_characteristics") %>% withSpinner()
          ),
          tabPanel(
            "Table characteristics weighted",
            DTOutput("tbl_characteristics_weighted") %>% withSpinner()
          )
        )
      ),
      ### nco ----
      tabItem(
        tabName ="nco",
        h3("Negative control outcomes"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "nco_cdm_name",
            label = "CDM name",
            choices = unique(nco$cdm_name),
            selected = unique(nco$cdm_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "nco_comparison_selector",
            label = "Comparison selector",
            choices = unique(nco$comparison_name),
            selected = nco$comparison_name[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "nco_censoring_method",
            label = "Censor selector",
            choices = unique(nco$censoring_method),
            selected = nco$censoring_method[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        tags$hr(),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Effect size plot",
            plotlyOutput("effect_size_plot", height = "600px", width = "800px") %>% withSpinner()
          ),
          tabPanel(
            "Calibration plot",
            plotlyOutput("calibration_plot", height = "600px", width = "600px") %>% withSpinner()
          )
        )
      ),
      ### estimates ----
      tabItem(
        tabName = "estimates",
        h3("See estimates"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "estimates_cdm_name",
            label = "Database name",
            choices = unique(estimates$cdm_name),
            selected = unique(estimates$cdm_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "estimates_comparison_selector",
            label = "Comparison selector",
            choices = unique(estimates$comparison_name),
            selected = estimates$comparison_name[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "estimates_censoring_method",
            label = "Censor selector",
            choices = unique(estimates$censoring_method),
            selected = estimates$censoring_method[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "estimates_adjustment",
            label = "Adjustment",
            choices = unique(estimates$adjustment),
            selected = estimates$adjustment[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "estimates_outcome_group",
            label = "Outcome group",
            choices = unique(estimates$outcome_group),
            selected = unique(estimates$outcome_group)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = TRUE
          )
        ),
        uiOutput("outcome_selection"),
        tags$hr(),
        uiOutput("color_selection"),
        tags$hr(),
        plotlyOutput("effect_size_outcomes", height = "600px", width = "1200px") %>% withSpinner()
      ),
      ### survival ----
      tabItem(
        tabName = "km",
        h3("Kaplan Meier plots"),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "survival_db",
            label = "Database name",
            choices = unique(survival_plot$cdm_name),
            selected = unique(survival_plot$cdm_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "survival_comparison",
            label = "Study",
            choices = unique(survival_plot$comparison_name),
            selected = unique(survival_plot$comparison_name)[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "survival_outcome",
            label = "Outcome",
            choices = unique(survival_plot$outcome_name),
            selected = survival_plot$outcome_name[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        div(
          style = "display: inline-block;vertical-align:top; width: 150px;",
          pickerInput(
            inputId = "survival_censoring_method",
            label = "Censor selector",
            choices = unique(survival_plot$censoring_method),
            selected = survival_plot$censoring_method[1],
            options = list(
              `actions-box` = TRUE,
              size = 10,
              `selected-text-format` = "count > 3"
            ),
            multiple = FALSE
          )
        ),
        plotlyOutput("km_plot", height = "600px", width = "1000px") #%>% withSpinner()
      )
    )
  )
)

# server shiny ----
server <- function(input, output, session) {
  ## databases ----
  output$cdm_snapshot <- renderDataTable({
    datatable(cdm_snapshot, options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE))
  })
  ## initial counts ----
  initial_counts_table <- reactive({
    initial_counts %>%
      filter(.data$cohort_table_name %in% input$initial_counts_table) %>%
      select("Cohort group" = "cohort_table_name", "cohort_name", starts_with(input$initial_counts_type)) %>%
      rename_with( 
        ~ gsub(paste0(input$initial_counts_type, "_"), "", .x), 
        starts_with(input$initial_counts_type)
      )
  })
  output$initial_counts_out <- renderDataTable({
    datatable(initial_counts_table(), options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE))
  })
  ## study attrition ----
  get_study_attrition <- reactive({
    study_attrition %>%
      filter(cdm_name == input$study_attrition_db) %>%
      arrange(reason_id) %>%
      select("number_observations", "reason", "excluded")
  })
  output$study_attrition <- renderDataTable({
    datatable(
      get_study_attrition(), 
      options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE)
    )
  })
  get_counts_before_sorting <- reactive({
    counts_before %>%
      filter(.data$cdm_name == input$study_attrition_db) %>%
      filter(.data[["Vaccine brand"]] %in% input$before_vax) %>%
      filter(.data$Study %in% input$before_study) %>%
      select(-"cdm_name")
  })
  output$counts_before_sorting <- renderDataTable({
    datatable(
      get_counts_before_sorting(), 
      options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE)
    )
  })
  get_attrition_index <- reactive({
    attrition_index %>%
      filter(.data$cdm_name == input$study_attrition_db) %>%
      filter(.data$Study %in% input$attrition_study) %>%
      filter(.data$Exposure %in% input$attrition_exposure) %>%
      filter(.data$`Vaccine brand` %in% input$attrition_brand) %>%
      select(-"cdm_name")
  })
  output$attrition_index_date <- renderDataTable({
    datatable(
      get_attrition_index(), 
      options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE)
    )
  })
  get_counts_after_sorting <- reactive({
    counts_after %>%
      filter(.data$cdm_name == input$study_attrition_db) %>%
      filter(.data[["Vaccine brand"]] %in% input$after_vax) %>%
      filter(.data$Study %in% input$after_study) %>%
      select(-"cdm_name")
  })
  output$counts_after_sorting <- renderDataTable({
    datatable(
      get_counts_after_sorting(), 
      options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE)
    )
  })
  ## outcome cohort counts ----
  outcome_counts_table <- reactive({
    outcome_count %>%
      filter(.data$outcome_group %in% input$outcome_count_group) %>%
      select("Outcome group" = "outcome_group", "cohort_name", starts_with(input$outcome_counts_type)) %>%
      rename_with( 
        ~ gsub(paste0(input$outcome_counts_type, "_"), "", .x), 
        starts_with(input$outcome_counts_type)
      )
  })
  output$outcome_counts_out <- renderDataTable({
    datatable(outcome_counts_table(), options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE))
  })
  ## outcome attrition ----
  outcome_attrition_table <- reactive({
    outcome_attrition %>%
      filter(.data$cdm_name %in% input$outcome_attrition_db) %>%
      filter(.data$cohort_name %in% input$outcome_attrition_name) %>%
      select(-c("cdm_name", "cohort_name"))
  })
  output$outcome_attrition <- renderDataTable({
    datatable(outcome_attrition_table(), options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE))
  })
  ## comparisons ----
  get_comparison_id <- reactive({
    comparisons %>%
      select("comparison_id", "comparison_name") %>%
      distinct() %>%
      filter(comparison_name == input$ps_comparison_name) %>%
      pull("comparison_id")
  })
  get_comparisons_table <- reactive({
    comparisons %>%
      filter(cdm_name %in% input$ps_cdm_name) %>%
      select("comparison_name", "number_exposures", "number_comparators")
  })
  output$comparison_details <- renderDataTable({
    datatable(
      get_comparisons_table(), 
      options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE)
    )
  })
  ## large scale asmd ----
  get_lsc_asmd <- reactive({
    asmd %>%
      filter(cdm_name == input$ps_cdm_name) %>%
      filter(comparison_id == get_comparison_id())
  })
  output$tbl_lsc <- renderDataTable({
    get_lsc_asmd() %>%
      select(-c("comparison_id", "cdm_name"))
  })
  ## ps plot ----
  output$plot_ps <- renderPlotly({
    x <- ps_distribution %>%
      filter(cdm_name == input$ps_cdm_name) %>%
      filter(comparison_id == get_comparison_id())
    validate(
      need(nrow(x)>0, "No data for this comparison")
    )
    x %>%
      select("group", "ps", "n") %>%
      group_by(group) %>%
      mutate(n = as.numeric(smooth(.data$n))) %>%
      mutate(freq = 100*.data$n / sum(.data$n)) %>% 
      ungroup() %>%
      ggplot(aes(x = ps, ymin = 0, y = freq, ymax = freq, color = group, fill = group)) +
      geom_line() +
      geom_ribbon(alpha = 0.5)
  })
  ## ps coefficients ----
  output$coef_ps <- renderDataTable({
    x <- ps_coefficients %>%
      filter(cdm_name == input$ps_cdm_name) %>%
      filter(comparison_id == get_comparison_id()) %>%
      select("covariate", "value")
    if (input$gp_considered == FALSE) {
      x <- x %>% filter(!grepl("gp", covariate))
    }
    datatable(x)
  })
  ## table characteristics ----
  get_characteristics_table <- reactive({
    tbl_char <- table_characteristics %>%
      filter(comparison_id == get_comparison_id()) %>%
      filter(cdm_name == input$ps_cdm_name) %>%
      select(-c("comparison_id", "cdm_name"))
    validate(
      need(nrow(tbl_char)>0, "No data for this comparison")
    )
    tbl_char <- tbl_char %>%
      filter(group %in% c("exposure", "comparator")) %>%
      pivot_wider(names_from = c("variable", "estimate"), values_from = "value") %>%
      mutate(`Number individuals` = nice(number_subjects_count, T)) %>%
      mutate(`Age median [Q25-Q75]` = paste0(nice(age_median), "[", nice(age_quantile25), "-", nice(age_quantile75), "]")) %>%
      mutate(`Prior history median [Q25-Q75]` = paste0(nice(prior_history_median), "[", nice(prior_history_quantile25), "-", nice(prior_history_quantile75), "]")) %>%
      mutate(`Numebr visits median [Q25-Q75]` = paste0(nice(number_visits_median), "[", nice(number_visits_quantile25), "-", nice(number_visits_quantile75), "]")) %>%
      mutate(`Number PCRs median [Q25-Q75]` = paste0(nice(number_pcrs_median), "[", nice(number_pcrs_quantile25), "-", nice(number_pcrs_quantile75), "]")) %>%
      mutate(`Age group` = as.character(NA)) %>%
      mutate(across(starts_with("age_group"), ~ paste0(nice(.x, T), " (", nice(100 * .x / number_subjects_count), "%)"))) %>%
      mutate(across(starts_with("region:"), ~ paste0(nice(.x, T), " (", nice(100 * .x / number_subjects_count), "%)"))) %>%
      mutate(across(all_of(covariates), ~ paste0(nice(.x, T), " (", nice(100 * .x / number_subjects_count), "%)"))) %>%
      mutate(`Sex Female (%)` = paste0(nice(sex_female_count, T), " (", nice(100 * sex_female_count / number_subjects_count), "%)")) %>%
      mutate(`Region` = as.character(NA)) %>%
      mutate(GP = as.character(NA)) %>%
      mutate("distinct counts" = gp_distinct_counts) %>%
      mutate("individuals median [Q25 -Q75]" = paste0(gp_median_individuals, " [", gp_quantile25_individuals, "-", gp_quantile75_individuals, "]")) %>%
      select(
        "group", "Number individuals", "Age median [Q25-Q75]", "Age group", starts_with("age_group"), "Sex Female (%)",
        "Prior history median [Q25-Q75]", "Numebr visits median [Q25-Q75]",
        "Number PCRs median [Q25-Q75]", "Region", starts_with("region:"), "GP", "distinct counts", "individuals median [Q25 -Q75]",
        !!covariates
      ) %>%
      rename_with(~ gsub("age_group_", "", gsub("_count", "", gsub(";", " to ", .))), starts_with("age_group")) %>%
      rename_with(~ gsub("region: ", "", gsub("_count", "", .)), starts_with("region: ")) %>%
      rename_with(~ gsub("_count", "", .), all_of(covariates)) %>%
      mutate(across(.fns = as.character)) %>%
      pivot_longer(!"group", names_to = "covariate", values_to = "value") %>%
      pivot_wider(names_from = group, values_from = value) %>%
      left_join(
        tbl_char %>%
          filter(group == "comparison") %>%
          mutate(variable = if_else(variable == "sex", "Sex Female (%)", variable)) %>%
          mutate(variable = if_else(variable == "age_group", "Age group", variable)) %>%
          mutate(variable = if_else(variable == "region", "Region", variable)) %>%
          mutate(variable = if_else(variable == "gp", "GP", variable)) %>%
          mutate(variable = if_else(variable == "prior_history", "Prior history median [Q25-Q75]", variable)) %>%
          mutate(variable = if_else(variable == "age", "Age median [Q25-Q75]", variable)) %>%
          mutate(variable = if_else(variable == "number_visits", "Numebr visits median [Q25-Q75]", variable)) %>%
          mutate(variable = if_else(variable == "number_pcrs", "Number PCRs median [Q25-Q75]", variable)) %>%
          select("covariate" = "variable", "asmd" = "value"),
        by = "covariate"
      )
    tbl_char
  })
  output$tbl_characteristics <- renderDataTable({
    datatable(get_characteristics_table(), options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE))
  })
  ## table characteristics weighted ----
  get_characteristics_table_weighted <- reactive({
    tbl_char <- table_characteristics_weighted %>%
      filter(comparison_id == get_comparison_id()) %>%
      filter(cdm_name == input$ps_cdm_name) %>%
      select(-c("comparison_id", "cdm_name"))
    validate(
      need(nrow(tbl_char)>0, "No data for this comparison")
    )
    tbl_char <- tbl_char %>%
      filter(group %in% c("exposure", "comparator")) %>%
      pivot_wider(names_from = c("variable", "estimate"), values_from = "value") %>%
      mutate(`Number individuals` = nice(number_subjects_count, T)) %>%
      mutate(`Age median [Q25-Q75]` = paste0(nice(age_median), "[", nice(age_quantile25), "-", nice(age_quantile75), "]")) %>%
      mutate(`Prior history median [Q25-Q75]` = paste0(nice(prior_history_median), "[", nice(prior_history_quantile25), "-", nice(prior_history_quantile75), "]")) %>%
      mutate(`Numebr visits median [Q25-Q75]` = paste0(nice(number_visits_median), "[", nice(number_visits_quantile25), "-", nice(number_visits_quantile75), "]")) %>%
      mutate(`Number PCRs median [Q25-Q75]` = paste0(nice(number_pcrs_median), "[", nice(number_pcrs_quantile25), "-", nice(number_pcrs_quantile75), "]")) %>%
      mutate(`Age group` = as.character(NA)) %>%
      mutate(across(starts_with("age_group"), ~ paste0(nice(.x, T), " (", nice(100 * .x / number_subjects_count), "%)"))) %>%
      mutate(across(starts_with("region:"), ~ paste0(nice(.x, T), " (", nice(100 * .x / number_subjects_count), "%)"))) %>%
      mutate(across(all_of(covariates), ~ paste0(nice(.x, T), " (", nice(100 * .x / number_subjects_count), "%)"))) %>%
      mutate(`Sex Female (%)` = paste0(nice(sex_female_count, T), " (", nice(100 * sex_female_count / number_subjects_count), "%)")) %>%
      mutate(`Region` = as.character(NA)) %>%
      mutate(GP = as.character(NA)) %>%
      mutate("distinct counts" = gp_distinct_counts) %>%
      mutate("individuals median [Q25 -Q75]" = paste0(gp_median_individuals, " [", gp_quantile25_individuals, "-", gp_quantile75_individuals, "]")) %>%
      select(
        "group", "Number individuals", "Age median [Q25-Q75]", "Age group", starts_with("age_group"), "Sex Female (%)",
        "Prior history median [Q25-Q75]", "Numebr visits median [Q25-Q75]",
        "Number PCRs median [Q25-Q75]", "Region", starts_with("region:"), "GP", "distinct counts", "individuals median [Q25 -Q75]",
        !!covariates
      ) %>%
      rename_with(~ gsub("age_group_", "", gsub("_count", "", gsub(";", " to ", .))), starts_with("age_group")) %>%
      rename_with(~ gsub("region: ", "", gsub("_count", "", .)), starts_with("region: ")) %>%
      rename_with(~ gsub("_count", "", .), all_of(covariates)) %>%
      mutate(across(.fns = as.character)) %>%
      pivot_longer(!"group", names_to = "covariate", values_to = "value") %>%
      pivot_wider(names_from = group, values_from = value) %>%
      left_join(
        tbl_char %>%
          filter(group == "comparison") %>%
          mutate(variable = if_else(variable == "sex", "Sex Female (%)", variable)) %>%
          mutate(variable = if_else(variable == "age_group", "Age group", variable)) %>%
          mutate(variable = if_else(variable == "region", "Region", variable)) %>%
          mutate(variable = if_else(variable == "gp", "GP", variable)) %>%
          mutate(variable = if_else(variable == "prior_history", "Prior history median [Q25-Q75]", variable)) %>%
          mutate(variable = if_else(variable == "age", "Age median [Q25-Q75]", variable)) %>%
          mutate(variable = if_else(variable == "number_visits", "Numebr visits median [Q25-Q75]", variable)) %>%
          mutate(variable = if_else(variable == "number_pcrs", "Number PCRs median [Q25-Q75]", variable)) %>%
          select("covariate" = "variable", "asmd" = "value"),
        by = "covariate"
      )
    tbl_char
  })
  output$tbl_characteristics_weighted <- renderDataTable({
    datatable(get_characteristics_table_weighted(), options = list(lengthChange = FALSE, searching = FALSE, ordering = FALSE, paging = FALSE))
  })
  ## nco ----
  get_nco_table <- reactive({
    nco %>%
      filter(cdm_name == input$nco_cdm_name) %>%
      filter(comparison_name == input$nco_comparison_selector) %>%
      filter(censoring_method == input$nco_censoring_method)
  })
  output$effect_size_plot <- renderPlotly({
    nco_data <- get_nco_table()
    validate(
      need(nrow(nco_data)>0, "At least one outcome must be selected")
    )
    nco_data <- nco_data %>%
      arrange(hr) %>%
      mutate(id = row_number()) %>%
      mutate(group = if_else(
        lower_hr > 1,
        "positively associated",
        if_else(
          upper_hr < 1,
          "negatively associated",
          "no association"
        )
      )) %>%
      mutate(group = if_else(is.na(group), "no association", group))
    p <- ggplot(data=nco_data, aes(y=id, x=hr, xmin=lower_hr, xmax=upper_hr, group = group, color = group)) +
      geom_point() + 
      geom_errorbarh(height=.1) +
      scale_y_continuous(name = "", breaks=nco_data$id, labels=nco_data$outcome_name) +
      scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),limits = c(0.1, 10)) +
      theme_bw()
    p
  })
  output$calibration_plot <- renderPlotly({
    nco_data <- get_nco_table()
    validate(
      need(nrow(nco_data)>0, "At least one outcome must be selected")
    )
    p <- plotCiCalibrationEffect(nco_data$coef, nco_data$se_coef, 0*nco_data$coef)
    p
  })
  ## estimates ----
  get_varying_variables <- reactive({
    choices <- NULL
    if (length(input$estimates_cdm_name) > 1) {
      choices <- c(choices, "cdm_name")
    }
    if (length(input$estimates_comparison_selector) > 1) {
      choices <- c(choices, "comparison_name")
    }
    if (length(input$estimates_censoring_method) > 1) {
      choices <- c(choices, "censoring_method")
    }
    if (length(input$estimates_adjustment) > 1) {
      choices <- c(choices, "adjustment")
    }
    if (length(input$estimates_outcomes) > 1) {
      choices <- c(choices, "outcome_name")
    }
    choices
  })
  output$outcome_selection <- renderUI({
    selections <- estimates %>%
      filter(cdm_name %in% input$estimates_cdm_name) %>%
      filter(comparison_name %in% input$estimates_comparison_selector) %>%
      filter(censoring_method %in% input$estimates_censoring_method) %>%
      filter(adjustment %in% input$estimates_adjustment) %>%
      filter(outcome_group %in% input$estimates_outcome_group) %>%
      pull("outcome_name") %>%
      unique()
    div(
      style = "display: inline-block;vertical-align:top; width: 150px;",
      pickerInput(
        inputId = "estimates_outcomes",
        label = "Outcomes",
        choices = selections,
        selected = selections,
        options = list(
          `actions-box` = TRUE,
          size = 10,
          `selected-text-format` = "count > 3"
        ),
        multiple = TRUE
      )
    )
  })
  output$color_selection <- renderUI({
    selections <- get_varying_variables()
    if (length(selections) > 0) {
      div(
        style = "display: inline-block;vertical-align:top; width: 150px;",
        pickerInput(
          inputId = "color_by",
          label = "Color by",
          choices = selections,
          selected = selections[1],
          options = list(
            `actions-box` = TRUE,
            size = 10,
            `selected-text-format` = "count > 3"
          ),
          multiple = FALSE
        )
      )
    }
  })
  get_estimates <- reactive({
    estimates %>%
      filter(cdm_name %in% input$estimates_cdm_name) %>%
      filter(comparison_name %in% input$estimates_comparison_selector) %>%
      filter(censoring_method %in% input$estimates_censoring_method) %>%
      filter(adjustment %in% input$estimates_adjustment) %>%
      filter(outcome_name %in% input$estimates_outcomes)
  })
  output$effect_size_outcomes <- renderPlotly({
    estimates_data <- get_estimates() %>%
      mutate(id = row_number())
    if (!is.null(input$color_by)) {
      estimates_data <- estimates_data %>%
        mutate(group = .data[[input$color_by]])
    } else {
      estimates_data <- estimates_data %>%
        mutate(group = "all_estimates")
    }
    validate(
      need(nrow(estimates_data)>0, "At least one outcome must be selected")
    )
    p <- ggplot(data=estimates_data, aes(y=id, x=hr, xmin=lower_hr, xmax=upper_hr, group = group, color = group)) +
      geom_point() + 
      geom_errorbarh(height=.1) +
      scale_y_continuous(name = "", breaks=estimates_data$id, labels=estimates_data$outcome_name) +
      scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 10),limits = c(0.1, 10)) +
      theme_bw()
    p
  })
  ## survival ----
  get_survival <- reactive({
    survival_plot %>%
      filter(.data$cdm_name == input$survival_db) %>%
      filter(.data$outcome_name == input$survival_outcome) %>%
      filter(.data$comparison_name == input$survival_comparison) %>%
      filter(.data$censoring_method == input$survival_censoring_method)
  })
  output$km_plot <- renderPlotly({
    survival_plot_data <- get_survival()
    fit <- survfit(Surv(time, status) ~ group, data=survival_plot_data, weights = weight)
    autoplot(fit, censor = FALSE)
  })
}

# run shiny ----
shinyApp(ui, server)
