# ANALYSIS RESULTS DATA TO REPORT ----
asmd <- readFiles("asmd")
cohort_details <- readFiles("cohort_details")
comparisons <- readFiles("comparison")
nco <- readFiles("negative_control_outcomes")
estimates <- readFiles("outcome_estimates")
study_attrition <- readFiles("study_attrition")
table_characteristics <- readFiles("table_characteristics_crude")
table_characteristics_weighted <- readFiles("table_characteristics_weighted")
survival_plot <- readFiles("survival_plot")
estimates_plot <- as_tibble(read.csv(here("data_manuscript.csv")))
ps_distribution <- readFiles("ps_distribution")


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
comparisons <- comparisons %>%
  filter(cdm_name != "SIDIAP" | grepl("\\(VE\\)", .data$comparison_name)) %>%
  filter(skip == 0)
study_attrition <- study_attrition %>%
  mutate(excluded = if_else(is.na(excluded), 0, excluded)) %>%
  mutate(number_observations = nice(number_observations)) %>%
  mutate(excluded = nice(excluded))
estimates <- estimates %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name")) %>%
  mutate(number_exposures = ifelse( number_exposures < 5,
                                    "< 5",
                                    nice(number_exposures)),
         number_comparators = ifelse( number_comparators < 5,
                                      "< 5",
                                      nice(number_comparators)))
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
  summarise(weight = sum(.data$count), .groups = "drop") %>%
  inner_join(
    comparisons %>%
      mutate(comp = substr(comparison_name, 4, nchar(comparison_name))) %>%
      filter(comp %in% c("any (VE)", "pfizer (VE)", "astrazeneca (VE)")) %>%
      mutate(study = paste0("Study ", study)), 
    by = c("cdm_name", "comparison_name"))  %>%
  mutate(comp = factor(comp, 
                       levels = c("any (VE)", "pfizer (VE)", "astrazeneca (VE)"), 
                       labels = c("Any vaccine", "BNT162b2", "ChAdOx1")),
         study = gsub("Study", "Cohort", .data$study)) %>% 
  filter(cdm_name != "UiO") %>%
  union_all(survival_plot %>%
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
              summarise(weight = sum(.data$count), .groups = "drop") %>%
              inner_join(
                comparisons %>%
                  mutate(comp = substr(comparison_name, 4, nchar(comparison_name))) %>%
                  filter(comp %in% c("pfizer - astrazeneca (CVE)")) %>%
                  mutate(study = paste0("Study ", study)), 
                by = c("cdm_name", "comparison_name"))  %>%
              filter(cdm_name %in% c("CORIVA", "UiO")) %>%
              mutate(study = gsub("Study", "Cohort", .data$study)))


covariates <- c("anxiety", "asthma", "chronic_kidney_disease", "chronic_liver_disease", "copd", "dementia", "depressive_disorder", "diabetes", "gerd", "heart_failure", "hiv", "hypertension", "hypothyroidism", "infertility", "inflammarory_bowel_disease", "malignant_neoplastic_disease", "myocardial_infarction", "osteoporosis", "pneumonia", "rheumatoid_arthritis", "stroke", "venous_thromboembolism")
covariates <- paste0(covariates, "_count")
asmd <- asmd %>%
  filter(!is.na(asmd_unadjusted) & !is.na(asmd_adjusted))
table_characteristics_weighted <- table_characteristics_weighted %>%
  inner_join(comparisons %>%
               select(comparison_id, cdm_name, study, exposure_name, comparator_name), 
             by = c("comparison_id", "cdm_name")) 
table_characteristics <- table_characteristics %>%
  inner_join(comparisons %>%
               select(comparison_id, cdm_name, study, exposure_name, comparator_name), 
             by = c("comparison_id", "cdm_name")) 


# DATA FOR TABLE AND FIGURES PLOTTING AND CAPTIONING ----
## Table One ----
captions <- comparisons %>%
  filter(skip != 1) %>%
  filter(grepl("(VE)", .data$comparison_name)) %>%
  select(cdm_name, exposure_name, comparator_name) %>%
  mutate(nice_exposure_name = ifelse(exposure_name == "any vaccine vaccinated",
                                     "any COVID-19 vaccine",
                                     ifelse(exposure_name == "astrazeneca vaccinated",
                                            "ChAdOx1 vaccine",
                                            ifelse(exposure_name == "pfizer vaccinated",
                                                   "BNT162b2 vaccine",
                                                   "no")))) %>%
  filter(nice_exposure_name != "no") %>%
  filter(!grepl(" vaccinated", .data$comparator_name)) %>%
  mutate(cdm_name_nice = ifelse(.data$cdm_name == "AURUM" | .data$cdm_name == "GOLD",
                            paste0("CPRD ", .data$cdm_name),
                            .data$cdm_name)) %>%
  mutate(caption_w_1 = paste0("Characteristics of weighted populations in ", cdm_name_nice),
         caption_u_1 = paste0("Characteristics of unweighted populations in ", cdm_name_nice),
         caption_2 = paste0(" database, stratified by staggered cohort and exposure status. Exposure is ", nice_exposure_name, ".")) %>%
  select(-nice_exposure_name) %>%
  distinct() %>%
  mutate(cve = FALSE) %>%
  union_all(comparisons %>%
              filter(skip != 1) %>%
              filter(grepl("(CVE)", .data$comparison_name)) %>%
              select(cdm_name, exposure_name, comparator_name) %>%
              filter(comparator_name %in% c("astrazeneca vaccinated")) %>%
              filter(exposure_name == "pfizer vaccinated") %>%
              filter(cdm_name != "SIDIAP") %>%
              mutate(cdm_name_nice = ifelse(.data$cdm_name == "AURUM" | .data$cdm_name == "GOLD",
                                            paste0("CPRD ", .data$cdm_name),
                                            .data$cdm_name)) %>%
              filter((cdm_name %in% c("AURUM", "GOLD") & comparator_name == "astrazeneca vaccinated")) %>%
              mutate(caption_w_1 = paste0("Characteristics of weighted populations in ", cdm_name_nice),
                     caption_u_1 = paste0("Characteristics of unweighted populations in ", cdm_name_nice),
                     caption_2 = " database, stratified by staggered cohort and vaccine brand.") %>%
              distinct() %>%
              mutate(cve = TRUE)
  ) %>%
  filter(!(cdm_name == "CORIVA" & exposure_name == "astrazeneca vaccinated")) %>%
  mutate(cdm_name = factor(cdm_name, levels = c("AURUM", "GOLD", "SIDIAP", "CORIVA"))) %>%
  arrange(cdm_name)

## Outcome tables ----
outcomes <- c("longcovid_any_symptom_90_365", "longcovid_any_symptom_28_365", "longcovid_post_acute_covid19_90_365", "next_post_acute_covid19")

captions_outcome <- captions %>%
  select(-caption_w_1, -caption_u_1, caption_2, -cdm_name, - cdm_name_nice) %>%
  distinct() %>%
  cross_join(expand_grid(
    outcome_name = outcomes,
    censoring_method = c("leave+vaccine", "leave")
  )) %>%
  mutate(outcome_name_nice = ifelse(outcome_name == outcomes[1],
                                    "any long COVID symptoms between 90 and 365 days after SARS-CoV-2",
                                    ifelse(outcome_name == outcomes[2],
                                           "any long COVID symptoms between 28 and 365 days after SARS-CoV-2",
                                           ifelse(outcome_name == outcomes[3],
                                                  "post-acute COVID-19 recorded between 90 and 365 days after SARS-CoV-2",
                                                  "post-acute COVID-19" ))),
         nice_exposure_name = ifelse(exposure_name == "any vaccine vaccinated",
                                     "any COVID-19 vaccine",
                                     ifelse(exposure_name == "astrazeneca vaccinated",
                                            "ChAdOx1 vaccine",
                                            ifelse(exposure_name == "pfizer vaccinated",
                                                   "BNT162b2 vaccine",
                                                   "no"))),
         censoring_name = ifelse(censoring_method == "leave",
                                 "No censoring at second vaccine dose for vaccinated group.",
                                 " "
         )) %>%
  mutate(
    caption_1 = paste0("Number of records for ", outcome_name_nice),
    caption_2 = ifelse(cve,
                          paste0(" across cohorts and databases, stratified by vaccine brand. ", censoring_name),
                          paste0(" across cohorts and databases, stratified by exposure status. Exposure is ", nice_exposure_name, ". ", censoring_name))) %>% 
  filter(!(cve & grepl("post-acute", .data$outcome_name_nice)))

border <- fp_border_default(width = 0.5, color = "gray")

## Kaplan Meier plots ----
captions_km <- captions_outcome %>%
  mutate(caption_1 = gsub("Number of records", "Kaplan-Meier plots depicting survival", .data$caption_1)) 

## Forest plot VE ----
forest_plot_caption <- expand_grid(
  outcome_name = outcomes,
  censoring = c("leave+vaccine", "leave")
) %>%
  inner_join(captions_km %>%
               select(outcome_name, outcome_name_nice) %>%
               distinct()) %>%
  mutate(censoring_name = ifelse(censoring== "leave",
                                 "No censoring at second vaccine dose for vaccinated group.",
                                 " ")) %>%
  mutate(caption_1 = paste0("Forest plots for vaccine effectiveness on preventing ", outcome_name_nice),
         caption_2 = paste0(" across cohorts and databases. ", censoring_name),
         censoring_filter = ifelse(
           censoring == "leave",
           "leave.png",
           "vaccine.png"
         )) 

forest_plot_caption_cve <- forest_plot_caption %>%
  mutate(caption_1 = gsub("for vaccine effectiveness", "for comparative effectiveness", .data$caption_1)) 
