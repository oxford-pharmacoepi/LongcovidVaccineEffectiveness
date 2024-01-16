library(dplyr)
library(tidyr)
library(EmpiricalCalibration)
library(readr)
library(survival)
library(ggplot2)
library(ggfortify)
library(here)
library(meta)
library(scales)

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
asmd <- readFiles("asmd")
cohort_details <- readFiles("cohort_details")
comparisons <- readFiles("comparison")
nco <- readFiles("negative_control_outcomes")
estimates <- readFiles("outcome_estimates")
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
comparisons <- comparisons %>%
  filter(cdm_name != "SIDIAP" | grepl("\\(VE\\)", .data$comparison_name)) %>%
  filter(skip == 0)
study_attrition <- study_attrition %>%
  mutate(excluded = if_else(is.na(excluded), 0, excluded)) %>%
  mutate(number_observations = nice(number_observations)) %>%
  mutate(excluded = nice(excluded))
nco <- nco %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name"))
estimates <- estimates %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name"))
covariates <- c("anxiety", "asthma", "chronic_kidney_disease", "chronic_liver_disease", "copd", "dementia", "depressive_disorder", "diabetes", "gerd", "heart_failure", "hiv", "hypertension", "hypothyroidism", "infertility", "inflammarory_bowel_disease", "malignant_neoplastic_disease", "myocardial_infarction", "osteoporosis", "pneumonia", "rheumatoid_arthritis", "stroke", "venous_thromboembolism")
covariates <- paste0(covariates, "_count")
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
asmd <- asmd %>%
  filter(!is.na(asmd_unadjusted) & !is.na(asmd_adjusted))

# Table 1 baseline characteristics ----
# Figure 1 asmd before vs after ----
## GOLD ----
data_plot <- asmd %>%
  inner_join(
    comparisons %>%
      filter(cdm_name == "GOLD") %>%
      mutate(comp = substr(comparison_name, 4, nchar(comparison_name))) %>%
      filter(comp %in% c("any (VE)", "pfizer (VE)", "astrazeneca (VE)", "pfizer - astrazeneca (CVE)")) %>%
      mutate(comp = factor(comp, c("any (VE)", "pfizer (VE)", "astrazeneca (VE)", "pfizer - astrazeneca (CVE)"))) %>%
      mutate(study = paste0("Cohort ", study)), 
    by = c("cdm_name", "comparison_id")
  ) %>%
  group_by(comparison_name) %>%
  mutate(
    number_unadjusted = paste0(
      sum(asmd_unadjusted > 0.1), 
      " (",
      round(100*sum(asmd_unadjusted > 0.1)/n()),
      "%)"
    )
  ) %>%
  ungroup()
unbalanced <- data_plot %>%
  filter(asmd_adjusted > 0.1)
data_label <- data_plot %>%
  select(study, comp, number_unadjusted) %>%
  distinct()
ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
  geom_point(alpha = 0.2, shape = 21, colour = "blue", fill = "blue", size = 1.5) +
  geom_hline(yintercept=0.1) +
  geom_vline(xintercept=0.1) +
  facet_grid(study ~ comp) + 
  geom_point(data = unbalanced, aes(x = asmd_unadjusted, y = asmd_adjusted), colour = "black", fill = NA, size = 2) +
  geom_text(data = unbalanced, aes(label=variable, x = asmd_unadjusted, y = asmd_adjusted), nudge_x = 0.1, colour = "black") +
  geom_label(data = data_label, aes(label=number_unadjusted, x=1.5, y =0.05), nudge_x = 0.05) + 
  xlab("ASMD before OW") +
  ylab("ASMD after OW")
ggsave(
  "Figure1_GOLD.png",
  plot = last_plot(),
  path = here("figures"),
  scale = 1,
  width = 3000,
  height = 2000,
  units = "px",
  dpi = 300
)
## AURUM ----
data_plot <- asmd %>%
  inner_join(
    comparisons %>%
      filter(cdm_name == "AURUM") %>%
      mutate(comp = substr(comparison_name, 4, nchar(comparison_name))) %>%
      filter(comp %in% c("any (VE)", "pfizer (VE)", "astrazeneca (VE)", "pfizer - astrazeneca (CVE)")) %>%
      mutate(comp = factor(comp, c("any (VE)", "pfizer (VE)", "astrazeneca (VE)", "pfizer - astrazeneca (CVE)"))) %>%
      mutate(study = paste0("Cohort ", study)), 
    by = c("cdm_name", "comparison_id")
  ) %>%
  group_by(comparison_name) %>%
  mutate(
    number_unadjusted = paste0(
      sum(asmd_unadjusted > 0.1), 
      " (",
      round(100*sum(asmd_unadjusted > 0.1)/n()),
      "%)"
    )
  ) %>%
  ungroup()
unbalanced <- data_plot %>%
  filter(asmd_adjusted > 0.1)
data_label <- data_plot %>%
  select(study, comp, number_unadjusted) %>%
  distinct()
ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
  geom_point(alpha = 0.2, shape = 21, colour = "blue", fill = "blue", size = 1.5) +
  geom_hline(yintercept=0.1) +
  geom_vline(xintercept=0.1) +
  facet_grid(study ~ comp) + 
  geom_point(data = unbalanced, aes(x = asmd_unadjusted, y = asmd_adjusted), colour = "black", fill = NA, size = 2) +
  geom_text(data = unbalanced, aes(label=variable, x = asmd_unadjusted, y = asmd_adjusted), nudge_x = 0.1, colour = "black") +
  geom_label(data = data_label, aes(label=number_unadjusted, x=1.3, y =0.145)) + 
  xlab("ASMD before OW") +
  ylab("ASMD after OW")
ggsave(
  "Figure1_AURUM.png",
  plot = last_plot(),
  path = here("figures"),
  scale = 1,
  width = 3000,
  height = 2000,
  units = "px",
  dpi = 300
)
## SIDIAP ----
data_plot <- asmd %>%
  inner_join(
    comparisons %>%
      filter(cdm_name == "SIDIAP") %>%
      mutate(comp = substr(comparison_name, 4, nchar(comparison_name))) %>%
      filter(comp %in% c("any (VE)", "pfizer (VE)", "astrazeneca (VE)")) %>%
      mutate(comp = factor(comp, c("any (VE)", "pfizer (VE)", "astrazeneca (VE)"))) %>%
      mutate(study = paste0("Cohort ", study)), 
    by = c("cdm_name", "comparison_id")
  ) %>%
  group_by(comparison_name) %>%
  mutate(
    number_unadjusted = paste0(
      sum(asmd_unadjusted > 0.1), 
      " (",
      round(100*sum(asmd_unadjusted > 0.1)/n()),
      "%)"
    )
  ) %>%
  ungroup()
unbalanced <- data_plot %>%
  filter(asmd_adjusted > 0.1)
data_label <- data_plot %>%
  select(study, comp, number_unadjusted) %>%
  distinct()
ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
  geom_point(alpha = 0.2, shape = 21, colour = "blue", fill = "blue", size = 1.5) +
  geom_hline(yintercept=0.1) +
  geom_vline(xintercept=0.1) +
  facet_grid(study ~ comp) + 
  geom_point(data = unbalanced, aes(x = asmd_unadjusted, y = asmd_adjusted), colour = "black", fill = NA, size = 2) +
  geom_text(data = unbalanced, aes(label=variable, x = asmd_unadjusted, y = asmd_adjusted), nudge_x = 0.1, colour = "black") +
  geom_label(data = data_label, aes(label=number_unadjusted, x=1.3, y =0.145)) + 
  xlab("ASMD before OW") +
  ylab("ASMD after OW")
ggsave(
  "Figure1_SIDIAP.png",
  plot = last_plot(),
  path = here("figures"),
  scale = 1,
  width = 3000,
  height = 2000,
  units = "px",
  dpi = 300
)

## CORIVA ----
data_plot <- asmd %>%
  inner_join(
    comparisons %>%
      filter(cdm_name == "CORIVA") %>%
      mutate(comp = substr(comparison_name, 4, nchar(comparison_name))) %>%
      filter(comp %in% c("any (VE)", "pfizer (VE)")) %>%
      mutate(comp = factor(comp, c("any (VE)", "pfizer (VE)"))) %>%
      mutate(study = paste0("Cohort ", study)), 
    by = c("cdm_name", "comparison_id")
  ) %>%
  group_by(comparison_name) %>%
  mutate(
    number_unadjusted = paste0(
      sum(asmd_unadjusted > 0.1), 
      " (",
      round(100*sum(asmd_unadjusted > 0.1)/n()),
      "%)"
    )
  ) %>%
  ungroup()
unbalanced <- data_plot %>%
  filter(asmd_adjusted > 0.1)
data_label <- data_plot %>%
  select(study, comp, number_unadjusted) %>%
  distinct()
ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
  geom_point(alpha = 0.2, shape = 21, colour = "blue", fill = "blue", size = 1.5) +
  geom_hline(yintercept=0.1) +
  geom_vline(xintercept=0.1) +
  facet_grid(study ~ comp) + 
  geom_point(data = unbalanced, aes(x = asmd_unadjusted, y = asmd_adjusted), colour = "black", fill = NA, size = 2) +
  geom_text(data = unbalanced, aes(label=variable, x = asmd_unadjusted, y = asmd_adjusted), nudge_x = 0.1, colour = "black") +
  geom_label(data = data_label, aes(label=number_unadjusted, x=0.25, y =0.05), nudge_x = 0.05) + 
  xlab("ASMD before OW") +
  ylab("ASMD after OW")
ggsave(
  "Figure1_CORIVA.png",
  plot = last_plot(),
  path = here("figures"),
  scale = 1,
  width = 3000,
  height = 2000,
  units = "px",
  dpi = 300
)
# Figure 2 vaccine effectiveness ----
modelMeta <- "random"
outcomes <- c(
  "longcovid_post_acute_covid19_28_365", "longcovid_post_acute_covid19_90_365", 
  "longcovid_any_symptom_28_365", "longcovid_any_symptom_90_365", 
  "next_post_acute_covid19"
)
estimates_plot <- estimates %>%
  mutate(comparison = substr(comparison_name, 4, nchar(comparison_name))) %>%
  filter(variable == "groupexposure") %>%
  filter(skip == 0) %>%
  filter(outcome_name %in% outcomes) %>%
  mutate(study = paste0("Study ", study)) %>%
  select(
    "cdm_name", "comparison", "adjustment", "censoring_method", "outcome_name", 
    "study", "coef", "se_coef", "hr", "lower_hr", "upper_hr"
  )
cdmName <- unique(estimates_plot$cdm_name)
comp <- unique(estimates_plot$comparison)
cens <- unique(estimates_plot$censoring_method)
ajust <- unique(estimates_plot$adjustment)
estimates_metaanalysis <- NULL
for (cdmNameK in cdmName) {
  for (compK in comp) {
    for (censK in cens) {
      for (ajustK in ajust) {
        for (outcome in outcomes) {
          x <- estimates_plot %>%
            filter(
              cdm_name == cdmNameK & comparison == compK & 
                censoring_method == censK & adjustment == ajustK &
                outcome_name == outcome
            )
          if (nrow(x) > 1) {
            xMA <- metagen(TE = x$coef, seTE = x$se_coef, sm = "HR")
            resultMeta <- tibble(
              cdm_name = cdmNameK, comparison = compK, adjustment = ajustK,
              censoring_method = censK, outcome_name = outcome,
              study = "Meta Analysis", 
              coef = eval(parse(text = paste0("xMA$TE.", modelMeta))),
              se_coef = eval(parse(text = paste0("xMA$seTE.", modelMeta))),
              hr = exp(eval(parse(text = paste0("xMA$TE.", modelMeta)))),
              lower_hr = exp(eval(parse(text = paste0("xMA$lower.", modelMeta)))),
              upper_hr = exp(eval(parse(text = paste0("xMA$upper.", modelMeta)))),
              heterogeneity = xMA$I2
            )
            estimates_metaanalysis <- bind_rows(estimates_metaanalysis, resultMeta)
          }
        }
      }
    }
  }
}
estimates_plot <- estimates_plot %>%
  bind_rows(estimates_metaanalysis)

write_csv(estimates_plot, here("data_manuscript.csv"))

for (censK in cens) {
  for (outcome in outcomes) {
    data <- estimates_plot %>%
      filter(censoring_method == censK) %>%
      filter(outcome_name == outcome) %>%
      filter(adjustment == "calibrated") %>%
      filter(comparison %in% c("any (VE)", "pfizer (VE)", "astrazeneca (VE)")) %>%
      mutate(comparison = case_when(
        comparison == "any (VE)" ~ "Any vaccine",
        comparison == "pfizer (VE)" ~ "BNT162b2",
        comparison == "astrazeneca (VE)" ~ "ChAdOx1"
      )) %>%
      mutate(
        comparison = factor(comparison, levels = c("Any vaccine", "BNT162b2", "ChAdOx1"))
      ) %>%
      mutate(
        cdm_name = factor(cdm_name, levels = c("AURUM", "GOLD", "SIDIAP", "CORIVA"))
      ) %>%
      mutate(adjustment_y = case_when(
        study == "Study 1" ~ 5,
        study == "Study 2" ~ 4,
        study == "Study 3" ~ 3,
        study == "Study 4" ~ 2,
        study == "Meta Analysis" ~ 0.5
      )) %>%
      mutate(meta = if_else(study == "Meta Analysis", "meta", "no meta")) 
    
    if (outcome == "next_post_acute_covid19" | 
        outcome == "longcovid_post_acute_covid19_90_365" | 
        outcome == "longcovid_post_acute_covid19_28_365") {
      data <- data %>%
        filter(cdm_name != "CORIVA") %>%
        filter(cdm_name != "UiO") %>%
        filter(cdm_name != "GOLD")
    }
    
    data %>%
      ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
      geom_point(aes(x=hr), size=1) +
      geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
      facet_grid(cdm_name ~ comparison) +
      # labs(title = outcome) +
      scale_y_continuous(limits = c(0, 5.5), breaks = c(0.5, 2, 3, 4, 5), labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
      scale_x_log10(limits = c(0.1, 4), breaks = c(0.1, 0.25, 0.50, 1, 2, 4), oob=scales::rescale_none) +
      geom_vline(xintercept = 1) +
      geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
      theme(
        axis.title.y = element_blank(),
        legend.title = element_blank()
      ) 
    
    ggsave(
      paste0("Figure2_", outcome, "_", censK,".png"),
      plot = last_plot(),
      path = here("figures"),
      scale = 1,
      width = 3000,
      height = 2000,
      units = "px",
      dpi = 300,
      device = "png"
    )
  }
}
for (censK in cens) {
  for (outcome in outcomes) {
    data_1 <- estimates_plot %>%
      filter(censoring_method == censK) %>%
      filter(outcome_name == outcome) %>%
      filter(adjustment == "calibrated") %>%
      filter(comparison == "pfizer - astrazeneca (CVE)") %>%
      mutate(comparison = "BNT162b2 - ChAdOx1") %>%
      filter(cdm_name %in% c("AURUM", "GOLD")) %>%
      mutate(
        cdm_name = factor(cdm_name, levels = c("AURUM", "GOLD"))
      ) %>%
      mutate(adjustment_y = case_when(
        study == "Study 1" ~ 5,
        study == "Study 2" ~ 4,
        study == "Study 3" ~ 3,
        study == "Study 4" ~ 2,
        study == "Meta Analysis" ~ 0.5
      )) %>%
      mutate(meta = if_else(study == "Meta Analysis", "meta", "no meta")) %>%
      filter(!is.na(se_coef))
    
    data_2 <- data_1 %>%
      filter(cdm_name == "null")
      
      # estimates_plot %>%
      # filter(censoring_method == censK) %>%
      # filter(outcome_name == outcome) %>%
      # filter(adjustment == "calibrated") %>%
      # filter(comparison == "pfizer - moderna (CVE)") %>%
      # mutate(comparison = "BNT162b2 - mRNA-1273") %>%
      # filter(cdm_name %in% c("UiO", "CORIVA")) %>%
      # mutate(
      #   cdm_name = factor(cdm_name, levels = c("CORIVA", "UiO"))
      # ) %>%
      # mutate(adjustment_y = case_when(
      #   study == "Study 1" ~ 5,
      #   study == "Study 2" ~ 4,
      #   study == "Study 3" ~ 3,
      #   study == "Study 4" ~ 2,
      #   study == "Meta Analysis" ~ 0.5
      # )) %>%
      # mutate(meta = if_else(study == "Meta Analysis", "meta", "no meta")) %>%
      # filter(!is.na(se_coef))
    
    if (outcome == "next_post_acute_covid19" | 
        outcome == "longcovid_post_acute_covid19_90_365" | 
        outcome == "longcovid_post_acute_covid19_28_365") {
      data_1 <- data_1 %>%
        filter(cdm_name != "CORIVA") %>%
        filter(cdm_name != "UiO") %>%
        filter(cdm_name != "GOLD")
      
      data_2 <- data_2 %>%
        filter(cdm_name == "lala") %>%
        filter(cdm_name != "UiO") %>%
        filter(cdm_name != "GOLD")
    }
    
    if (nrow(data_1) > 0 & nrow(data_2) > 0) {
      gg_1 <- data_1 %>%
        ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
        geom_point(aes(x=hr), size=1) +
        geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
        facet_grid(cdm_name ~ comparison) +
        
        scale_y_continuous(limits = c(0, 5.5), breaks = c(0.5, 2, 3, 4, 5), 
                           labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
        # scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5), oob=scales::rescale_none) +
        scale_x_log10(limits = c(0.1, 2), breaks = c(0.1, 0.25, 0.50, 1, 2), oob=scales::rescale_none) +
        geom_vline(xintercept = 1) +
        geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
        theme(
          axis.title.y = element_blank(),
          legend.position = "none"
        )
      
      gg_2 <- data_2 %>%
        ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
        geom_point(aes(x=hr), size=1) +
        geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
        facet_grid(cdm_name ~ comparison) +
        
        scale_y_continuous(limits = c(0, 5.5), breaks = c(0.5, 2, 3, 4, 5), 
                           labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
        scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5), oob=scales::rescale_none) +
        geom_vline(xintercept = 1) +
        geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
        theme(
          axis.title.y = element_blank(),
          legend.position = "none"
        )
      
      gg <- ggarrange(gg_1, gg_2, ncol = 2)
      
      ggsave(
        paste0("Figure3_", outcome, "_", censK,".png"),
        plot = gg,
        path = here("figures"),
        scale = 1,
        width = 2000,
        height = 1000,
        units = "px",
        dpi = 300
      )
      
    
    } else if (nrow(data_1) == 0) {
      gg <- data_2 %>%
        ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
        geom_point(aes(x=hr), size=1) +
        geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
        facet_grid(cdm_name ~ comparison) +
        
        scale_y_continuous(limits = c(0, 5.5), breaks = c(0.5, 2, 3, 4, 5), 
                           labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
        scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5), oob=scales::rescale_none) +
        geom_vline(xintercept = 1) +
        geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
        theme(
          axis.title.y = element_blank(),
          legend.position = "none"
        )
      
      ggsave(
        paste0("Figure3_", outcome, "_", censK,".png"),
        plot = gg,
        path = here("figures"),
        scale = 1,
        width = 1200,
        height = 750,
        units = "px",
        dpi = 300
      )
      
    } else if (nrow(data_2) == 0) {
      gg <- data_1 %>%
        ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
        geom_point(aes(x=hr), size=1) +
        geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
        facet_grid(cdm_name ~ comparison) +
        
        scale_y_continuous(limits = c(0, 5.5), breaks = c(0.5, 2, 3, 4, 5), 
                           labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
        # scale_x_continuous(limits = c(0, 3), breaks = seq(0, 3, by = 0.5), oob=scales::rescale_none) +
        scale_x_log10(limits = c(0.1, 4), breaks = c(0.1, 0.25, 0.50, 1, 2,4), oob=scales::rescale_none) +
        geom_vline(xintercept = 1) +
        geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
        theme(
          axis.title.y = element_blank(),
          legend.position = "none"
        )
      
      ggsave(
        paste0("Figure3_", outcome, "_", censK,".png"),
        plot = gg,
        path = here("figures"),
        scale = 1,
        width = 1200,
        height = 750,
        units = "px",
        dpi = 300
      )
      }
    
    
   
    
    # rm(data)
  }
}
# Figure 3 comparative vaccie effectiveness ----
# Figure 4 kaplan meier plots ----
# source(here("prep_data_lc.R"))
# for(ii in 1:nrow(captions_km)) {
#   data <- survival_plot %>%
#     filter(censoring_method == captions_km$censoring_method[ii]) %>%
#     filter(outcome_name == captions_km$outcome_name[ii]) %>%
#     filter(exposure_name == captions_km$exposure_name[ii]) %>%
#     filter(comparator_name == captions_km$comparator_name[ii])
# 
#   if (captions_km$cve[ii]) {
#     data <- data %>%
#       mutate(cohort = ifelse(grepl("unvaccinated", .data$group), "ChAdOx1", "BNT162b2"))
#   } else {
#     data <- data %>%
#       mutate(cohort = ifelse(grepl("unvaccinated", .data$group), "Unvaccinated", "Vaccinated"))
#   }
# 
#   if(captions_outcome$exposure_name[ii] == "astrazeneca vaccinated" | captions_km$cve[ii]) {
#   data <- data %>%
#     filter(cdm_name != "CORIVA")
#   }
# 
#   if(captions_km$outcome_name[ii] %in% c("longcovid_post_acute_covid19_90_365", "next_post_acute_covid19")) {
#     data <- data %>%
#     filter(cdm_name != "CORIVA") %>%
#       filter(cdm_name != "GOLD")
#     if (captions_km$cve[ii]) {
#       if (nrow(data) > 1) {
#         data <- data[1,]
#       }
#     }
#   }
# 
#   if (nrow(data) > 2) {
#   fit <- survfit(Surv(time, status) ~ cohort + cdm_name + study, data = data, weights = weight)
# 
#   surv_data <- surv_summary(fit, data) %>%
#     mutate(cdm_name = factor(cdm_name, level = c("AURUM","GOLD", "SIDIAP", "CORIVA")))
# 
#   ylim_lower <- surv_data %>%
#     group_by(cdm_name) %>%
#     select(surv) %>%
#     filter(surv == min(surv)) %>%
#     distinct()
# 
#   y <- list()
#   for (num in 1:length(unique(surv_data$cdm_name))) {
#     y[[num]] <- eval(parse(text =
#                            paste0('cdm_name =="', ylim_lower$cdm_name[num], '"~ scale_y_continuous(limits = c(', ylim_lower$surv[num], ', 1), breaks = c(', paste(matlab::linspace(ylim_lower$surv[num],1,4), collapse= ", "), '), labels = c(', paste(sprintf(matlab::linspace(ylim_lower$surv[num],1,4), fmt = '%#.4f'), collapse = ","), '))')))
# 
#   }
#   ggsurv <- ggsurvplot(fit,
#                        data = data,
#                        censor = FALSE,
#                        color = "cohort",
#                        palette = c("#8fb996", "#457b9d"))
# 
# 
#   # cat(gsub(".  ",".",paste0("Figure ", num_fig, ": ", captions_km$caption[ii])))
# 
#   ggsurv$plot +
#           theme_gray() +
#           theme(legend.position = "right",
#                 legend.title = element_blank())+
#           facet_grid(cdm_name ~ study, scales = "free_y") +
#           facetted_pos_scales(y = y)
#   
# 
#   name <- paste0("KM_",
#                  captions_km$outcome_name[ii], 
#                  "_", 
#                  captions_km$censoring_method[ii],
#                  "_",
#                  gsub(" ", "_", captions_km$exposure_name[ii]))
#   
#   if (captions_km$cve[ii]) {
#     name <- paste0(name, "_cve")
#   }
#   ggsave(
#     paste0(name, ".png"),
#     plot = last_plot(),
#     path = here("figures"),
#     scale = 1,
#     width = 3000,
#     height = 1500,
#     units = "px",
#     dpi = 300
#   )
#   }
# }
# Figure 5 cumulative incidence ----
# Figure 6 negative control outcomes ----
