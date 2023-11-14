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

# rename aurum and gold
renameAG <- function(x) {
  x %>% 
    mutate(cdm_name = if_else(cdm_name == "AURUM", "CPRD AURUM", cdm_name)) %>%
    mutate(cdm_name = if_else(cdm_name == "GOLD", "CPRD GOLD", cdm_name)) 
}
asmd <- renameAG(asmd)
cohort_details <- renameAG(cohort_details)
comparisons <- renameAG(comparisons)
nco <- renameAG(nco)
estimates <- renameAG(estimates)
study_attrition <- renameAG(study_attrition)
table_characteristics <- renameAG(table_characteristics)
table_characteristics_weighted <- renameAG(table_characteristics_weighted)
survival_plot <- renameAG(survival_plot)

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
  filter(!(cdm_name == "SIDIAP" & grepl("astrazeneca", comparison_name) & study %in% c(1, 4))) %>%
  filter(!(cdm_name == "CORIVA" & grepl("astrazeneca", comparison_name) & study %in% c(2, 3, 4))) %>%
  filter(skip == 0)
study_attrition <- study_attrition %>%
  mutate(excluded = if_else(is.na(excluded), 0, excluded)) %>%
  mutate(number_observations = nice(number_observations)) %>%
  mutate(excluded = nice(excluded))
nco <- nco %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name"))
estimates <- estimates %>%
  inner_join(comparisons, by = c("comparison_id", "cdm_name")) %>%
  filter(model == "finegray")
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
      filter(cdm_name == "CPRD GOLD") %>%
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
p <- ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
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
  "Figure1_GOLD.png", device = "png",
  plot = p,
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
      filter(cdm_name == "CPRD AURUM") %>%
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
p <-ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
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
  "Figure1_AURUM.png", device = "png",
  plot = p,
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
p <- ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
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
  "Figure1_SIDIAP.png", device = "png",
  plot = p,
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
p <- ggplot(data = data_plot, aes(x = asmd_unadjusted, y = asmd_adjusted)) +
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
  "Figure1_CORIVA.png", device = "png",
  plot = p,
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
    p <- estimates_plot %>%
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
      ) 
    if (!grepl("post_acute_covid19", outcome)) {
      p <- p %>%
        mutate(
          cdm_name = factor(cdm_name, levels = c("CPRD AURUM", "CPRD GOLD", "SIDIAP", "CORIVA"))
        )
    } else {
      p <- p %>%
        filter(cdm_name %in% c("CPRD AURUM", "SIDIAP")) %>%
        mutate(
          cdm_name = factor(cdm_name, levels = c("CPRD AURUM", "SIDIAP"))
        )
    }
    p <- p %>%
      mutate(adjustment_y = case_when(
        study == "Study 1" ~ 5,
        study == "Study 2" ~ 4,
        study == "Study 3" ~ 3,
        study == "Study 4" ~ 2,
        study == "Meta Analysis" ~ 0.5
      )) %>%
      mutate(meta = if_else(study == "Meta Analysis", "meta", "no meta")) %>%
      ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
      geom_point(aes(x=hr, shape = adjustment), size=1) +
      geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
      facet_grid(cdm_name ~ comparison) +
      labs(x = "Hazard ratio") +
      scale_y_continuous(breaks = c(0.5, 2, 3, 4, 5), labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
      scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 1, 2), trans = "log10") + 
      coord_cartesian(ylim = c(0, 5.5), xlim = c(0.1, 2)) +
      geom_vline(xintercept = 1) +
      geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
      theme(
        axis.title.y = element_blank(), legend.position = "none"
      )
    ggsave(
      paste0("Figure2_", outcome, "_", censK,".png"), device = "png",
      plot = p,
      path = here("figures"),
      scale = 1,
      width = 3000,
      height = 2000,
      units = "px",
      dpi = 300
    )
  }
}
for (censK in cens) {
  for (outcome in outcomes) {
    p <- estimates_plot %>%
      filter(censoring_method == censK) %>%
      filter(outcome_name == outcome) %>%
      filter(adjustment == "calibrated") %>%
      filter(comparison == "pfizer - astrazeneca (CVE)") %>%
      mutate(comparison = "BNT162b2 - ChAdOx1") %>%
      filter(cdm_name %in% c("CPRD AURUM", "CPRD GOLD")) %>%
      mutate(
        cdm_name = factor(cdm_name, levels = c("CPRD AURUM", "CPRD GOLD"))
      ) %>%
      mutate(adjustment_y = case_when(
        study == "Study 1" ~ 5,
        study == "Study 2" ~ 4,
        study == "Study 3" ~ 3,
        study == "Study 4" ~ 2,
        study == "Meta Analysis" ~ 0.5
      )) %>%
      mutate(meta = if_else(study == "Meta Analysis", "meta", "no meta")) %>%
      ggplot(aes(y = adjustment_y, colour = study, fill = study)) +
      geom_point(aes(x=hr, shape = adjustment), size=1) +
      geom_linerange(aes(xmin=lower_hr, xmax=upper_hr)) +
      facet_grid(cdm_name ~ comparison) +
      labs(x = "Hazard ratio") +
      scale_y_continuous(breaks = c(0.5, 2, 3, 4, 5), labels = c("meta-analysis", "cohort 4", "cohort 3", "cohort 2", "cohort 1")) +
      scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 1, 2), trans = "log10") + 
      coord_cartesian(ylim = c(0, 5.5), xlim = c(0.1, 3)) +
      geom_vline(xintercept = 1) +
      geom_hline(yintercept = 1.25, linetype="dashed", color = "black") +
      theme(
        axis.title.y = element_blank(), legend.position = "none"
      )
    ggsave(
      paste0("Figure3_", outcome, "_", censK,".png"), device = "png",
      plot = p,
      path = here("figures"),
      scale = 1,
      width = 1500,
      height = 1000,
      units = "px",
      dpi = 300
    )
  }
}

