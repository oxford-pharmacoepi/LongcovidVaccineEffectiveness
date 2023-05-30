addAge <- function(x, cdm, ageAt = "cohort_start_date", defaultMonth = 1, defaultDay = 1, imposeMonth = TRUE, imposeDay = TRUE, compute = TRUE) {
  defaultMonth <- as.integer(defaultMonth)
  defaultDay <- as.integer(defaultDay)
  
  person <- cdm[["person"]] %>%
    dplyr::rename("subject_id" = "person_id") %>%
    dplyr::inner_join(
      x %>%
        dplyr::select("subject_id", dplyr::all_of(ageAt)) %>%
        dplyr::distinct(),
      by = "subject_id"
    )
  
  if (imposeMonth == TRUE) {
    person <- person %>%
      dplyr::mutate(month_of_birth = .env$defaultMonth)
  } else {
    person <- person %>%
      dplyr::mutate(month_of_birth = dplyr::if_else(
        is.na(.data$month_of_birth),
        .env$defaultMonth,
        .data$month_of_birth
      ))
  }
  
  if (imposeDay == TRUE) {
    person <- person %>%
      dplyr::mutate(day_of_birth = .env$defaultDay)
  } else {
    person <- person %>%
      dplyr::mutate(day_of_birth = dplyr::if_else(
        is.na(.data$day_of_birth),
        .env$defaultDay,
        .data$day_of_birth
      ))
  }
  
  person <- person %>%
    dplyr::filter(!is.na(.data$year_of_birth)) %>%
    dplyr::mutate(year_of_birth1 = as.character(as.integer(.data$year_of_birth))) %>%
    dplyr::mutate(month_of_birth1 = as.character(as.integer(.data$month_of_birth))) %>%
    dplyr::mutate(day_of_birth1 = as.character(as.integer(.data$day_of_birth))) %>%
    dplyr::mutate(birth_date = as.Date(paste0(
      .data$year_of_birth1, "-",
      .data$month_of_birth1, "-",
      .data$day_of_birth1
    ))) %>%
    dplyr::mutate(age = floor(dbplyr::sql(sqlGetAge(
      dialect = CDMConnector::dbms(cdm),
      dob = "birth_date",
      dateOfInterest = ageAt
    )))) %>%
    dplyr::select("subject_id", dplyr::all_of(ageAt), "age") %>%
    dplyr::right_join(x, by = c("subject_id", ageAt)) %>%
    dplyr::select(dplyr::all_of(colnames(x)), "age")
  if (isTRUE(compute)) {
    person <- person %>% dplyr::compute()
  }
  return(person)
}

sqlGetAge <- function(dialect, dob, dateOfInterest) {
  SqlRender::translate(
    SqlRender::render("((YEAR(@date_of_interest) * 10000 + MONTH(@date_of_interest) * 100 +
                      DAY(@date_of_interest)-(YEAR(@dob)* 10000 + MONTH(@dob) * 100 + DAY(@dob))) / 10000)",
                      dob = dob,
                      date_of_interest = dateOfInterest
    ),
    targetDialect = dialect
  )
}

addAgeGroup <- function(x, ageGroup = NULL, name = "age_group", compute = TRUE) {
  ageGroup <- lapply(ageGroup, function(xx) {
    xx[1] <- ifelse(is.na(xx[1]), 0, xx[1])
    xx[2] <- ifelse(is.na(xx[2]), 150, xx[2])
    nam <- paste0(xx[1], ";", xx[2])
    xx <- dplyr::tibble(age_min = xx[1], age_max = xx[2], !!name := nam)
    return(xx)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate(to_join = 1)
  ageGroup <- ageGroup %>%
    dplyr::inner_join(
      dplyr::tibble(
        to_join = 1,
        age = seq(min(ageGroup$age_min), max(ageGroup$age_max))
      ),
      multiple = "all",
      by = "to_join"
    ) %>%
    dplyr::filter(.data$age >= .data$age_min & .data$age <= .data$age_max) %>%
    dplyr::select(dplyr::all_of(c("age", name)))
  x <- x %>%
    dplyr::left_join(ageGroup, by = "age", copy = TRUE)
  if (isTRUE(compute)) {
    x <- x %>% dplyr::compute()
  }
  return(x)
}

addPriorHistory <- function(x, cdm, priorHistoryAt = "cohort_start_date", compute = TRUE) {
  x <- cdm[["observation_period"]] %>%
    dplyr::select(
      "subject_id" = "person_id", "observation_period_start_date"
    ) %>%
    dplyr::inner_join(
      x %>%
        dplyr::select("subject_id", dplyr::all_of(priorHistoryAt)) %>%
        dplyr::distinct(),
      by = "subject_id"
    ) %>%
    dplyr::mutate(prior_history = CDMConnector::datediff(
      start = "observation_period_start_date",
      end = !!priorHistoryAt
    )) %>%
    dplyr::right_join(
      x,
      by = c("subject_id", priorHistoryAt)
    ) %>%
    dplyr::select(dplyr::all_of(colnames(x)), "prior_history")
  if (isTRUE(compute)) {
    x <- x %>% dplyr::compute()
  }
  return(x)
}

addSex <- function(x, cdm, compute = TRUE) {
  x <- cdm[["person"]] %>%
    dplyr::rename("subject_id" = "person_id") %>%
    dplyr::inner_join(
      x %>% dplyr::select("subject_id") %>% dplyr::distinct(),
      by = c("subject_id")
    ) %>%
    dplyr::mutate(sex = dplyr::case_when(
      .data$gender_concept_id == 8507 ~ "Male",
      .data$gender_concept_id == 8532 ~ "Female",
      TRUE ~ as.character(NA)
    )) %>%
    dplyr::select("subject_id", "sex") %>%
    dplyr::right_join(x, by = "subject_id") %>%
    dplyr::select(dplyr::all_of(colnames(x)), "sex")
  if (isTRUE(compute)) {
    x <- x %>% dplyr::compute()
  }
  return(x)
}

addNumberVisit <- function(x, cdm, window = c(NA, NA), name = "number_visit", compute = TRUE) {
  result <- cdm[["visit_occurrence"]] %>%
    dplyr::select(
      "subject_id" = "person_id", "visit_concept_id", "visit_start_date"
    ) %>%
    dplyr::inner_join(
      x %>%
        dplyr::select("subject_id", "cohort_start_date") %>%
        dplyr::distinct(),
      by = "subject_id"
    )
  if (!is.na(window[1])) {
    result <- result %>%
      dplyr::filter(
        CDMConnector::dateadd("cohort_start_date", window[1]) <=
          .data$visit_start_date
      )
  }
  if (!is.na(window[2])) {
    result <- result %>%
      dplyr::filter(
        CDMConnector::dateadd("cohort_start_date", window[2]) >=
          .data$visit_start_date
      )
  }
  result <- result %>%
    dplyr::group_by(
      .data$subject_id, .data$cohort_start_date
    ) %>%
    dplyr::summarise(!!name := dplyr::n(), .groups = "drop") %>%
    dplyr::compute() %>%
    dplyr::right_join(
      x,
      by = c("subject_id", "cohort_start_date")
    ) %>%
    dplyr::mutate(!!name := dplyr::if_else(is.na(.data[[name]]), 0, .data[[name]])) %>%
    dplyr::select(dplyr::all_of(c(colnames(x), name)))
  if (isTRUE(compute)) {
    result <- result %>% dplyr::compute()
  }
  return(result)
}

addNumberEvent <- function(x, cdm, eventTableName, filter = NULL, window = c(NA, NA), name = "number_event", eventAt = "cohort_start_date", eventDate = "cohort_start_date", compute = TRUE) {
  events <- cdm[[eventTableName]]
  if (!is.null(filter)) {
    namesFilter <- names(filter)
    for (k in 1:length(filter)) {
      events <- events %>%
        dplyr::filter_at(
          dplyr::vars(dplyr::all_of(namesFilter[k])), 
          dplyr::any_vars(. == !!filter[[k]])
        )
    }
  }
  if ("person_id" %in% colnames(events)) {
    events <- events %>%
      dplyr::rename("subject_id" = "person_id")
  }
  events <- events %>%
    dplyr::select("subject_id","event_date" = dplyr::all_of(eventDate)) %>%
    dplyr::inner_join(
      x %>%
        dplyr::select("subject_id", "event_at" = dplyr::all_of(eventAt)) %>%
        dplyr::distinct(),
      by = "subject_id"
    ) %>%
    dplyr::mutate(date_dif = !!CDMConnector::datediff("event_at", "event_date"))
  if (!is.na(window[1])) {
    events <- events %>%
      dplyr::filter(.data$date_dif >= !!window[1])
  }
  if (!is.na(window[2])) {
    events <- events %>%
      dplyr::filter(.data$date_dif <= !!window[2])
  }
  events <- events %>%
    dplyr::group_by(.data$subject_id, .data$event_at) %>%
    dplyr::summarise(!!name := dplyr::n(), .groups = "drop") %>%
    dplyr::rename(!!eventAt := "event_at") %>%
    dplyr::right_join(
      x,
      by = c("subject_id", eventAt)
    ) %>%
    dplyr::mutate(dplyr::across(
      dplyr::all_of(name), 
      ~ dplyr::if_else(is.na(.x), 0, .x)
    )) %>%
    dplyr::select(dplyr::all_of(c(colnames(x), name))) 
  if (isTRUE(compute)) {
    events <- events %>% dplyr::compute()
  }
  return(events)
}

addEvent <- function(x, cdm, eventTableName, eventId = NULL, window = c(NA, NA), name = "event", eventDate = "cohort_start_date", eventAt = "cohort_start_date", order = "first", compute = TRUE) {
  eventTable <- cdm[[eventTableName]]
  if (!is.null(eventId)) {
    eventTable <- eventTable %>% dplyr::filter(.data$cohort_definition_id %in% .env$eventId)
  }
  if ("person_id" %in% colnames(eventTable)) {
    eventTable <- eventTable %>% dplyr::rename("subject_id" = "person_id")
  }
  xx <- x %>%
    dplyr::select(dplyr::all_of(c("subject_id", eventAt))) %>%
    dplyr::distinct() %>%
    dplyr::inner_join(
      eventTable %>%
        dplyr::select("subject_id", "event_date" = dplyr::all_of(eventDate)),
      by = "subject_id"
    ) %>%
    dplyr::mutate(dif_time = !!CDMConnector::datediff(eventAt, "event_date"))
  if (!is.na(window[1])) {
    xx <- xx %>%
      dplyr::filter(.data$dif_time >= !!window[1])
  }
  if (!is.na(window[2])) {
    xx <- xx %>%
      dplyr::filter(.data$dif_time <= !!window[2])
  }
  xx <- xx %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("subject_id", eventAt))))
  if (order == "first") {
    xx <- xx %>%
      dplyr::summarise(
        !!name := min(.data$event_date, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (order == "last") {
    xx <- xx %>%
      dplyr::summarise(
        !!name := max(.data$event_date, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (order == "all") {
    xx <- xx %>%
      dbplyr::window_order(.data$event_date) %>%
      dplyr::mutate(nam = dplyr::row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(nam = paste0(.env$name, "_", .data$nam)) %>%
      tidyr::pivot_wider(names_from = "nam", values_from = "event_date")
  }
  xx <- x %>% dplyr::left_join(xx, by = c("subject_id", eventAt))
  if (isTRUE(compute)) {
    xx <- xx %>% dplyr::compute()
  }
  return(xx)
}

addRegionalData <- function(x, cdm) {
  database <- cdmName(cdm)
  result <- x %>%
    dplyr::select("subject_id", "cohort_start_date") %>%
    dplyr::distinct()
  if (database %in% c("AURUM", "GOLD")) {
    regionalData <- read_csv(here("Data", paste0("data_", database, ".csv")), show_col_types = FALSE) %>%
      mutate(dose1_fraction = if_else(is.na(dose1_fraction), 0, dose1_fraction)) %>%
      mutate(dose2_fraction = if_else(is.na(dose2_fraction), 0, dose2_fraction)) %>%
      rename("cohort_start_date" = "date")
    result <- result %>%
      inner_join(
        cdm$person %>% select("subject_id" = "person_id", "care_site_id"),
        by = "subject_id"
      ) %>%
      inner_join(cdm$care_site, by = "care_site_id") %>%
      inner_join(cdm$location, by = "location_id") %>%
      select("subject_id", "cohort_start_date", "region" = "location_source_value") %>%
      mutate(region = tolower(region)) %>%
      mutate(region = if_else(
        region %in% c("south east coast", "south central"),
        "south east",
        region
      ))
    if (database == "GOLD") {
      result <- result %>%
        mutate(region = if_else(
          region == "yorkshire  & the humber",
          "yorkshire and the humber",
          region
        ))
    }
  } else if (database == "SIDIAP") {
    regionalData <- read_csv(here("Data", paste0("data_", database, ".csv")), show_col_types = FALSE) %>%
      mutate(dose1_fraction = if_else(is.na(dose1_fraction), 0, dose1_fraction)) %>%
      mutate(dose2_fraction = if_else(is.na(dose2_fraction), 0, dose2_fraction)) %>%
      rename("cohort_start_date" = "date") %>%
      mutate(region = tolower(region))
    result <- result %>%
      inner_join(
        cdm$person %>% 
          select("subject_id" = "person_id", "location_id") %>%
          inner_join(
            cdm$location %>%
              select("location_id", "region" = "location_source_value"), 
            by = "location_id"
          ) %>%
          select(-"location_id"),
        by = "subject_id"
      ) %>%
      mutate(region = tolower(region)) 
  }
  
  if (database %in% c("UiO", "CORIVA")) {
    result <- result %>%
      mutate(
        region = "all country",
        new_cases_7d_100k = 0,
        new_tests_7d_100k = 0,
        dose1_fraction = 0,
        dose2_fraction = 0
      )
  } else {
    result <- result %>%
      left_join(regionalData, by = c("region", "cohort_start_date"), copy = TRUE) 
  }
  
  result <- result %>%
    right_join(x, by = c("subject_id", "cohort_start_date")) %>%
    compute()
  
  if (result %>%
      filter(
        is.na(new_cases_7d_100k) | is.na(new_tests_7d_100k) |
        is.na(dose1_fraction) | is.na(dose2_fraction)
      ) %>%
      tally() %>%
      pull() > 0) {
    stop("missing regional data for some subjects")
  }
  return(result)
}

addGPIdentifier <- function(x, cdm) {
  database <- cdmName(cdm)
  result <- x %>%
    dplyr::select("subject_id") %>%
    dplyr::distinct()
  if (database %in% c("GOLD", "AURUM")) {
    result <- result %>%
      inner_join(
        cdm$person %>% select("subject_id" = "person_id", "gp" = "care_site_id"),
        by = "subject_id"
      )
  } else if (database %in% c("UiO", "CORIVA")) {
    result <- result %>%
      mutate(gp = "unique gp")
  } else if (database == "SIDIAP") {
    result <- result %>%
      inner_join(
        cdm$person %>% 
          select("subject_id" = "person_id", "location_id") %>%
          inner_join(
            cdm$location %>%
              select("location_id", "gp" = "location_source_value"), 
            by = "location_id"
          ) %>%
          select(-"location_id"),
        by = "subject_id"
      ) %>%
      mutate(gp = tolower(gp))
  }
  result <- result %>%
    right_join(x, by = "subject_id") %>%
    compute()
  if (result %>% filter(is.na(gp)) %>% tally() %>% pull() > 0) {
    stop("missing gp identifier for some subjects")
  }
  return(result)
}

asmdBinary <- function(x, variables = NULL, groupName = "group", weight = NULL) {
  if (is.null(variables)) {
    variables <- colnames(x)[!(colnames(x) %in% c(groupName, weight))]
  }
  if (is.null(weight)) {
    x <- x %>% mutate(weight = 1)
  } else {
    x <- x %>% rename("weight" = dplyr::all_of(weight))
  }
  x <- x %>% rename("group" = dplyr::all_of(groupName))
  for (variable in variables) {
    lab <- unique(x[[variable]])
    if (!all(lab %in% c(0, 1))) {
      x <- dplyr::mutate(x, !!variable := dplyr::if_else(.data[[variable]] == .env$lab[1], 0, 1))
    }
  }
  lab <- unique(x$group)
  x %>%
    dplyr::mutate(group = dplyr::if_else(.data$group == .env$lab[1], 1, 2)) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(dplyr::across(
      dplyr::all_of(variables),
      list(
        mean = function(x) {Hmisc::wtd.mean(x, .data$weight)},
        var = function(x) {Hmisc::wtd.var(x, .data$weight)}
      ),
      .names = "{.col} {.fn}"
    )) %>%
    tidyr::pivot_longer(!"group", names_to = "variable.func") %>%
    tidyr::separate("variable.func", c("variable", "func"), " ") %>%
    tidyr::pivot_wider(names_from = c("func", "group"), values_from = "value") %>%
    dplyr::mutate(asmd = abs(.data$mean_1-.data$mean_2)/sqrt((.data$var_1+.data$var_2)/2)) %>%
    dplyr::select("variable", "asmd") %>%
    mutate(asmd_type = "binary")
  
}

asmdContinuous <- function(x, variables = NULL, groupName = "group", weight = NULL) {
  if (is.null(variables)) {
    variables <- colnames(x)[!(colnames(x) %in% c(groupName, weight))]
  }
  if (is.null(weight)) {
    x <- x %>% mutate(weight = 1)
  } else {
    x <- x %>% rename("weight" = dplyr::all_of(weight))
  }
  x <- x %>% rename("group" = dplyr::all_of(groupName))
  lab <- unique(x$group)
  x %>%
    dplyr::mutate(group = dplyr::if_else(.data$group == .env$lab[1], 1, 2)) %>%
    dplyr::group_by(.data$group) %>%
    dplyr::summarise(dplyr::across(
      dplyr::all_of(variables),
      list(
        mean = function(x) {Hmisc::wtd.mean(x, .data$weight)},
        var = function(x) {Hmisc::wtd.var(x, .data$weight)}
      ),
      .names = "{.col} {.fn}"
    )) %>%
    tidyr::pivot_longer(!"group", names_to = "variable.func") %>%
    tidyr::separate("variable.func", c("variable", "func"), " ") %>%
    tidyr::pivot_wider(names_from = c("func", "group"), values_from = "value") %>%
    dplyr::mutate(asmd = abs(.data$mean_1-.data$mean_2)/sqrt(.data$var_1+.data$var_2)) %>%
    dplyr::select("variable", "asmd") %>%
    mutate(asmd_type = "continuous")
  
}

asmdCategorical <- function(x, variables = NULL, groupName = "group", weight = NULL) {
  if (is.null(variables)) {
    variables <- colnames(x)[!(colnames(x) %in% c(groupName, weight))]
  }
  if (is.null(weight)) {
    x <- x %>% dplyr::mutate(weight = 1)
  } else {
    x <- x %>% dplyr::rename("weight" = dplyr::all_of(weight))
  }
  x <- x %>% 
    dplyr::rename("group" = dplyr::all_of(groupName)) %>%
    dplyr::select("group", "weight", dplyr::all_of(variables))
  lab <- unique(x$group)
  if (length(lab) != 2) {
    stop("Number of labels in group column different from 2.")
  }
  x <- x %>% 
    dplyr::mutate(group = dplyr::if_else(.data$group == .env$lab[1], 1, 0))
  denominator <- x %>% 
    dplyr::group_by(.data$group) %>%
    dplyr::tally(wt = .data$weight, name = "denominator")
  result <- NULL
  for (k in 1:length(variables)) {
    y <- x %>% 
      dplyr::rename("label" = dplyr::all_of(variables[k])) %>%
      dplyr::group_by(.data$group,.data$label) %>%
      dplyr::tally(wt = .data$weight) %>%
      dplyr::right_join(denominator, by = "group") %>%
      dplyr::mutate(percentage = dplyr::if_else(
        is.na(.data$n), 0, .data$n/.data$denominator
      )) %>%
      dplyr::select("label", "group", "percentage") %>%
      tidyr::pivot_wider(names_from = "group", values_from = "percentage", values_fill = 0)
    TT <- y[["1"]]
    CC <- y[["0"]]
    result <- result %>% dplyr::union_all(dplyr::tibble(
      variable = variables[k], 
      asmd = asmdFromPercentage(TT, CC)
    ))
  }
  result <- result %>% mutate(asmd_type = "categorical")
  return(result)
}

asmdFromPercentage <- function(TT, CC) {
  if (length(TT) == 1) {
    return(NA)
  } else {
    TT <- TT[-1]
    CC <- CC[-1]
    n <- length(TT)
    vect <- TT - CC
    TT1 <- matrix(rep(TT, n), nrow = n, byrow = TRUE)
    TT2 <- matrix(rep(TT, n), nrow = n, byrow = FALSE)
    CC1 <- matrix(rep(CC, n), nrow = n, byrow = TRUE)
    CC2 <- matrix(rep(CC, n), nrow = n, byrow = FALSE)
    S <- (TT1*TT2 + CC1*CC2) / 2
    diag(S) <- (TT*(1-TT) + CC*(1-CC)) / 2
    asmd <- as.numeric(sqrt(vect %*% solve(S) %*% vect))
    return(asmd)
  }
}

getCohortCount <- function(cohort) {
  cohort %>%
    group_by(cohort_definition_id) %>%
    summarise(
      number_records = dplyr::n(),
      number_subjects = dplyr::n_distinct(.data$subject_id),
      .groups = "drop"
    ) %>%
    collect() %>%
    arrange(cohort_definition_id)
}

addOverlap <- function(x, cdm, cohortName, id, name, window = c(NA, 0)) {
  xx <- x %>%
    select("subject_id", "cohort_start_date") %>%
    inner_join(
      cdm[[cohortName]] %>%
        filter(cohort_definition_id %in% !!id) %>%
        inner_join(
          dplyr::tibble(
            cohort_definition_id = id,
            overlap_name = name
          ),
          by = "cohort_definition_id",
          copy = TRUE
        ) %>%
        select(
          "subject_id", 
          "overlap_name", 
          "overlap_start" = "cohort_start_date", 
          "overlap_end" = "cohort_end_date"
        ),
      by = "subject_id"
    )
  if (!is.na(window[1])) {
    xx <- xx %>%
      filter(overlap_end >= CDMConnector::dateadd("cohort_start_date", window[1]))
  }
  if (!is.na(window[2])) {
    xx <- xx %>%
      filter(overlap_start <= CDMConnector::dateadd("cohort_start_date", window[2]))
  }
  xx <- xx %>%
    select("subject_id", "cohort_start_date", "overlap_name") %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = "overlap_name", values_from = "value", values_fill = 0) %>%
    right_join(x, by = c("subject_id", "cohort_start_date"))
  for (nam in name) {
    if (nam %in% colnames(xx)) {
      xx <- xx %>%
        mutate(!!nam := if_else(is.na(.data[[nam]]), 0, .data[[nam]]))
    } else {
      xx <- mutate(xx, !!nam := 0) 
    }
  }
  return(xx %>% compute())
}

addMultipleEvent <- function(x, cdm, eventTableName, eventId = NULL, window = c(NA, NA), name = "event", eventDate = "cohort_start_date", eventAt = "cohort_start_date", order = "first", compute = TRUE) {
  eventTable <- cdm[[eventTableName]]
  if (!is.null(eventId)) {
    eventTable <- eventTable %>% dplyr::filter(.data$cohort_definition_id %in% .env$eventId)
  }
  if ("person_id" %in% colnames(eventTable)) {
    eventTable <- eventTable %>% dplyr::rename("subject_id" = "person_id")
  }
  xx <- x %>%
    dplyr::select(dplyr::all_of(c("subject_id", eventAt))) %>%
    dplyr::inner_join(
      eventTable %>%
        dplyr::select("event_id" = "cohort_definition_id", "subject_id", "event_date" = dplyr::all_of(eventDate)),
      by = "subject_id"
    ) %>%
    dplyr::mutate(dif_time = !!CDMConnector::datediff(eventAt, "event_date"))
  if (!is.na(window[1])) {
    xx <- xx %>%
      dplyr::filter(.data$dif_time >= !!window[1])
  }
  if (!is.na(window[2])) {
    xx <- xx %>%
      dplyr::filter(.data$dif_time <= !!window[2])
  }
  xx <- xx %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c("event_id", "subject_id", eventAt))))
  if (order == "first") {
    xx <- xx %>%
      dplyr::summarise(
        event_date = min(.data$event_date, na.rm = TRUE),
        .groups = "drop"
      )
  } else if (order == "last") {
    xx <- xx %>%
      dplyr::summarise(
        event_date = max(.data$event_date, na.rm = TRUE),
        .groups = "drop"
      )
  }
  xx <- xx %>%
    dplyr::inner_join(
      dplyr::tibble(
        event_id = .env$eventId,
        event_name = .env$name
      ),
      by = "event_id",
      copy = TRUE
    ) %>%
    select(-"event_id") %>%
    tidyr::pivot_wider(names_from = "event_name", values_from = "event_date")
  for (nam in name) {
    if (!(nam %in% colnames(xx))) {
      xx <- xx %>% mutate(!!nam := as.Date(NA))
    }
  }
  xx <- x %>% dplyr::left_join(xx, by = c("subject_id", eventAt))
  if (isTRUE(compute)) {
    xx <- xx %>% dplyr::compute()
  }
  return(xx)
}

firstEventId <- function(x, columns) {
  x <- x %>%
    mutate("first_event" = reduce(select(., all_of(columns)), pmin, na.rm = TRUE)) %>%
    mutate(event_id = as.numeric(NA))
  for (k in 1:length(columns)) {
    x <- x %>% 
      mutate(event_id = if_else(
        is.na(.data$event_id) & .data$first_event == .data[[!!columns[k]]], 
        as.numeric(.env$k), 
        .data$event_id
      ))
  }
  return(x)
}

createSurvivalTable <- function(x, events, censoring, covariates = NULL, start = "cohort_start_date") {
  if (length(censoring) == 1) {
    x <- mutate(x, "censor" = .data[[censoring]])
  } else {
    x <- mutate(x, "censor" = reduce(select(x, all_of(censoring)), pmin, na.rm = TRUE))
  }
  x <- firstEventId(x, c(events, "censor")) %>%
    mutate(event_id = if_else(event_id > !!length(events), 0, event_id)) %>%
    mutate(fine_gray = factor(event_id, 0:length(events), labels=c("censor", events))) %>%
    mutate(cox = as.numeric(if_else(event_id == 1, 1, 0))) %>%
    mutate(time = .data$first_event - .data[[start]]) %>%
    select(all_of(c("group", "weight", "time", "cox", "fine_gray", covariates)))
}

outcomeModel <- function(x, covariates = NULL, unadjusted = FALSE, model = "finegray") {
  
  if (model == "cox") {
    runCox <- TRUE
    runFinegray <- FALSE
  } else if (model == "finegray") {
    runCox <- FALSE
    runFinegray <- TRUE
  } else {
    runCox <- TRUE
    runFinegray <- TRUE
    model <- c("finegray", "cox")
  }
  
  covariates <- c("group", covariates)
  covariates <- paste0(covariates, collapse = " + ")
  finegrayFormula <- as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", covariates))
  coxFormula <- as.formula(paste0("Surv(time, cox) ~ ", covariates))
  
  numberEventsComparator <- x %>%
    filter(group == "comparator") %>%
    pull("cox") %>% 
    sum()
  numberEventsExposure <- x %>%
    filter(group == "exposure") %>%
    pull("cox") %>% 
    sum()
  
  if (numberEventsComparator + numberEventsExposure == 0) {
    
    if (unadjusted) {
      adjustment <- c("overlap weighting", "crude")
    } else {
      adjustment <- "overlap weighting"
    }
    
    result <- expand_grid(
      variable = str_split_1(covariates, " \\+ "),
      model = model,
      adjustment = adjustment,
      coef = as.numeric(NA),
      hr = as.numeric(NA),
      se_coef = as.numeric(NA),
      z = as.numeric(NA),
      p = as.numeric(NA),
      lower_hr = as.numeric(NA),
      upper_hr = as.numeric(NA)
    ) %>%
      mutate(variable = if_else(variable == "group", "groupexposure", variable))
    
  } else {
    
    result <- NULL
    
    # adjusted
    if (runFinegray) {
      fgData <- finegray(Surv(time, fine_gray) ~ ., data=x, weights = weight)
      fgRegression <- coxph(finegrayFormula, weight=fgwt, data=fgData)
      result <- result %>%
        union_all(
          summary(fgRegression)$coefficients %>%
            as_tibble(rownames = "variable") %>%
            select(
              "variable", "coef", "hr" = "exp(coef)", "se_coef" = "se(coef)", "z",
              "p" = "Pr(>|z|)"
            ) %>%
            left_join(
              summary(fgRegression)$conf.int %>%
                as_tibble(rownames = "variable") %>%
                select("variable", "lower_hr" = "lower .95", "upper_hr" = "upper .95"),
              by = "variable"
            ) %>%
            mutate(model = "finegray", adjustment = "overlap weighting")
        )
    }
    
    if (runCox) {
      coxRegression <- coxph(coxFormula, weight=weight, data=x)
      result <- result %>%
        union_all(
          summary(coxRegression)$coefficients %>%
            as_tibble(rownames = "variable") %>%
            select(
              "variable", "coef", "hr" = "exp(coef)", "se_coef" = "se(coef)", "z",
              "p" = "Pr(>|z|)"
            ) %>%
            left_join(
              summary(coxRegression)$conf.int %>%
                as_tibble(rownames = "variable") %>%
                select("variable", "lower_hr" = "lower .95", "upper_hr" = "upper .95"),
              by = "variable"
            ) %>%
            mutate(model = "cox", adjustment = "overlap weighting")
        )
    }
    
    # unadjusted
    if (unadjusted) {
      
      x <- mutate(x, weight = 1)
      
      if (runFinegray) {
        fgData <- finegray(Surv(time, fine_gray) ~ ., data=x, weights = weight)
        fgRegression <- coxph(finegrayFormula, weight=fgwt, data=fgData)
        result <- result %>%
          union_all(
            summary(fgRegression)$coefficients %>%
              as_tibble(rownames = "variable") %>%
              select(
                "variable", "coef", "hr" = "exp(coef)", "se_coef" = "se(coef)", "z",
                "p" = "Pr(>|z|)"
              ) %>%
              left_join(
                summary(fgRegression)$conf.int %>%
                  as_tibble(rownames = "variable") %>%
                  select("variable", "lower_hr" = "lower .95", "upper_hr" = "upper .95"),
                by = "variable"
              ) %>%
              mutate(model = "finegray", adjustment = "crude")
          )
      }
      
      if (runCox) {
        coxRegression <- coxph(coxFormula, weight=weight, data=x)
        result <- result %>%
          union_all(
            summary(coxRegression)$coefficients %>%
              as_tibble(rownames = "variable") %>%
              select(
                "variable", "coef", "hr" = "exp(coef)", "se_coef" = "se(coef)", "z",
                "p" = "Pr(>|z|)"
              ) %>%
              left_join(
                summary(coxRegression)$conf.int %>%
                  as_tibble(rownames = "variable") %>%
                  select("variable", "lower_hr" = "lower .95", "upper_hr" = "upper .95"),
                by = "variable"
              ) %>%
              mutate(model = "cox", adjustment = "crude")
          )
      }
      
    }
    
  }
  
  result <- result %>%
    mutate(number_events_comparator = numberEventsComparator) %>%
    mutate(number_events_exposure = numberEventsExposure)
  
  return(result)
}

collectCohortSummary <- function(cdm, cohortName = NULL) {
  checkmate::assertClass(cdm, "cdm_reference")
  checkmate::assertCharacter(cohortName, any.missing = FALSE, null.ok = TRUE)
  if (is.null(cohortName)) {
    cohortName <- names(cdm)[lapply(names(cdm), function(x) {
      "GeneratedCohortSet" %in% class(cdm[[x]])
    }) %>%
      unlist()]
    if (length(cohortName) == 0) {
      stop("No cohorts in cdm object")
    }
  }
  notCohort <- cohortName[lapply(cohortName, function(x) {
    !("GeneratedCohortSet" %in% class(cdm[[x]]))
  }) %>%
    unlist()]
  if (length(notCohort) != 0) {
    notCohort <- paste0(notCohort, collapse = ", ")
    stop(paste0(notCohort, " are not valid GeneratedCohortSet objects"))
  }
  
  result <- NULL
  for (name in cohortName) {
    set <- cohortSet(cdm[[name]]) %>%
      select("cohort_definition_id", "cohort_name")
    attrition <- cohortAttrition(cdm[[name]]) %>%
      select(
        "cohort_definition_id", "number_records", "number_subjects", "reason", 
        "reason_id", "excluded_records", "excluded_subjects"
      ) %>%
      mutate(across(
        all_of(c(
          "number_records", "number_subjects", "excluded_records", 
          "excluded_subjects"
        )), 
        ~ as.numeric(.)
      ))
    result <- result %>%
      dplyr::union_all(
        set %>%
          dplyr::inner_join(attrition, by = "cohort_definition_id") %>%
          dplyr::collect() %>%
          dplyr::arrange(.data$cohort_definition_id, .data$reason_id) %>%
          dplyr::mutate(cohort_table_name = .env$name)
      )
  }
  result <- result %>% mutate(cdm_name = CDMConnector::cdmName(cdm))
  return(result)
}

insertTable <- function(x,
                        cdm, 
                        name,
                        overwrite = TRUE) {
  con <- attr(cdm, "dbcon")
  writeSchema <- attr(cdm, "write_schema")
  checkTableExist <- name %in% CDMConnector::listTables(con, writeSchema)
  if (checkTableExist) {
    if (overwrite) {
      DBI::dbRemoveTable(con, CDMConnector:::inSchema(writeSchema, name))
    } else {
      stop(paste0("'", name, "' table already exists."))
    }
  }
  DBI::dbCreateTable(con, CDMConnector:::inSchema(writeSchema, name), x)
  DBI::dbAppendTable(con, CDMConnector:::inSchema(writeSchema, name), x)
  if (methods::is(con, "duckdb_connection")) {
    ref <- dplyr::tbl(con, paste(c(writeSchema, name), collapse = "."))
  } else if (length(writeSchema) == 2) {
    ref <- dplyr::tbl(con,
                      dbplyr::in_catalog(writeSchema[[1]], writeSchema[[2]], name))
  } else if (length(writeSchema) == 1) {
    ref <- dplyr::tbl(con, dbplyr::in_schema(writeSchema, name))
  } else {
    ref <- dplyr::tbl(con, name)
  }
  return(ref)
}

getCovariateId <- function(x) {
  if (grepl(")_window:", x)) {
    idWindow <- str_locate(x, "_window:")
    windowName <- substr(x, idWindow[2]+1, nchar(x))
    windowId <- switch(windowName, "-Any;-366" = "3", "-365;-31" = "2", "-30;-1" = "1")
    pId <- str_locate_all(x, "\\(")[[1]]
    conceptId <- substr(x, pId[nrow(pId),1]+1, idWindow[1]-2)
    return(paste0("f", conceptId, "_", windowId))
  } else {
    return(NA)
  }
}

getEstimates <- function(x, outcomes, censors, covariates, unadjusted = FALSE, model = "finegray") {
  result <- NULL
  survivalPlot <- NULL
  for (outcome_name in outcomes) {
    for (censoring_method in censors) {
      censor <- switch(
        censoring_method,
        "leave" = "leave_db",
        "leave+covid" = c("leave_db", "next_covid"),
        "leave+vaccine" = c("leave_db", "next_vaccine"),
        "leave+covid+vaccine" = c("leave_db", "next_covid", "next_vaccine")
      )
      survivalTable <- createSurvivalTable(x, c(outcome_name, "death"), censor, covariates)
      result <- result %>%
        union_all(
          outcomeModel(survivalTable, covariates, unadjusted, model) %>%
            mutate(
              outcome_name = outcome_name,
              censoring_method = censoring_method,
              covariates = paste0(covariates, collapse = "; ")
            )
        )
      survivalPlot <- survivalPlot %>%
        dplyr::union_all(
            survivalTable %>%
          dplyr::mutate(time = as.numeric(time)) %>%
          dplyr::group_by(group, time, fine_gray) %>%
          dplyr::summarise(count = sum(.data$weight), .groups = "drop") %>% 
          dplyr::rename("event" = "fine_gray") %>%
          dplyr::mutate(
            outcome_name = outcome_name,
            censoring_method = censoring_method,
            covariates = paste0(covariates, collapse = "; ")
          )
        )
    }
  }
  x <- list("result" = result, "survival_plot" = survivalPlot)
  return(x)
}

filterInObservation <- function(x, cdm, date) {
  x %>%
    dplyr::select("person_id") %>%
    dplyr::distinct() %>%
    dplyr::inner_join(
      cdm$observation_period %>%
        dplyr::select(
          "person_id", "observation_period_start_date",
          "observation_period_end_date"
        ),
      by = "person_id"
    ) %>%
    dplyr::filter(.data$observation_period_start_date <= .env$date) %>%
    dplyr::filter(.data$observation_period_end_date > .env$date) %>%
    dplyr::select("person_id") %>%
    dplyr::distinct() %>%
    dplyr::anti_join(
      cdm$death %>%
        dplyr::filter(.data$death_date <= .env$date) %>%
        dplyr::select("person_id"),
      by = "person_id"
    ) %>%
    dplyr::inner_join(x, by = "person_id") %>%
    dplyr::compute()
}
