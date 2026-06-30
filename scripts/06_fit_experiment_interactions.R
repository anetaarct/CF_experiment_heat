suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(readr)
  library(performance)
})

source("scripts/01_prepare_data.R")

dir.create("models/experiment_interaction", recursive = TRUE, showWarnings = FALSE)

dat <- readRDS("data/prepared_nestlings.rds")

response_labels <- c(
  D8_MASS = "Day-8 body mass",
  D12_MASS = "Day-12 body mass",
  Early.g = "Early growth",
  Late.g = "Late growth",
  Growth = "Whole-period growth",
  D12_SMI = "Day-12 scaled mass index",
  D12_RTARS = "Day-12 tarsus length",
  D12_SURV = "Day-12 survival"
)

later_responses <- c(
  "D8_MASS",
  "D12_MASS",
  "Early.g",
  "Late.g",
  "Growth",
  "D12_SMI",
  "D12_RTARS"
)

extract_experiment_interaction <- function(fit, response, response_label, model_name) {
  coef_tab <- summary(fit)$coefficients$cond
  term <- grep("EXP\\.INC.*:EXP\\.NEST|EXP\\.NEST.*:EXP\\.INC", rownames(coef_tab), value = TRUE)

  if (length(term) == 0) {
    return(tibble(
      response = response,
      response_label = response_label,
      model_name = model_name,
      term = "EXP.INC:EXP.NEST",
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    ))
  }

  term <- term[1]
  tibble(
    response = response,
    response_label = response_label,
    model_name = model_name,
    term = term,
    estimate = coef_tab[term, "Estimate"],
    std_error = coef_tab[term, "Std. Error"],
    statistic = coef_tab[term, "z value"],
    p_value = coef_tab[term, "Pr(>|z|)"]
  )
}

fit_experiment_interaction <- function(data, response) {
  response_label <- response_labels[[response]]
  model_data <- data %>%
    filter(
      !is.na(.data[[response]]),
      !is.na(EXP.INC),
      !is.na(EXP.NEST),
      !is.na(SEX),
      !is.na(BS_sc)
    ) %>%
    droplevels()

  additive_formula <- as.formula(paste0(
    response,
    " ~ EXP.INC + EXP.NEST + SEX + BS_sc + (1|F_RING) + (1|YEAR)"
  ))

  interaction_formula <- as.formula(paste0(
    response,
    " ~ EXP.INC * EXP.NEST + SEX + BS_sc + (1|F_RING) + (1|YEAR)"
  ))

  fits <- list(
    additive_gaussian = glmmTMB(additive_formula, data = model_data),
    interaction_gaussian = glmmTMB(interaction_formula, data = model_data)
  )

  comparison <- performance::compare_performance(
    fits$additive_gaussian,
    fits$interaction_gaussian
  ) %>%
    as.data.frame()

  names(comparison)[names(comparison) == "Name"] <- "candidate_model"
  comparison$response <- response
  comparison$response_label <- response_label
  comparison$n_observations <- nrow(model_data)

  saveRDS(fits, file.path("models/experiment_interaction", paste0("experiment_interaction_models_", response, ".rds")))
  write_csv(comparison, file.path("models/experiment_interaction", paste0("experiment_interaction_comparison_", response, ".csv")))

  for (nm in names(fits)) {
    capture.output(
      summary(fits[[nm]]),
      file = file.path("models/experiment_interaction", paste0("summary_", response, "_", nm, ".txt"))
    )
  }

  extract_experiment_interaction(fits$interaction_gaussian, response, response_label, "interaction_gaussian")
}

interaction_terms <- lapply(later_responses, function(response) {
  fit_experiment_interaction(dat, response)
})

surv_dat <- dat %>%
  filter(
    !is.na(D12_SURV),
    !is.na(EXP.INC),
    !is.na(EXP.NEST),
    !is.na(SEX),
    !is.na(BS_sc)
  ) %>%
  mutate(D12_SURV_bin = as.numeric(as.character(D12_SURV))) %>%
  filter(D12_SURV_bin %in% c(0, 1)) %>%
  droplevels()

surv_additive <- glmmTMB(
  D12_SURV_bin ~ EXP.INC + EXP.NEST + SEX + BS_sc + (1|F_RING) + (1|YEAR),
  family = binomial(),
  data = surv_dat
)

surv_interaction <- glmmTMB(
  D12_SURV_bin ~ EXP.INC * EXP.NEST + SEX + BS_sc + (1|F_RING) + (1|YEAR),
  family = binomial(),
  data = surv_dat
)

surv_comparison <- performance::compare_performance(surv_additive, surv_interaction) %>%
  as.data.frame()
names(surv_comparison)[names(surv_comparison) == "Name"] <- "candidate_model"
surv_comparison$response <- "D12_SURV"
surv_comparison$response_label <- response_labels[["D12_SURV"]]
surv_comparison$n_observations <- nrow(surv_dat)

saveRDS(
  list(additive_binomial = surv_additive, interaction_binomial = surv_interaction),
  "models/experiment_interaction/experiment_interaction_models_D12_SURV.rds"
)
write_csv(surv_comparison, "models/experiment_interaction/experiment_interaction_comparison_D12_SURV.csv")
capture.output(summary(surv_additive), file = "models/experiment_interaction/summary_D12_SURV_additive_binomial.txt")
capture.output(summary(surv_interaction), file = "models/experiment_interaction/summary_D12_SURV_interaction_binomial.txt")

interaction_terms <- bind_rows(
  interaction_terms,
  extract_experiment_interaction(surv_interaction, "D12_SURV", response_labels[["D12_SURV"]], "interaction_binomial")
)

write_csv(interaction_terms, "models/experiment_interaction/experiment_interaction_terms.csv")

message("Experiment-interaction models finished. Outputs saved in models/experiment_interaction/.")
