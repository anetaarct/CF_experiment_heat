suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(readr)
  library(emmeans)
})

dir.create(file.path("models", "group_contrasts"), showWarnings = FALSE, recursive = TRUE)

nestlings <- readRDS("data/prepared_nestlings.rds") %>%
  mutate(GROUP = factor(GROUP, levels = c("CONCON", "CONEXP", "EXPCON", "EXPEXP")))

hatchability <- readRDS("data/prepared_hatchability_eggs.rds") %>%
  mutate(GROUP = factor(GROUP, levels = c("CONCON", "CONEXP", "EXPCON", "EXPEXP")))

response_labels <- c(
  HATCHED = "Hatchability",
  D2_MASS = "Body mass at day 2",
  D8_MASS = "Body mass at day 8",
  D12_MASS = "Body mass at day 12",
  Early.g = "Early growth",
  Late.g = "Late growth",
  Growth = "Whole-period growth",
  D12_SMI = "Scaled mass index at day 12",
  D12_RTARS = "Tarsus length at day 12",
  D12_SURV = "Survival to day 12"
)

fit_specs <- tibble::tribble(
  ~response, ~family, ~data_type,
  "HATCHED", "binomial", "egg",
  "D2_MASS", "gaussian", "nestling",
  "D8_MASS", "gaussian", "nestling",
  "D12_MASS", "gaussian", "nestling",
  "Early.g", "gaussian", "nestling",
  "Late.g", "gaussian", "nestling",
  "Growth", "gaussian", "nestling",
  "D12_SMI", "gaussian", "nestling",
  "D12_RTARS", "gaussian", "nestling",
  "D12_SURV", "binomial", "nestling"
)

fit_group_model <- function(response, family, data_type) {
  label <- response_labels[[response]]

  if (data_type == "egg") {
    model_data <- hatchability %>%
      filter(!is.na(.data[[response]]), !is.na(GROUP), !is.na(CS_sc), !is.na(F_RING), !is.na(YEAR)) %>%
      droplevels()
    formula <- as.formula(paste0(response, " ~ GROUP + CS_sc + (1|F_RING) + (1|YEAR)"))
  } else if (family == "binomial") {
    model_data <- nestlings %>%
      filter(!is.na(.data[[response]]), !is.na(GROUP), !is.na(SEX), !is.na(BS_sc), !is.na(F_RING), !is.na(YEAR)) %>%
      mutate("{response}" := factor(.data[[response]], levels = c(0, 1), labels = c("No", "Yes"))) %>%
      droplevels()
    formula <- as.formula(paste0(response, " ~ GROUP + SEX + BS_sc + (1|F_RING) + (1|YEAR)"))
  } else {
    model_data <- nestlings %>%
      filter(!is.na(.data[[response]]), !is.na(GROUP), !is.na(SEX), !is.na(BS_sc), !is.na(F_RING), !is.na(YEAR)) %>%
      droplevels()
    formula <- as.formula(paste0(response, " ~ GROUP + SEX + BS_sc + (1|F_RING) + (1|YEAR)"))
  }

  fit <- if (family == "binomial") {
    glmmTMB(formula, family = binomial(), data = model_data)
  } else {
    glmmTMB(formula, data = model_data)
  }

  emm <- emmeans(fit, ~ GROUP, type = if (family == "binomial") "response" else "link")
  means <- as.data.frame(emm) %>%
    as_tibble() %>%
    mutate(
      response = response,
      response_label = label,
      n = nrow(model_data),
      .before = GROUP
    )

  contrasts <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    as_tibble() %>%
    mutate(
      response = response,
      response_label = label,
      n = nrow(model_data),
      .before = contrast
    )

  saveRDS(fit, file.path("models", "group_contrasts", paste0("group_model_", response, ".rds")))
  capture.output(summary(fit), file = file.path("models", "group_contrasts", paste0("summary_group_model_", response, ".txt")))

  list(means = means, contrasts = contrasts)
}

results <- lapply(seq_len(nrow(fit_specs)), function(i) {
  fit_group_model(fit_specs$response[i], fit_specs$family[i], fit_specs$data_type[i])
})

group_means <- bind_rows(lapply(results, `[[`, "means"))
group_contrasts <- bind_rows(lapply(results, `[[`, "contrasts"))

write_csv(group_means, file.path("models", "group_contrasts", "group_estimated_means.csv"))
write_csv(group_contrasts, file.path("models", "group_contrasts", "group_pairwise_contrasts.csv"))

message("Treatment-group contrast sensitivity outputs saved in models/group_contrasts/.")
