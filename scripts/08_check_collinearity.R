suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(glmmTMB)
  library(performance)
})

dir.create(file.path("models", "collinearity"), recursive = TRUE, showWarnings = FALSE)

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

temperature_correlation <- dat %>%
  summarise(
    n_complete = sum(complete.cases(incubation_temperature_sc, nestling_temperature_sc)),
    pearson_r = cor(incubation_temperature_sc, nestling_temperature_sc, use = "complete.obs")
  )

write_csv(temperature_correlation, "models/collinearity/temperature_correlation.csv")

continuous_responses <- c(
  "D8_MASS",
  "D12_MASS",
  "Early.g",
  "Late.g",
  "Growth",
  "D12_SMI",
  "D12_RTARS"
)

as_vif_table <- function(vif_object, response, response_label) {
  as.data.frame(vif_object) %>%
    as_tibble() %>%
    mutate(
      response = response,
      response_label = response_label,
      .before = 1
    )
}

fit_temperature_vif <- function(response) {
  model_data <- dat %>%
    filter(
      !is.na(.data[[response]]),
      !is.na(EXP.INC),
      !is.na(EXP.NEST),
      !is.na(incubation_temperature_sc),
      !is.na(nestling_temperature_sc),
      !is.na(SEX),
      !is.na(BS_sc)
    ) %>%
    droplevels()

  model_formula <- as.formula(paste0(
    response,
    " ~ EXP.INC + EXP.NEST + incubation_temperature_sc * nestling_temperature_sc + SEX + BS_sc + (1|F_RING) + (1|YEAR)"
  ))

  fit <- glmmTMB(model_formula, data = model_data)
  as_vif_table(
    performance::check_collinearity(fit),
    response,
    response_labels[[response]]
  )
}

vif_tables <- lapply(continuous_responses, fit_temperature_vif)

surv_dat <- dat %>%
  filter(
    !is.na(D12_SURV),
    !is.na(EXP.INC),
    !is.na(EXP.NEST),
    !is.na(incubation_temperature_sc),
    !is.na(nestling_temperature_sc),
    !is.na(SEX),
    !is.na(BS_sc)
  ) %>%
  mutate(D12_SURV_bin = as.numeric(as.character(D12_SURV))) %>%
  filter(D12_SURV_bin %in% c(0, 1)) %>%
  droplevels()

surv_fit <- glmmTMB(
  D12_SURV_bin ~ EXP.INC + EXP.NEST + incubation_temperature_sc * nestling_temperature_sc + SEX + BS_sc + (1|F_RING) + (1|YEAR),
  family = binomial(),
  data = surv_dat
)

vif_tables <- bind_rows(
  vif_tables,
  as_vif_table(
    performance::check_collinearity(surv_fit),
    "D12_SURV",
    response_labels[["D12_SURV"]]
  )
)

write_csv(vif_tables, "models/collinearity/temperature_model_vif.csv")

print(temperature_correlation)
print(vif_tables)

message("Collinearity diagnostics saved in models/collinearity/.")
