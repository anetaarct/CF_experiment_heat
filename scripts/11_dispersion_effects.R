suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

dir.create(file.path("models", "dispersion_effects"), showWarnings = FALSE, recursive = TRUE)

comparison <- read_csv("models/model_comparisons.csv", show_col_types = FALSE)

response_labels <- comparison %>%
  distinct(response, response_label)

extract_dispersion <- function(response, model_name) {
  candidate_path <- file.path("models", paste0("candidate_models_", response, ".rds"))
  if (!file.exists(candidate_path)) return(tibble())

  fits <- readRDS(candidate_path)
  if (!model_name %in% names(fits)) return(tibble())

  fit <- fits[[model_name]]
  disp_coef <- summary(fit)$coefficients$disp
  if (is.null(disp_coef)) return(tibble())

  tibble::as_tibble(disp_coef, rownames = "term") %>%
    rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `z value`,
      p_value = `Pr(>|z|)`
    ) %>%
    mutate(response = response, model_name = model_name, .before = term)
}

scale_models <- comparison %>%
  filter(grepl("scale", model_name), !is.na(AIC)) %>%
  select(response, model_name)

dispersion_effects <- bind_rows(
  lapply(seq_len(nrow(scale_models)), function(i) {
    extract_dispersion(scale_models$response[i], scale_models$model_name[i])
  })
) %>%
  left_join(response_labels, by = "response") %>%
  relocate(response_label, .after = response)

best_models <- comparison %>%
  filter(!is.na(AIC)) %>%
  group_by(response) %>%
  slice_min(AIC, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(response, best_model = model_name, best_model_has_dispersion = grepl("scale", model_name))

dispersion_effects <- dispersion_effects %>%
  left_join(best_models, by = "response")

write_csv(dispersion_effects, file.path("models", "dispersion_effects", "dispersion_coefficients.csv"))
write_csv(best_models, file.path("models", "dispersion_effects", "dispersion_best_model_flags.csv"))

message("Dispersion-effect tables saved in models/dispersion_effects/.")
