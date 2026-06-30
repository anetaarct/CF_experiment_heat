suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(readr)
  library(performance)
  library(DHARMa)
})

source("scripts/01_prepare_data.R")

dir.create("models", showWarnings = FALSE)
dir.create(file.path("models", "summary_pdf"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("models", "diagnostics"), showWarnings = FALSE, recursive = TRUE)

dat <- readRDS("data/prepared_nestlings.rds")

response_labels <- c(
  D2_MASS = "Day-2 body mass",
  D8_MASS = "Day-8 body mass",
  D12_MASS = "Day-12 body mass",
  Early.g = "Early growth",
  Late.g = "Late growth",
  Growth = "Whole-period growth",
  D12_SMI = "Day-12 scaled mass index",
  D12_RTARS = "Day-12 tarsus length",
  D12_SURV = "Day-12 survival"
)

model_specs <- tibble::tribble(
  ~response, ~label, ~use_nestling_phase,
  "D2_MASS", "Day-2 body mass", FALSE,
  "D8_MASS", "Day-8 body mass", TRUE,
  "D12_MASS", "Day-12 body mass", TRUE,
  "Early.g", "Early growth", TRUE,
  "Late.g", "Late growth", TRUE,
  "Growth", "Whole-period growth", TRUE,
  "D12_SMI", "Day-12 scaled mass index", TRUE,
  "D12_RTARS", "Day-12 tarsus length", TRUE
)

as_clean_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

write_text_pdf <- function(lines, file, title) {
  grDevices::pdf(file, width = 8.27, height = 11.69, paper = "a4")
  on.exit(grDevices::dev.off(), add = TRUE)
  lines <- enc2utf8(lines)
  lines <- c(title, "", lines)
  per_page <- 52
  pages <- split(lines, ceiling(seq_along(lines) / per_page))
  for (page in pages) {
    graphics::plot.new()
    graphics::par(mar = c(0.5, 0.5, 0.5, 0.5))
    y <- seq(0.97, 0.04, length.out = length(page))
    graphics::text(0.02, y, labels = page, adj = c(0, 1), family = "mono", cex = 0.62)
  }
}

capture_model_artifacts <- function(fit, response, model_name, label) {
  stem <- paste(response, model_name, sep = "_")
  summary_lines <- capture.output(summary(fit))
  txt_file <- file.path("models", paste0("summary_", stem, ".txt"))
  pdf_file <- file.path("models", "summary_pdf", paste0("summary_", stem, ".pdf"))
  writeLines(summary_lines, txt_file, useBytes = TRUE)
  write_text_pdf(summary_lines, pdf_file, paste(label, "-", model_name))

  diag_pdf <- file.path("models", "diagnostics", paste0("diagnostic_", stem, ".pdf"))
  grDevices::pdf(diag_pdf, width = 8, height = 8)
  try({
    residuals_sim <- DHARMa::simulateResiduals(fit)
    plot(residuals_sim)
  }, silent = TRUE)
  grDevices::dev.off()

  invisible(list(summary_txt = txt_file, summary_pdf = pdf_file, diagnostic_pdf = diag_pdf))
}

fixed_effect_table <- function(fit, response, model_name, label, selected) {
  coef_mat <- summary(fit)$coefficients$cond
  tibble::as_tibble(coef_mat, rownames = "term") %>%
    rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `z value`,
      p_value = `Pr(>|z|)`
    ) %>%
    mutate(
      response = response,
      response_label = label,
      model_name = model_name,
      selected_model = selected,
      .before = term
    )
}

fit_continuous_set <- function(response, label, use_nestling_phase) {
  model_data <- dat %>% filter(!is.na(.data[[response]])) %>% droplevels()

  if (use_nestling_phase) {
    mean_formula <- as.formula(paste0(
      response,
      " ~ EXP.INC + EXP.NEST + SEX + BS_sc + (1|F_RING) + (1|YEAR)"
    ))
    disp_formula <- ~ EXP.INC + EXP.NEST
  } else {
    mean_formula <- as.formula(paste0(
      response,
      " ~ EXP.INC + SEX + BS_sc + (1|F_RING) + (1|YEAR)"
    ))
    disp_formula <- ~ EXP.INC
  }

  fits <- list(
    gaussian_mean = glmmTMB(mean_formula, data = model_data),
    gaussian_scale = glmmTMB(mean_formula, dispformula = disp_formula, data = model_data),
    student_mean = glmmTMB(mean_formula, family = t_family(), data = model_data),
    student_scale = glmmTMB(mean_formula, family = t_family(), dispformula = disp_formula, data = model_data)
  )

  comparison <- tibble(
    response = response,
    response_label = label,
    model_name = names(fits),
    model_family = c("Gaussian", "Gaussian", "Student-t", "Student-t"),
    dispersion_formula = c("~ 1", deparse(disp_formula), "~ 1", deparse(disp_formula)),
    AIC = as.numeric(AIC(fits$gaussian_mean, fits$gaussian_scale, fits$student_mean, fits$student_scale)$AIC)
  ) %>%
    arrange(AIC) %>%
    mutate(delta_AIC = if (all(is.na(AIC))) NA_real_ else AIC - min(AIC, na.rm = TRUE))

  selected_model <- "gaussian_mean"

  comparison <- comparison %>% mutate(selected_model = model_name == selected_model)
  best_fit <- fits[[selected_model]]

  saveRDS(fits, file.path("models", paste0("candidate_models_", response, ".rds")))
  saveRDS(best_fit, file.path("models", paste0("best_model_", response, ".rds")))
  write_csv(comparison, file.path("models", paste0("model_comparison_", response, ".csv")))

  artifact_rows <- list()
  coef_rows <- list()
  for (nm in names(fits)) {
    artifacts <- capture_model_artifacts(fits[[nm]], response, nm, label)
    artifact_rows[[nm]] <- tibble(
      response = response,
      response_label = label,
      model_name = nm,
      selected_model = nm == selected_model,
      summary_txt = artifacts$summary_txt,
      summary_pdf = artifacts$summary_pdf,
      diagnostic_pdf = artifacts$diagnostic_pdf
    )
    coef_rows[[nm]] <- fixed_effect_table(fits[[nm]], response, nm, label, nm == selected_model)
  }

  list(
    comparison = comparison,
    artifacts = bind_rows(artifact_rows),
    coefficients = bind_rows(coef_rows),
    selected = tibble(
      response = response,
      response_label = label,
      selected_model = selected_model,
      formula = paste(deparse(formula(best_fit)), collapse = " "),
      dispersion_formula = paste(deparse(best_fit$modelInfo$allForm$dispformula %||% ~1), collapse = " "),
      n = nrow(model_data)
    )
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

results <- lapply(seq_len(nrow(model_specs)), function(i) {
  fit_continuous_set(
    model_specs$response[i],
    model_specs$label[i],
    model_specs$use_nestling_phase[i]
  )
})

surv_dat <- dat %>%
  filter(!is.na(D12_SURV)) %>%
  mutate(D12_SURV = factor(D12_SURV, levels = c(0, 1), labels = c("Dead", "Survived"))) %>%
  droplevels()

model_surv <- glmmTMB(
  D12_SURV ~ EXP.INC + EXP.NEST + SEX + BS_sc + (1|F_RING) + (1|YEAR),
  family = binomial(),
  data = surv_dat
)

saveRDS(model_surv, "models/best_model_D12_SURV.rds")
surv_artifacts <- capture_model_artifacts(model_surv, "D12_SURV", "binomial", response_labels[["D12_SURV"]])
surv_coef <- fixed_effect_table(model_surv, "D12_SURV", "binomial", response_labels[["D12_SURV"]], TRUE)
surv_selected <- tibble(
  response = "D12_SURV",
  response_label = response_labels[["D12_SURV"]],
  selected_model = "binomial",
  formula = paste(deparse(formula(model_surv)), collapse = " "),
  dispersion_formula = "not used",
  n = nrow(surv_dat)
)
surv_artifact_row <- tibble(
  response = "D12_SURV",
  response_label = response_labels[["D12_SURV"]],
  model_name = "binomial",
  selected_model = TRUE,
  summary_txt = surv_artifacts$summary_txt,
  summary_pdf = surv_artifacts$summary_pdf,
  diagnostic_pdf = surv_artifacts$diagnostic_pdf
)
surv_comparison <- tibble(
  response = "D12_SURV",
  response_label = response_labels[["D12_SURV"]],
  model_name = "binomial",
  model_family = "Binomial",
  dispersion_formula = "not used",
  AIC = AIC(model_surv),
  delta_AIC = 0,
  selected_model = TRUE
)

model_comparisons <- bind_rows(lapply(results, `[[`, "comparison"), surv_comparison)
model_artifacts <- bind_rows(lapply(results, `[[`, "artifacts"), surv_artifact_row)
model_coefficients <- bind_rows(lapply(results, `[[`, "coefficients"), surv_coef)
selected_models <- bind_rows(lapply(results, `[[`, "selected"), surv_selected)

write_csv(model_comparisons, "models/model_comparisons.csv")
write_csv(model_artifacts, "models/model_artifacts.csv")
write_csv(model_coefficients, "models/model_coefficients.csv")
write_csv(selected_models, "models/selected_models.csv")
saveRDS(selected_models, "models/selected_models.rds")

message("Model fitting finished. Summaries, PDFs, diagnostics, and tables saved in models/.")
