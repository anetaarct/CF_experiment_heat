suppressPackageStartupMessages({
  library(glmmTMB)
  library(dplyr)
  library(ggplot2)
  library(ggeffects)
  library(readr)
})

dir.create(file.path("figures", "experiment_interactions"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path("figures", "temperature_interactions"), recursive = TRUE, showWarnings = FALSE)

dat <- readRDS("data/prepared_nestlings.rds")

trait_labels <- c(
  D8_MASS = "Body mass day 8",
  D12_MASS = "Body mass day 12",
  D12_RTARS = "Tarsus length day 12",
  Early.g = "Early growth",
  Late.g = "Late growth",
  Growth = "Whole-period growth",
  D12_SMI = "Scaled mass index day 12",
  D12_SURV = "Survival to day 12"
)

responses <- c(
  "D8_MASS",
  "D12_MASS",
  "Early.g",
  "Late.g",
  "Growth",
  "D12_SMI",
  "D12_RTARS",
  "D12_SURV"
)

`%||%` <- function(x, y) if (is.null(x)) y else x

theme_manuscript <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.caption = element_text(size = 9, color = "grey35", hjust = 0)
    )
}

load_interaction_fit <- function(path) {
  fits <- readRDS(path)
  interaction_name <- grep("^interaction", names(fits), value = TRUE)[1]
  fits[[interaction_name]]
}

plot_experiment_interaction <- function(response) {
  model_file <- file.path("models", "experiment_interaction", paste0("experiment_interaction_models_", response, ".rds"))
  if (!file.exists(model_file)) return(invisible(NULL))

  fit <- load_interaction_fit(model_file)
  y_label <- trait_labels[[response]] %||% response

  if (response == "D12_SURV") {
    raw_data <- dat %>%
      filter(!is.na(D12_SURV), !is.na(EXP.INC), !is.na(EXP.NEST)) %>%
      mutate(D12_SURV_num = as.numeric(as.character(D12_SURV)))
    pred <- ggpredict(fit, terms = c("EXP.INC", "EXP.NEST"))

    p <- ggplot() +
      geom_jitter(
        data = raw_data,
        aes(x = EXP.INC, y = D12_SURV_num),
        width = 0.12,
        height = 0.04,
        alpha = 0.16,
        size = 0.9,
        color = "grey35"
      ) +
      geom_line(
        data = pred,
        aes(x = x, y = predicted, color = group, group = group),
        linewidth = 0.85
      ) +
      geom_pointrange(
        data = pred,
        aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group, group = group),
        linewidth = 0.5,
        position = position_dodge(width = 0.18)
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        x = "Incubation treatment",
        y = "Predicted probability of survival to day 12",
        color = "Nestling treatment",
        caption = "Pale jittered points show original 0/1 outcomes; dark points and lines show population-level model predictions with 95% confidence intervals."
      )
  } else {
    raw_data <- dat %>% filter(!is.na(.data[[response]]), !is.na(EXP.INC), !is.na(EXP.NEST))
    pred <- ggpredict(fit, terms = c("EXP.INC", "EXP.NEST"))

    p <- ggplot() +
      geom_jitter(
        data = raw_data,
        aes(x = EXP.INC, y = .data[[response]]),
        width = 0.13,
        height = 0,
        alpha = 0.14,
        size = 0.9,
        color = "grey35"
      ) +
      geom_line(
        data = pred,
        aes(x = x, y = predicted, color = group, group = group),
        linewidth = 0.85
      ) +
      geom_pointrange(
        data = pred,
        aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group, group = group),
        linewidth = 0.5,
        position = position_dodge(width = 0.18)
      ) +
      labs(
        x = "Incubation treatment",
        y = paste("Predicted", y_label),
        color = "Nestling treatment",
        caption = "Pale points show original observations; dark points and lines show population-level model predictions with 95% confidence intervals."
      )
  }

  p <- p +
    scale_color_manual(values = c("#476F67", "#C98C1E", "#6E93A2", "#8FAEA3")) +
    theme_manuscript()

  stem <- file.path("figures", "experiment_interactions", paste0("experiment_interaction_", response))
  ggsave(paste0(stem, ".png"), p, width = 7.2, height = 4.8, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 7.2, height = 4.8)
  invisible(p)
}

plot_temperature_interaction <- function(response) {
  model_file <- file.path("models", "sensitivity_temperature_interaction", paste0("temperature_interaction_models_", response, ".rds"))
  if (!file.exists(model_file)) return(invisible(NULL))

  fit <- load_interaction_fit(model_file)
  y_label <- trait_labels[[response]] %||% response

  if (response == "D12_SURV") {
    raw_data <- dat %>%
      filter(!is.na(D12_SURV), !is.na(incubation_temperature_sc), !is.na(nestling_temperature_sc)) %>%
      mutate(D12_SURV_num = as.numeric(as.character(D12_SURV)))
    pred <- ggpredict(fit, terms = c("incubation_temperature_sc [all]", "nestling_temperature_sc [-1,0,1]"))

    p <- ggplot() +
      geom_jitter(
        data = raw_data,
        aes(x = incubation_temperature_sc, y = D12_SURV_num),
        width = 0.03,
        height = 0.04,
        alpha = 0.14,
        size = 0.8,
        color = "grey35"
      ) +
      geom_ribbon(
        data = pred,
        aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
        alpha = 0.16,
        color = NA
      ) +
      geom_line(
        data = pred,
        aes(x = x, y = predicted, color = group),
        linewidth = 0.9
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        x = "Incubation temperature (standardized log mean)",
        y = "Predicted probability of survival to day 12",
        color = "Nestling temperature",
        fill = "Nestling temperature",
        caption = "Lines show population-level predictions across incubation temperature at low, average, and high nestling-period temperature; ribbons are 95% confidence intervals."
      )
  } else {
    raw_data <- dat %>% filter(!is.na(.data[[response]]), !is.na(incubation_temperature_sc), !is.na(nestling_temperature_sc))
    pred <- ggpredict(fit, terms = c("incubation_temperature_sc [all]", "nestling_temperature_sc [-1,0,1]"))

    p <- ggplot() +
      geom_point(
        data = raw_data,
        aes(x = incubation_temperature_sc, y = .data[[response]]),
        alpha = 0.14,
        size = 0.8,
        color = "grey35"
      ) +
      geom_ribbon(
        data = pred,
        aes(x = x, ymin = conf.low, ymax = conf.high, fill = group),
        alpha = 0.16,
        color = NA
      ) +
      geom_line(
        data = pred,
        aes(x = x, y = predicted, color = group),
        linewidth = 0.9
      ) +
      labs(
        x = "Incubation temperature (standardized log mean)",
        y = paste("Predicted", y_label),
        color = "Nestling temperature",
        fill = "Nestling temperature",
        caption = "Lines show population-level predictions across incubation temperature at low, average, and high nestling-period temperature; ribbons are 95% confidence intervals and pale points are original observations."
      )
  }

  p <- p +
    scale_color_manual(values = c("#476F67", "#C98C1E", "#6E93A2")) +
    scale_fill_manual(values = c("#476F67", "#C98C1E", "#6E93A2")) +
    theme_manuscript()

  stem <- file.path("figures", "temperature_interactions", paste0("temperature_interaction_", response))
  ggsave(paste0(stem, ".png"), p, width = 7.2, height = 4.8, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 7.2, height = 4.8)
  invisible(p)
}

for (response in responses) {
  plot_experiment_interaction(response)
  plot_temperature_interaction(response)
}

message("Secondary figures saved in figures/experiment_interactions and figures/temperature_interactions.")
