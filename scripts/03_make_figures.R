suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(ggeffects)
})

dir.create("figures", showWarnings = FALSE)
dir.create(file.path("figures", "model_effects"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path("figures", "raw_distributions"), showWarnings = FALSE, recursive = TRUE)

dat <- readRDS("data/prepared_nestlings.rds")
selected_models <- read_csv("models/selected_models.csv", show_col_types = FALSE)

trait_labels <- c(
  D2_MASS = "Body mass day 2",
  D8_MASS = "Body mass day 8",
  D12_MASS = "Body mass day 12",
  D12_RTARS = "Tarsus length day 12",
  Early.g = "Early growth",
  Late.g = "Late growth",
  Growth = "Whole-period growth",
  D12_SMI = "Scaled mass index day 12",
  D12_SURV = "Survival to day 12"
)

continuous_responses <- setdiff(selected_models$response, "D12_SURV")

plot_raw_distribution <- function(response) {
  plot_data <- dat %>% filter(!is.na(.data[[response]]))
  y_label <- trait_labels[[response]] %||% response

  if (response == "D2_MASS") {
    p <- ggplot(plot_data, aes(x = EXP.INC, y = .data[[response]], fill = EXP.INC)) +
      geom_boxplot(alpha = 0.72, outlier.alpha = 0.25) +
      geom_jitter(width = 0.12, alpha = 0.22, size = 1) +
      labs(x = "Incubation treatment", y = y_label)
  } else {
    p <- ggplot(plot_data, aes(x = EXP.INC, y = .data[[response]], fill = EXP.INC)) +
      geom_boxplot(alpha = 0.72, outlier.alpha = 0.25) +
      geom_jitter(width = 0.12, alpha = 0.22, size = 1) +
      facet_wrap(~ EXP.NEST) +
      labs(x = "Incubation treatment", y = y_label, fill = "Incubation")
  }

  p <- p +
    scale_fill_manual(values = c("#476F67", "#C98C1E", "#8FAEA3", "#6E93A2")) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")

  stem <- file.path("figures", "raw_distributions", paste0("raw_", response))
  ggsave(paste0(stem, ".png"), p, width = 7.2, height = 4.8, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 7.2, height = 4.8)
  invisible(p)
}

plot_model_effect <- function(response) {
  fit <- readRDS(file.path("models", paste0("best_model_", response, ".rds")))
  y_label <- trait_labels[[response]] %||% response

  if (response == "D2_MASS") {
    raw_data <- dat %>% filter(!is.na(.data[[response]]), !is.na(EXP.INC))
    pred <- ggpredict(fit, terms = "EXP.INC")

    p <- ggplot() +
      geom_jitter(
        data = raw_data,
        aes(x = EXP.INC, y = .data[[response]]),
        width = 0.12,
        height = 0,
        alpha = 0.16,
        size = 1,
        color = "grey35"
      ) +
      geom_pointrange(
        data = pred,
        aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
        size = 0.62,
        linewidth = 0.55,
        color = "#1F3A3D"
      ) +
      labs(
        x = "Incubation treatment",
        y = paste("Predicted", y_label),
        caption = "Pale points show individual observations; dark points and vertical intervals show model-predicted means with 95% confidence intervals."
      )
  } else if (response == "D12_SURV") {
    raw_prop <- dat %>%
      filter(!is.na(D12_SURV), !is.na(EXP.INC), !is.na(EXP.NEST)) %>%
      mutate(D12_SURV_num = as.numeric(as.character(D12_SURV))) %>%
      group_by(EXP.INC, EXP.NEST) %>%
      summarise(observed_survival = mean(D12_SURV_num, na.rm = TRUE), .groups = "drop")

    pred <- ggpredict(fit, terms = c("EXP.INC", "EXP.NEST"))

    p <- ggplot() +
      geom_col(
        data = raw_prop,
        aes(x = EXP.INC, y = observed_survival, fill = EXP.NEST),
        position = position_dodge(width = 0.45),
        width = 0.32,
        alpha = 0.18,
        color = NA
      ) +
      geom_line(
        data = pred,
        aes(x = x, y = predicted, color = group, group = group),
        linewidth = 0.85
      ) +
      geom_pointrange(
        data = pred,
        aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high, color = group, group = group),
        size = 0.58,
        linewidth = 0.5,
        position = position_dodge(width = 0.18)
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      labs(
        x = "Incubation treatment",
        y = "Predicted survival",
        color = "Nestling treatment",
        fill = "Observed nestling treatment",
        caption = "Pale bars show observed survival proportions; dark points/lines and vertical intervals show model-predicted survival with 95% confidence intervals."
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
        alpha = 0.13,
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
        size = 0.58,
        linewidth = 0.5,
        position = position_dodge(width = 0.18)
      ) +
      labs(
        x = "Incubation treatment",
        y = y_label,
        color = "Nestling-stage treatment",
        caption = "Pale points show individual observations. Dark points and lines show model-predicted means with 95% confidence intervals."
      )
  }

  p <- p +
    scale_color_manual(values = c("#476F67", "#C98C1E", "#6E93A2", "#8FAEA3")) +
    scale_fill_manual(values = c("#476F67", "#C98C1E", "#6E93A2", "#8FAEA3")) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      plot.caption = element_text(size = 9, color = "grey35", hjust = 0)
    )

  stem <- file.path("figures", "model_effects", paste0("model_effect_", response))
  ggsave(paste0(stem, ".png"), p, width = 7.2, height = 4.8, dpi = 300)
  ggsave(paste0(stem, ".pdf"), p, width = 7.2, height = 4.8)
  invisible(p)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

for (response in continuous_responses) {
  plot_raw_distribution(response)
  plot_model_effect(response)
}

message("Figures saved in figures/raw_distributions and figures/model_effects.")
