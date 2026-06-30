suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(glmmTMB)
  library(ggeffects)
  library(patchwork)
})

source("scripts/01_prepare_data.R")

dir.create(file.path("figures", "publication"), showWarnings = FALSE, recursive = TRUE)

dat <- readRDS("data/prepared_nestlings.rds") %>%
  filter(!is.na(D12_SURV)) %>%
  mutate(
    D12_SURV = as.numeric(as.character(D12_SURV)),
    EXP.INC = droplevels(EXP.INC),
    EXP.NEST = droplevels(EXP.NEST),
    SEX = droplevels(SEX),
    YEAR = droplevels(YEAR),
    F_RING = droplevels(F_RING)
  ) %>%
  filter(D12_SURV %in% c(0, 1))

survival_model <- glmmTMB(
  D12_SURV ~ EXP.INC + EXP.NEST + incubation_temperature_sc + nestling_temperature_sc +
    SEX + BS_sc + (1 | F_RING) + (1 | YEAR),
  family = binomial(),
  data = dat
)

saveRDS(survival_model, file.path("models", "publication_model_D12_SURV.rds"))

prediction_theme <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black"),
    legend.position = "bottom",
    legend.title = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
    plot.tag = element_text(face = "bold", size = 13),
    plot.margin = margin(6, 8, 6, 8)
  )

panel_a_pred <- as.data.frame(ggpredict(
  survival_model,
  terms = "incubation_temperature_sc [n=100]",
  type = "fixed"
))

panel_a <- ggplot() +
  geom_jitter(
    data = dat %>% filter(!is.na(incubation_temperature_sc)),
    aes(x = incubation_temperature_sc, y = D12_SURV),
    width = 0,
    height = 0.045,
    alpha = 0.18,
    size = 0.8,
    color = "grey35"
  ) +
  geom_ribbon(
    data = panel_a_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#476F67",
    alpha = 0.22
  ) +
  geom_line(
    data = panel_a_pred,
    aes(x = x, y = predicted),
    color = "#1F3A3D",
    linewidth = 0.95
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Incubation provisioning rate (standardized log mean)",
    y = "Predicted probability of survival to day 12"
  ) +
  prediction_theme

panel_b_pred <- as.data.frame(ggpredict(
  survival_model,
  terms = "nestling_temperature_sc [n=100]",
  type = "fixed"
))

panel_b <- ggplot() +
  geom_jitter(
    data = dat %>% filter(!is.na(nestling_temperature_sc)),
    aes(x = nestling_temperature_sc, y = D12_SURV),
    width = 0,
    height = 0.045,
    alpha = 0.18,
    size = 0.8,
    color = "grey35"
  ) +
  geom_ribbon(
    data = panel_b_pred,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    fill = "#6E93A2",
    alpha = 0.22
  ) +
  geom_line(
    data = panel_b_pred,
    aes(x = x, y = predicted),
    color = "#264B5A",
    linewidth = 0.95
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Nestling provisioning rate (standardized log mean)",
    y = "Predicted probability of survival to day 12"
  ) +
  prediction_theme

inc_pred <- as.data.frame(ggpredict(survival_model, terms = "EXP.INC", type = "fixed")) %>%
  mutate(panel = "Incubation treatment")
nest_pred <- as.data.frame(ggpredict(survival_model, terms = "EXP.NEST", type = "fixed")) %>%
  mutate(panel = "Nestling treatment")
panel_c_pred <- bind_rows(inc_pred, nest_pred) %>%
  mutate(panel = factor(panel, levels = c("Incubation treatment", "Nestling treatment")))

inc_raw <- dat %>%
  group_by(EXP.INC) %>%
  summarise(observed = mean(D12_SURV), .groups = "drop") %>%
  transmute(x = as.character(EXP.INC), observed, panel = "Incubation treatment")
nest_raw <- dat %>%
  group_by(EXP.NEST) %>%
  summarise(observed = mean(D12_SURV), .groups = "drop") %>%
  transmute(x = as.character(EXP.NEST), observed, panel = "Nestling treatment")
panel_c_raw <- bind_rows(inc_raw, nest_raw) %>%
  mutate(panel = factor(panel, levels = c("Incubation treatment", "Nestling treatment")))

panel_c <- ggplot(panel_c_pred, aes(x = x, y = predicted)) +
  geom_point(
    data = panel_c_raw,
    aes(x = x, y = observed),
    inherit.aes = FALSE,
    color = "grey35",
    alpha = 0.28,
    size = 2.4
  ) +
  geom_pointrange(
    aes(ymin = conf.low, ymax = conf.high),
    color = "#1F3A3D",
    linewidth = 0.55,
    size = 0.65
  ) +
  facet_wrap(~ panel, nrow = 1) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(
    x = "Experimental treatment",
    y = "Predicted probability of survival to day 12"
  ) +
  prediction_theme

survival_publication_figure <- (panel_a | panel_b) / panel_c +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag.position = c(0.02, 0.98))

pdf_file <- file.path("figures", "publication", "D12_SURV_publication_predictions.pdf")
png_file <- file.path("figures", "publication", "D12_SURV_publication_predictions.png")

ggsave(pdf_file, survival_publication_figure, width = 9.2, height = 8.2, device = cairo_pdf)
ggsave(png_file, survival_publication_figure, width = 9.2, height = 8.2, dpi = 600, bg = "white")

message("Saved publication survival figure to:")
message(pdf_file)
message(png_file)