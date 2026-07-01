suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

dir.create("models", showWarnings = FALSE)
dir.create(file.path("models", "manipulation_check"), showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE)
dir.create(file.path("figures", "manipulation_check"), showWarnings = FALSE, recursive = TRUE)

broods <- read_excel(
  "CF_exp_all_2026.xlsx",
  sheet = "CF_exp.broods_analysis",
  na = c("NA", "", "na", "NaN")
) %>%
  mutate(
    TEMP.BOX = suppressWarnings(as.numeric(TEMP.BOX)),
    EXP.INC = factor(EXP.INC, levels = c("Control", "Heated")),
    YEAR = factor(YEAR)
  ) %>%
  filter(!is.na(TEMP.BOX), !is.na(EXP.INC), !is.na(YEAR)) %>%
  droplevels()

manipulation_model <- lm(TEMP.BOX ~ EXP.INC + YEAR, data = broods)

means <- broods %>%
  group_by(EXP.INC) %>%
  summarise(
    n = n(),
    mean_temp = mean(TEMP.BOX),
    se_temp = sd(TEMP.BOX) / sqrt(n()),
    .groups = "drop"
  )

drop_tests <- drop1(manipulation_model, test = "F")
test_table <- as.data.frame(drop_tests) %>%
  tibble::rownames_to_column("term") %>%
  filter(term %in% c("EXP.INC", "YEAR")) %>%
  transmute(
    term,
    df = Df,
    sum_of_squares = `Sum of Sq`,
    F_value = `F value`,
    p_value = `Pr(>F)`
  )

treatment_estimate <- tibble(
  term = "EXP.INCHeated",
  estimate = coef(summary(manipulation_model))["EXP.INCHeated", "Estimate"],
  std_error = coef(summary(manipulation_model))["EXP.INCHeated", "Std. Error"],
  statistic = coef(summary(manipulation_model))["EXP.INCHeated", "t value"],
  p_value = coef(summary(manipulation_model))["EXP.INCHeated", "Pr(>|t|)"],
  residual_df = df.residual(manipulation_model),
  n = nrow(broods)
)

write_csv(means, "models/manipulation_check/nestbox_temperature_means.csv")
write_csv(test_table, "models/manipulation_check/nestbox_temperature_tests.csv")
write_csv(treatment_estimate, "models/manipulation_check/nestbox_temperature_treatment_estimate.csv")
capture.output(
  summary(manipulation_model),
  file = "models/manipulation_check/nestbox_temperature_lm_summary.txt"
)

pred <- means %>%
  mutate(
    ymin = mean_temp - se_temp,
    ymax = mean_temp + se_temp
  )

p <- ggplot() +
  geom_jitter(
    data = broods,
    aes(x = EXP.INC, y = TEMP.BOX),
    width = 0.12,
    height = 0,
    alpha = 0.18,
    size = 1.2,
    color = "grey35"
  ) +
  geom_pointrange(
    data = pred,
    aes(x = EXP.INC, y = mean_temp, ymin = ymin, ymax = ymax),
    color = "#1F3A3D",
    linewidth = 0.55,
    size = 0.72
  ) +
  labs(
    x = "Incubation-period treatment",
    y = "Nest-box temperature (°C)"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("figures/manipulation_check/nestbox_temperature_manipulation.png", p, width = 6.2, height = 4.2, dpi = 300)
ggsave("figures/manipulation_check/nestbox_temperature_manipulation.pdf", p, width = 6.2, height = 4.2)

message("Temperature manipulation check saved in models/manipulation_check and figures/manipulation_check.")
