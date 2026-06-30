suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

source_file <- normalizePath("CF_exp_all_2026.xlsx", winslash = "\\", mustWork = TRUE)
nestling_sheet <- "CF_exp.nestlings_analysis"
source_file_for_readxl <- file.path(tempdir(), basename(source_file))
if (file.exists(source_file_for_readxl)) {
  unlink(source_file_for_readxl)
}
ps_quote <- function(x) paste0("'", gsub("'", "''", x, fixed = TRUE), "'")
copy_command <- paste(
  "Copy-Item -LiteralPath", ps_quote(source_file),
  "-Destination", ps_quote(source_file_for_readxl),
  "-Force"
)
system2("powershell", c("-NoProfile", "-Command", copy_command))
if (!file.exists(source_file_for_readxl) || file.info(source_file_for_readxl)$size == 0) {
  stop("Could not copy the Excel workbook to a temporary file for readxl.")
}

dir.create("data", showWarnings = FALSE)

raw_nestlings <- read_excel(
  source_file_for_readxl,
  sheet = nestling_sheet,
  na = c("NA", "", "na", "NaN")
)

analysis_vars <- c(
  "F_RING", "LD", "INC", "CS", "BS", "HD", "SEX", "YEAR",
  "EXP.INC", "EXP.NEST",
  "D2_MASS", "D8_MASS", "D12_MASS", "D12_RTARS", "Early.g", "Late.g", "Growth",
  "LOGMEAN_INCUBATION", "LOGMEAN_NESTLING", "D12_SURV", "GROUP", "BOX", "TEMP.BOX",
  "TEMP.BOX.ENV", "BLOCK"
)

missing_vars <- setdiff(analysis_vars, names(raw_nestlings))
if (length(missing_vars) > 0) {
  stop("Missing expected columns: ", paste(missing_vars, collapse = ", "))
}

numeric_vars <- c(
  "LD", "INC", "CS", "BS", "HD", "Growth", "Early.g", "Late.g",
  "D2_MASS", "D8_MASS", "D12_RTARS", "D12_MASS",
  "LOGMEAN_INCUBATION", "LOGMEAN_NESTLING", "TEMP.BOX", "TEMP.BOX.ENV"
)

factor_vars <- c("F_RING", "SEX", "YEAR", "EXP.INC", "EXP.NEST", "D12_SURV", "GROUP", "BOX", "BLOCK")

nestlings_typed <- raw_nestlings %>%
  select(all_of(analysis_vars)) %>%
  mutate(
    across(any_of(c("SEX", "EXP.INC", "EXP.NEST", "GROUP")), ~ str_squish(as.character(.x))),
    across(any_of(numeric_vars), ~ suppressWarnings(as.numeric(.x))),
    across(any_of(factor_vars), as.factor)
  ) %>%
  filter(SEX %in% c("F", "M")) %>%
  mutate(SEX = droplevels(SEX))

brood_order <- nestlings_typed %>%
  distinct(F_RING, YEAR, BOX, LD, .keep_all = FALSE) %>%
  arrange(F_RING, YEAR, LD, BOX) %>%
  group_by(F_RING) %>%
  mutate(female_brood_order = row_number()) %>%
  ungroup()

repeated_females <- brood_order %>%
  count(F_RING, name = "n_broods") %>%
  filter(n_broods > 1)

retained_broods <- brood_order %>%
  filter(female_brood_order == 1) %>%
  transmute(F_RING, YEAR, BOX, LD, retained_brood = TRUE)

prepared_nestlings <- nestlings_typed %>%
  inner_join(retained_broods, by = c("F_RING", "YEAR", "BOX", "LD")) %>%
  select(-retained_brood) %>%
  droplevels() %>%
  mutate(
    HD_sc = as.numeric(scale(HD)),
    BS_sc = as.numeric(scale(BS)),
    incubation_temperature_sc = as.numeric(scale(LOGMEAN_INCUBATION)),
    nestling_temperature_sc = as.numeric(scale(LOGMEAN_NESTLING))
  )

smi_base <- prepared_nestlings %>%
  drop_na(D12_MASS, D12_RTARS) %>%
  filter(D12_MASS > 0, D12_RTARS > 0) %>%
  mutate(
    log_mass_n = log(D12_MASS),
    log_tars_n = log(D12_RTARS)
  )

b_sma_chicks <- sd(smi_base$log_mass_n, na.rm = TRUE) / sd(smi_base$log_tars_n, na.rm = TRUE)
L0_chicks <- mean(smi_base$D12_RTARS, na.rm = TRUE)

prepared_nestlings <- prepared_nestlings %>%
  mutate(
    D12_SMI = if_else(
      !is.na(D12_MASS) & !is.na(D12_RTARS) & D12_MASS > 0 & D12_RTARS > 0,
      D12_MASS * (L0_chicks / D12_RTARS)^b_sma_chicks,
      NA_real_
    )
  )

prepared_data_checks <- tibble(
  item = c(
    "Raw nestling records",
    "Records after sex filter",
    "Records after repeated-female brood rule",
    "Study years",
    "Unique females",
    "Unique broods",
    "Females with repeated broods removed",
    "D2 mass available",
    "D8 mass available",
    "D12 mass available",
    "D12 tarsus available",
    "Early growth available",
    "Late growth available",
    "Whole-period growth available",
    "D12 SMI available",
    "D12 survival available"
  ),
  value = c(
    nrow(raw_nestlings),
    nrow(nestlings_typed),
    nrow(prepared_nestlings),
    paste(sort(unique(prepared_nestlings$YEAR)), collapse = ", "),
    n_distinct(prepared_nestlings$F_RING),
    n_distinct(paste(prepared_nestlings$F_RING, prepared_nestlings$YEAR, prepared_nestlings$BOX, prepared_nestlings$LD)),
    nrow(repeated_females),
    sum(!is.na(prepared_nestlings$D2_MASS)),
    sum(!is.na(prepared_nestlings$D8_MASS)),
    sum(!is.na(prepared_nestlings$D12_MASS)),
    sum(!is.na(prepared_nestlings$D12_RTARS)),
    sum(!is.na(prepared_nestlings$Early.g)),
    sum(!is.na(prepared_nestlings$Late.g)),
    sum(!is.na(prepared_nestlings$Growth)),
    sum(!is.na(prepared_nestlings$D12_SMI)),
    sum(!is.na(prepared_nestlings$D12_SURV))
  )
)

treatment_counts <- prepared_nestlings %>%
  count(EXP.INC, EXP.NEST, GROUP, name = "nestling_records") %>%
  left_join(
    prepared_nestlings %>%
      distinct(F_RING, YEAR, BOX, LD, EXP.INC, EXP.NEST, GROUP) %>%
      count(EXP.INC, EXP.NEST, GROUP, name = "nests"),
    by = c("EXP.INC", "EXP.NEST", "GROUP")
  ) %>%
  select(EXP.INC, EXP.NEST, GROUP, nests, nestling_records)

analysis_columns <- c(
  "F_RING", "YEAR", "SEX",
  "EXP.INC", "EXP.NEST", "GROUP",
  "D2_MASS", "D8_MASS", "D12_MASS", "D12_RTARS",
  "Early.g", "Late.g", "Growth", "D12_SMI", "D12_SURV",
  "BS_sc", "incubation_temperature_sc", "nestling_temperature_sc"
)

analysis_nestlings <- prepared_nestlings %>%
  select(all_of(analysis_columns))

saveRDS(analysis_nestlings, "data/prepared_nestlings.rds")
write_csv(analysis_nestlings, "data/prepared_nestlings.csv")
saveRDS(repeated_females, "data/repeated_females.rds")
saveRDS(list(b_sma_chicks = b_sma_chicks, L0_chicks = L0_chicks), "data/smi_parameters.rds")
write_csv(prepared_data_checks, "data/prepared_data_checks.csv")
write_csv(treatment_counts, "data/treatment_counts.csv")
write_csv(repeated_females, "data/repeated_females.csv")

message("Prepared nestling data saved to data/prepared_nestlings.rds and data/prepared_nestlings.csv")
