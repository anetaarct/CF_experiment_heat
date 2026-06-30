
library(readxl)      
library(dplyr)       
library(tidyr)       
library(glmmTMB)
library(car)
library(DHARMa)
library(performance)
library(MuMIn)
library(ggplot2)
library(ggeffects)
library(sjPlot)


CF <- read_excel("C:/Users/IMBIO4/Downloads/CF_exp_all_2026.xlsx", 
                 sheet = "CF_exp.nestlings_analysis")



str(CF)


DATA_CF <- CF %>%
  select(F_RING, LD, INC, CS, BS, HD, SEX, YEAR, 
         EXP.INC, EXP.NEST,
         D2_MASS, D8_MASS, D12_MASS, D12_RTARS, Early.g, Late.g, Growth,
         LOGMEAN_INCUBATION, LOGMEAN_NESTLING, D12_SURV, GROUP, BOX, TEMP.BOX, TEMP.BOX.ENV, BLOCK)

DATA_CF <- DATA_CF %>%
  mutate(across(c(LD, INC, CS, BS, HD, Growth, Early.g, Late.g, 
                  D2_MASS, D8_MASS, D12_RTARS, D12_MASS, 
                  LOGMEAN_INCUBATION, LOGMEAN_NESTLING, TEMP.BOX, TEMP.BOX.ENV), as.numeric))

DATA_CF <- DATA_CF %>%
  mutate(across(c(F_RING, SEX, YEAR, EXP.INC, EXP.NEST, D12_SURV, GROUP, BOX, BLOCK), as.factor))

str(DATA_CF)

DATA_CF$HD_sc <- scale (DATA_CF$HD)
DATA_CF$BS_sc <- scale (DATA_CF$BS)
DATA_CF$incubation_temperature_sc <- scale(DATA_CF$LOGMEAN_INCUBATION)
DATA_CF$nestling_temperature_sc <- scale(DATA_CF$LOGMEAN_NESTLING)

cor.test(DATA_CF$HD, DATA_CF$LOGMEAN_NESTLING)
cor.test(DATA_CF$HD, DATA_CF$LOGMEAN_INCUBATION)

DATA_CF <- DATA_CF %>%
  filter(SEX %in% c("F", "M")) %>% 
  mutate(SEX = droplevels(SEX))
brood_order <- DATA_CF %>%
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

DATA_CF <- DATA_CF %>%
  inner_join(retained_broods, by = c("F_RING", "YEAR", "BOX", "LD")) %>%
  select(-retained_brood) %>%
  droplevels()


table(DATA_CF$SEX)


################# D2 MASS                                 


CF_D2 <- DATA_CF %>%
  filter(!is.na(D2_MASS))

print(nrow(CF_D2))
print(nrow(DATA_CF))

ggplot(CF_D2, aes(x = EXP.INC, y = D2_MASS, fill = EXP.INC)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "D2 MASS") +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(CF_D2, aes(D2_MASS, fill = EXP.INC)) + geom_histogram() + facet_wrap(~EXP.INC)


## D2 MASS
## model without dispformula
## gaussian distribution

CF_d2mass1 <- glmmTMB(D2_MASS ~ EXP.INC + incubation_temperature_sc + SEX + BS_sc + (1|F_RING) + (1|YEAR), 
                     data = CF_D2)

summary(CF_d2mass1)

## D2 MASS
## model with dispformula
## gaussian distribution

CF_d2mass2 <- glmmTMB(D2_MASS ~ EXP.INC + incubation_temperature_sc + SEX + BS_sc + (1|F_RING) + (1|YEAR), 
                      dispformula = ~ EXP.INC,
                      data = CF_D2)

summary(CF_d2mass2)


## D2 MASS
## model without dispformula
## t-student distribution

CF_d2mass3 <- glmmTMB(D2_MASS ~ EXP.INC + incubation_temperature_sc + SEX + BS_sc + (1|F_RING) + (1|YEAR),
                      family = t_family(),
                      data = CF_D2)

summary(CF_d2mass3) 


## D2 MASS
## model with dispformula
## t-student distribution

CF_d2mass4 <- glmmTMB(D2_MASS ~ EXP.INC + incubation_temperature_sc + SEX + BS_sc + (1|F_RING) + (1|YEAR),
                      family = t_family(),
                      dispformula = ~ EXP.INC,
                      data = CF_D2)

summary(CF_d2mass4) 


performance:: compare_performance(CF_d2mass1, CF_d2mass2, CF_d2mass3, CF_d2mass4)

## best model CF_d2mass3

summary(CF_d2mass3) 
res_d2mass3 <- simulateResiduals(CF_d2mass3)
plot(res_d2mass3)
check_collinearity(CF_d2mass3)
testDispersion(res_d2mass3)
testOutliers(res_d2mass3)
testUniformity(res_CF_d2mass3)
confint(CF_d2mass3)


emm_d2 <- emmeans(CF_d2mass3, ~ EXP.INC)
summary(emm_d2)
pairs(emm_d2)



########################################
##########################################

##d8 mass

CF_D8 <- DATA_CF %>%
  filter(!is.na(D8_MASS))

print(nrow(CF_D8))
print(nrow(DATA_CF))


ggplot(CF_D8, aes(x = GROUP, y = D8_MASS, fill = GROUP)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "D8 MASS") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_D8, aes(D8_MASS, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)


summary(CF_D8$D8_MASS)
tapply(CF_D8$D8_MASS, CF_D8$GROUP, summary)


## D8 MASS
## gaussian distribution

CF_d8mass1 <- glmmTMB(D8_MASS ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     data = CF_D8)

summary(CF_d8mass1)

## D8 MASS
## gaussian distribution
## model wih dispformula

CF_d8mass2 <- glmmTMB(D8_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      dispformula = ~ EXP.INC + EXP.NEST,
                      data = CF_D8)

summary(CF_d8mass2)

## D8 MASS
## t-student distribution

CF_d8mass3 <- glmmTMB(D8_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      family = t_family(),
                      data = CF_D8)

summary(CF_d8mass3)

## D8 MASS
## t-student distribution
## model wih dispformula

CF_d8mass4 <- glmmTMB(D8_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      family = t_family(),
                      dispformula = ~ EXP.INC + EXP.NEST,
                      data = CF_D8)

summary(CF_d8mass4)

performance:: compare_performance(CF_d8mass1, CF_d8mass2, CF_d8mass3, CF_d8mass4)


## best model CF_d8mass3


summary(CF_d8mass3) 
res_d8mass3 <- simulateResiduals(CF_d8mass3)
plot(res_d8mass3)
check_collinearity(CF_d8mass3)
testDispersion(res_d8mass3)
testOutliers(res_d8mass3)
testUniformity(res_CF_d8mass3)
confint(CF_d8mass3)



emm_d8 <- emmeans(CF_d8mass3, ~ EXP.INC)
summary(emm_d8)
pairs(emm_d8)


######################################################################
######################################################################
## d12

CF_D12 <- DATA_CF %>%
  filter(!is.na(D12_MASS))

print(nrow(CF_D12))
print(nrow(DATA_CF))


ggplot(CF_D12, aes(x = GROUP, y = D12_MASS, fill = GROUP)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "D12 MASS") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_D12, aes(D12_MASS, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)

summary(CF_D12$D12_MASS)
tapply(CF_D12$D12_MASS, CF_D12$GROUP, summary)

sum(is.na(CF_D12$F_RING))


## D12 MASS
## gaussian distribution

CF_d12mass1 <- glmmTMB(D12_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      data = CF_D12)

summary(CF_d12mass1)

## D12 MASS
## gaussian distribution
## model wih dispformula

CF_d12mass2 <- glmmTMB(D12_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      dispformula = ~ EXP.INC + EXP.NEST,
                      data = CF_D12)

summary(CF_d12mass2)

## D12 MASS
## t-student distribution

CF_d12mass3 <- glmmTMB(D12_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      family = t_family(),
                      data = CF_D12)

summary(CF_d12mass3)

## D12 MASS
## t-student distribution
## model wih dispformula

CF_d12mass4 <- glmmTMB(D12_MASS ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      family = t_family(),
                      dispformula = ~ EXP.INC + EXP.NEST,
                      data = CF_D12)

summary(CF_d12mass4)

performance:: compare_performance(CF_d12mass1, CF_d12mass2, CF_d12mass3, CF_d12mass4)


## CF_d12mass1 best option

summary(CF_d12mass1) 
res_d12mass1 <- simulateResiduals(CF_d12mass1)
plot(res_d12mass1)
check_collinearity(CF_d12mass1)
testDispersion(res_d12mass1)
testOutliers(res_d12mass1)
testUniformity(res_CF_d12mass1)
confint(CF_d12mass1)



emm_d12 <- emmeans(CF_d12mass1, ~ EXP.INC)
summary(emm_d12)
pairs(emm_d12)



#################################################################
#################################################################
##early body mass gain

CF_Early.g <- DATA_CF %>%
  filter(!is.na(Early.g))

print(nrow(CF_Early.g))
print(nrow(DATA_CF))


ggplot(CF_Early.g, aes(x = GROUP, y = Early.g, fill = GROUP)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "Early g") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_Early.g, aes(Early.g, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)


summary(CF_Early.g$Early.g)
tapply(CF_Early.g$Early.g, CF_Early.g$GROUP, summary)


## EARLY GROWTH
## gaussian distribution

CF_early1 <- glmmTMB(Early.g ~ EXP.INC  + incubation_temperature_sc +  
                         EXP.NEST + nestling_temperature_sc + 
                         SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                       data = CF_Early.g)

summary(CF_early1)

## EARLY GROWTH
## gaussian distribution
## model wih dispformula

CF_early2 <- glmmTMB(Early.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     dispformula = ~ EXP.INC + EXP.NEST,
                     data = CF_Early.g)

summary(CF_early2)


## EARLY GROWTH
## t-student distribution

CF_early3 <- glmmTMB(Early.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     family = t_family(),
                     data = CF_Early.g)

summary(CF_early3)



## EARLY GROWTH
## t-student distribution
## model wih dispformula

CF_early4 <- glmmTMB(Early.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     family = t_family(),
                     dispformula = ~ EXP.INC + EXP.NEST,
                     data = CF_Early.g)

summary(CF_early4)


performance:: compare_performance(CF_early1, CF_early2, CF_early3, CF_early4)

## CF_early4 best model

summary(CF_early4) 
res_early4 <- simulateResiduals(CF_early4)
plot(res_early4)
check_collinearity(CF_early4)
testDispersion(CF_early4)
testOutliers(res_early4)
testUniformity(res_CF_early4)
confint(CF_early4)



emm_early <- emmeans(CF_early4, ~ EXP.INC)
summary(emm_early)
pairs(emm_early)

###########################################################################
###########################################################################
## late body mass gain

CF_Late.g <- DATA_CF %>%
  filter(!is.na(Late.g))

print(nrow(CF_Late.g))
print(nrow(DATA_CF))


ggplot(CF_Late.g, aes(x = GROUP, y = Late.g, fill = GROUP)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "Late g") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_Late.g, aes(Late.g, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)


## LATE GROWTH
## gaussian distribution

CF_late1 <- glmmTMB(Late.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     data = CF_Late.g)

summary(CF_late1)

## LATE GROWTH
## gaussian distribution
## model wih dispformula

CF_late2 <- glmmTMB(Late.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     dispformula = ~ EXP.INC + EXP.NEST,
                     data = CF_Late.g)

summary(CF_late2)


## LATE GROWTH
## t-student distribution

CF_late3 <- glmmTMB(Late.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     family = t_family(),
                     data = CF_Late.g)

summary(CF_late3)



## LATE GROWTH
## t-student distribution
## model wih dispformula

CF_late4 <- glmmTMB(Late.g ~ EXP.INC  + incubation_temperature_sc +  
                       EXP.NEST + nestling_temperature_sc + 
                       SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                     family = t_family(),
                     dispformula = ~ EXP.INC + EXP.NEST,
                     data = CF_Late.g)

summary(CF_late4)


performance:: compare_performance(CF_late1, CF_late2, CF_late3, CF_late4)

## best model CF_late4

summary(CF_late4) 
res_CF_late4 <- simulateResiduals(CF_late4)
plot(res_CF_late4)
check_collinearity(CF_late4)
testDispersion(res_CF_late4)
testOutliers(res_CF_late4)
testUniformity(res_CF_late4)
confint(CF_late4)


emm_late <- emmeans(CF_late4, ~ EXP.INC)
summary(emm_late)
pairs(emm_late)





#####################################################################
##################################################################
##entire growth


CF_g <- DATA_CF %>%
  filter(!is.na(Growth))

print(nrow(CF_g ))
print(nrow(DATA_CF))


ggplot(CF_g, aes(x = GROUP, y = Growth, fill = GROUP)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "growth") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_g, aes(Growth, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)


## GROWTH
## gaussian distribution

CF_growth1 <- glmmTMB(Growth ~ EXP.INC  + incubation_temperature_sc +  
                      EXP.NEST + nestling_temperature_sc + 
                      SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                    data = CF_g)

summary(CF_growth1)

## GROWTH
## gaussian distribution
## model wih dispformula

CF_growth2 <- glmmTMB(Growth ~ EXP.INC  + incubation_temperature_sc +  
                      EXP.NEST + nestling_temperature_sc + 
                      SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                    dispformula = ~ EXP.INC + EXP.NEST,
                    data = CF_g)

summary(CF_growth2)


## GROWTH
## t-student distribution

CF_growth3 <- glmmTMB(Growth ~ EXP.INC  + incubation_temperature_sc +  
                      EXP.NEST + nestling_temperature_sc + 
                      SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                    family = t_family(),
                    data = CF_g)

summary(CF_growth3)



## GROWTH
## t-student distribution
## model wih dispformula

CF_growth4 <- glmmTMB(Growth ~ EXP.INC  + incubation_temperature_sc +  
                      EXP.NEST + nestling_temperature_sc + 
                      SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                    family = t_family(),
                    dispformula = ~ EXP.INC + EXP.NEST,
                    data = CF_g)

summary(CF_growth4)


performance:: compare_performance(CF_growth1, CF_growth2, CF_growth3, CF_growth4)


## CF_growth3 best model

summary(CF_growth3) 
res_CF_growth3 <- simulateResiduals(CF_growth3)
plot(res_CF_growth3)
check_collinearity(CF_growth3)
testDispersion(res_CF_growth3)
testOutliers(res_CF_growth3)
testUniformity(res_CF_growth3)
confint(CF_growth3)


emm_growth <- emmeans(CF_growth3, ~ EXP.INC)
summary(emm_growth)
pairs(emm_growth)



#################################################################
################################################################
## SCALED MASS INDEX


CF_smi <- DATA_CF %>%
  drop_na(D12_MASS, D12_RTARS)

sum(is.na(CF_smi$D12_MASS))
sum(is.na(CF_smi$D12_RTARS))

CF_smi$log_mass_n <- log(CF_smi$D12_MASS)
CF_smi$log_tars_n <- log(CF_smi$D12_RTARS)


sd_log_mass_n <- sd(CF_smi$log_mass_n, na.rm = TRUE)
sd_log_tars_n <- sd(CF_smi$log_tars_n, na.rm = TRUE)
b_sma_chicks  <- sd_log_mass_n / sd_log_tars_n

L0_chicks <- mean(CF_smi$D12_RTARS, na.rm = TRUE)

CF_smi$D12_SMI <- CF_smi$D12_MASS * (L0_chicks / CF_smi$D12_RTARS)^b_sma_chicks


ggplot(CF_smi, aes(x = GROUP, y = D12_SMI, fill = GROUP)) +
  geom_boxplot() +
  labs(title = "SMI NESTLINGS",
       x = "Group",
       y = "SMI") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_smi, aes(D12_SMI, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)

## SMI
## gaussian distribution

CF_SMI1 <- glmmTMB(D12_SMI ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      data = CF_smi)

summary(CF_SMI1)

## SMI
## gaussian distribution
## model wih dispformula

CF_SMI2 <- glmmTMB(D12_SMI ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      dispformula = ~ EXP.INC + EXP.NEST,
                      data = CF_smi)

summary(CF_SMI2)


## SMI
## t-student distribution

CF_SMI3 <- glmmTMB(D12_SMI ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      family = t_family(),
                      data = CF_smi)

summary(CF_SMI3)



## SMI
## t-student distribution
## model wih dispformula

CF_SMI4 <- glmmTMB(D12_SMI ~ EXP.INC  + incubation_temperature_sc +  
                        EXP.NEST + nestling_temperature_sc + 
                        SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                      family = t_family(),
                      dispformula = ~ EXP.INC + EXP.NEST,
                      data = CF_smi)

summary(CF_SMI4)


performance:: compare_performance(CF_SMI1, CF_SMI2, CF_SMI3, CF_SMI4)

## CF_SMI3 best model

summary(CF_SMI3) 
res_CF_SMI3 <- simulateResiduals(CF_SMI3)
plot(res_CF_SMI3)
check_collinearity(CF_SMI3)
testDispersion(res_CF_SMI3)
testOutliers(res_CF_SMI3)
testUniformity(res_CF_SMI3)
confint(CF_SMI3)


emm_SMI <- emmeans(CF_SMI3, ~ EXP.INC)
summary(emm_SMI)
pairs(emm_SMI)

##############################################################################
###########################################################################
## TARSUS

CF_TARS <- DATA_CF %>%
  filter(!is.na(D12_RTARS))

print(nrow(CF_TARS))
print(nrow(DATA_CF))


ggplot(CF_TARS, aes(x = GROUP, y = D12_RTARS, fill = GROUP)) +
  geom_boxplot() +
  labs(x = "GROUP",
       y = "TARSUS") +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(CF_TARS, aes(D12_RTARS, fill = GROUP)) + geom_histogram() + facet_wrap(~GROUP)



## tarsus
## gaussian distribution

CF_tars1 <- glmmTMB(D12_RTARS ~ EXP.INC  + incubation_temperature_sc +  
                     EXP.NEST + nestling_temperature_sc + 
                     SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                   data = CF_TARS )

summary(CF_tars1)

## tarsus
## gaussian distribution
## model wih dispformula

CF_tars2 <- glmmTMB(D12_RTARS ~ EXP.INC  + incubation_temperature_sc +  
                     EXP.NEST + nestling_temperature_sc + 
                     SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                   dispformula = ~ EXP.INC + EXP.NEST,
                   data = CF_TARS)

summary(CF_tars2)


## tarsus
## t-student distribution

CF_tars3 <- glmmTMB(D12_RTARS ~ EXP.INC  + incubation_temperature_sc +  
                     EXP.NEST + nestling_temperature_sc + 
                     SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                   family = t_family(),
                   data = CF_TARS)

summary(CF_tars3)



## tarsus
## t-student distribution
## model wih dispformula

CF_tars4 <- glmmTMB(D12_RTARS ~ EXP.INC  + incubation_temperature_sc +  
                     EXP.NEST + nestling_temperature_sc + 
                     SEX + BS_sc  +(1|F_RING) + (1|YEAR), 
                   family = t_family(),
                   dispformula = ~ EXP.INC + EXP.NEST,
                   data = CF_TARS)

summary(CF_tars4)


performance:: compare_performance(CF_tars1, CF_tars2, CF_tars3, CF_tars4)


## CF_tars3  best model

summary(CF_tars3) 
res_CF_tars3 <- simulateResiduals(CF_tars3)
plot(res_CF_tars3)
check_collinearity(CF_tars3)
testDispersion(res_CF_tars3)
testOutliers(res_CF_tars3)
testUniformity(res_CF_tars3)
confint(CF_tars3)


emm_tars <- emmeans(CF_tars3, ~ EXP.INC)
summary(emm_tars)
pairs(emm_tars)

## SURVIVAL

## D12 SURVIVAL

DATA_CF$D12_SURV <- factor(DATA_CF$D12_SURV, 
                           levels = c(0, 1), 
                           labels = c("Dead", "Survived"))



model_surv <- glmmTMB(D12_SURV ~ EXP.INC + EXP.NEST+ incubation_temperature_sc 
                       + nestling_temperature_sc +SEX+ + BS_sc 
                       + (1|F_RING) + (1|YEAR),
                       family = binomial(),
                       data = DATA_CF)

summary(model_surv)

res_surv <- simulateResiduals(model_surv)
plot(res_surv)

testUniformity(res_surv)
testDispersion(res_surv)
testOutliers(res_surv)

check_collinearity(model_surv)


ci_m_surv <- confint(model_surv)
ci_m_surv

emm_surv <- emmeans(model_surv, ~ EXP.INC, type = "response")
summary(emm_surv)

pairs(emm_surv)

sessionInfo()

