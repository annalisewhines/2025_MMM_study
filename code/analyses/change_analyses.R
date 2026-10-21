library(rstanarm)
library(bayestestR)
library(dplyr)
library(lubridate)
library(tidyr)
library(here)

# Set working directory
here::i_am("2025_MMM_study/code/analyses/change_analyses.R")

# Load data
mhq_analysis    <- readRDS(here("2025_MMM_study", "data", "processed_dat", "main_analysis", "mhq_analysis.RDS"))
combined_params <- readRDS(here("2025_MMM_study", "data", "processed_dat", "main_analysis", "combined_params.RDS"))

## Data preparation

# Create mhq_analysis_wide
mhq_analysis_wide <- mhq_analysis %>%
  select(subj_id, cycle_phase, measure, total) %>%
  pivot_wider(names_from = c(cycle_phase, measure), 
              values_from = total,
              names_sep = "_")

# Mental health residualised changes
# PHQ
phq_resid_model <- lm(men_phq ~ ovu_phq, data = mhq_analysis_wide)
mhq_analysis_wide$phq_resid_change <- residuals(phq_resid_model)

# SHAPS
shaps_resid_model <- lm(men_shaps ~ ovu_shaps, data = mhq_analysis_wide)
mhq_analysis_wide$shaps_resid_change <- residuals(shaps_resid_model)

# STAI
stai_resid_model <- lm(men_stai ~ ovu_stai, data = mhq_analysis_wide)
mhq_analysis_wide$stai_resid_change <- residuals(stai_resid_model)

# PMS
pms_resid_model <- lm(men_pms ~ ovu_pms, data = mhq_analysis_wide)
mhq_analysis_wide$pms_resid_change <- residuals(pms_resid_model)

# Create params_wide
params_wide <- combined_params %>%
  select(subj_id, phase, kR, a, kE) %>%
  pivot_wider(names_from = phase, values_from = c(kR, a, kE))

# Create residualised parameter changes
kR_resid_model <- lm(kR_perimen ~ kR_ovulation, data = params_wide)
params_wide$kR_resid_change <- residuals(kR_resid_model)

a_resid_model <- lm(a_perimen ~ a_ovulation, data = params_wide)
params_wide$a_resid_change <- residuals(a_resid_model)

kE_resid_model <- lm(kE_perimen ~ kE_ovulation, data = params_wide)
params_wide$kE_resid_change <- residuals(kE_resid_model)

# Combine residualised changes with age and group
change_analysis <- params_wide %>%
  select(subj_id, kR_resid_change, a_resid_change, kE_resid_change) %>%
  left_join(mhq_analysis_wide %>% 
              select(subj_id, phq_resid_change, shaps_resid_change, 
                     stai_resid_change, pms_resid_change), 
            by = "subj_id") %>%
  left_join(combined_params %>% 
              filter(phase == "perimen") %>%
              select(subj_id, age, group), 
            by = "subj_id")

change_analysis %>%
  summarise(
    phq_sd = sd(phq_resid_change, na.rm = TRUE),
    shaps_sd = sd(shaps_resid_change, na.rm = TRUE),
    stai_sd = sd(stai_resid_change, na.rm = TRUE),
    pms_sd = sd(pms_resid_change, na.rm = TRUE),
    kR_sd = sd(kR_resid_change, na.rm = TRUE),
    a_sd = sd(a_resid_change, na.rm = TRUE),
    kE_sd = sd(kE_resid_change, na.rm = TRUE)
  )

## kR predicting mental health

# PHQ
kR_phq <- stan_glm(
  phq_resid_change ~ kR_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(kR_phq))
print(bayestestR::hdi(kR_phq), digits = 3)

# SHAPS
kR_shaps <- stan_glm(
  shaps_resid_change ~ kR_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(kR_shaps))
print(bayestestR::hdi(kR_shaps), digits = 3)

# STAI
kR_stai <- stan_glm(
  stai_resid_change ~ kR_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(kR_stai))
print(bayestestR::hdi(kR_stai), digits = 3)

# PMS
kR_pms <- stan_glm(
  pms_resid_change ~ kR_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(kR_pms))
print(bayestestR::hdi(kR_pms), digits = 3)


## a predicting mental health

# PHQ
a_phq <- stan_glm(
  phq_resid_change ~ a_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(a_phq))
print(bayestestR::hdi(a_phq), digits = 3)

# SHAPS
a_shaps <- stan_glm(
  shaps_resid_change ~ a_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(a_shaps))
print(bayestestR::hdi(a_shaps), digits = 3)

# STAI
a_stai <- stan_glm(
  stai_resid_change ~ a_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(a_stai))
print(bayestestR::hdi(a_stai), digits = 3)

# PMS
a_pms <- stan_glm(
  pms_resid_change ~ a_resid_change + age + group,
  data = change_analysis,
  iter = 2000, seed = 123
)
print(fixef(a_pms))
print(bayestestR::hdi(a_pms), digits = 3)
