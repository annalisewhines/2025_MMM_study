library(rstanarm)
library(bayestestR)
library(dplyr)
library(lubridate)
library(tidyr)
library(here)

# Set working directory
here::i_am("2025_MMM_study/code/analyses/task_analyses.R")
setwd(here::here())

# Source functions
source(here("2025_MMM_study", "code", "functions", "helper_fun.R"))
source(here("2025_MMM_study" ,"code", "functions", "plot_funs.R"))

# Read data
m3_indiv_params_men <- readRDS(here("2025_MMM_study", "data", "model_fits", "menstrual", "m3_indiv_params_men.RDS"))
m3_indiv_params_ovu <- readRDS(here("2025_MMM_study", "data", "model_fits", "ovulation", "m3_indiv_params_ovu.RDS"))
main_mhq_data <- readRDS(here("2025_MMM_study", "data", "processed_dat", "main_analysis", "main_mhq_data.RDS"))

# Analysis settings
sens_cycle <- FALSE
sens_time  <- FALSE

# Combine individual parameters from both groups
combined_params <- rbind(
  data.frame(m3_indiv_params_men$individual, phase = "perimen"),
  data.frame(m3_indiv_params_ovu$individual, phase = "ovulation")
)

# Filter to selected analysis IDs if sensitivity analysis
if (sens_cycle) {
  combined_params <- combined_params %>% filter(subj_id %in% analysis_id_cycle$subj_id)
} else if (sens_time) {
  combined_params <- combined_params %>% filter(subj_id %in% analysis_id_time$subj_id)
}

# Merge with demographics to add age and group (A or B)
combined_params <- combined_params %>%
  left_join(main_mhq_data$screen_demo %>% select(subj_id, age, group),
            by = "subj_id")

# Reorder columns for clarity
combined_params <- combined_params %>%
  select(subj_id, phase, group, age, kE, kR, a)

# Check that age is numeric
combined_params$age <- as.numeric(combined_params$age)

# Participant characteristic table info
task_descrip <- bind_rows(
  combined_params %>%
    group_by(phase) %>%
    summarise(
      group = "All",
      kE_mean = round(mean(kE, na.rm = TRUE), 2),
      kE_sd = round(sd(kE, na.rm = TRUE), 2),
      kE_min = round(min(kE, na.rm = TRUE), 2),
      kE_max = round(max(kE, na.rm = TRUE), 2),
      kR_mean = round(mean(kR, na.rm = TRUE), 2),
      kR_sd = round(sd(kR, na.rm = TRUE), 2),
      kR_min = round(min(kR, na.rm = TRUE), 2),
      kR_max = round(max(kR, na.rm = TRUE), 2),
      a_mean = round(mean(a, na.rm = TRUE), 2),
      a_sd = round(sd(a, na.rm = TRUE), 2),
      a_min = round(min(a, na.rm = TRUE), 2),
      a_max = round(max(a, na.rm = TRUE), 2),
      n = n()
    ),
  combined_params %>%
    group_by(group, phase) %>%
    summarise(
      kE_mean = round(mean(kE, na.rm = TRUE), 2),
      kE_sd = round(sd(kE, na.rm = TRUE), 2),
      kE_min = round(min(kE, na.rm = TRUE), 2),
      kE_max = round(max(kE, na.rm = TRUE), 2),
      kR_mean = round(mean(kR, na.rm = TRUE), 2),
      kR_sd = round(sd(kR, na.rm = TRUE), 2),
      kR_min = round(min(kR, na.rm = TRUE), 2),
      kR_max = round(max(kR, na.rm = TRUE), 2),
      a_mean = round(mean(a, na.rm = TRUE), 2),
      a_sd = round(sd(a, na.rm = TRUE), 2),
      a_min = round(min(a, na.rm = TRUE), 2),
      a_max = round(max(a, na.rm = TRUE), 2),
      n = n(),
      .groups = "drop"
    )
)

# Rescale the parameters and change scores between 0 and 1
combined_params <- combined_params %>%
  mutate(
    kE_scaled        = rescale(kE),
    kR_scaled        = rescale(kR),
    a_scaled         = rescale(a)
  )

# Effort sensitivity
kE_mlm <- stan_glmer(
  kE_scaled ~ phase + age + group + (1 | subj_id),
  data = combined_params,
  iter = 2000,
  seed = 123
)
effort_mean <- fixef(kE_mlm)
effort_mean
print(bayestestR::hdi(kE_mlm), digits = 3)

# Reward sensitivity
kR_mlm <- stan_glmer(
  kR_scaled ~ phase + age + group + (1 | subj_id),
  data = combined_params,
  iter = 2000,
  seed = 123
)
reward_mean <- fixef(kR_mlm)
reward_mean
print(bayestestR::hdi(kR_mlm), digits = 3)

# Acceptance bias
a_mlm <- stan_glmer(
  a_scaled ~ phase + age + group + (1 | subj_id),
  data = combined_params,
  iter = 2000,
  seed = 123
)
a_mean <- fixef(a_mlm)
a_mean
print(bayestestR::hdi(a_mlm), digits = 3)

# Acceptance Bias plot
a_plot <- raincloud_plot(
  dat = combined_params %>%
    mutate(phase = case_when(
      phase == "ovulation" ~ "ovu",
      phase == "perimen"   ~ "men"
    )) %>%
    mutate(phase = factor(phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Acceptance Bias",
  predictor_var = "phase",
  outcome_var = "a",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("acceptance_bias_raincloud.png", plot = a_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

kR_plot <- raincloud_plot(
  dat = combined_params %>%
    mutate(phase = case_when(
      phase == "ovulation" ~ "ovu",
      phase == "perimen"   ~ "men"
    )) %>%
    mutate(phase = factor(phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Reward Sensitivity",
  predictor_var = "phase",
  outcome_var = "kR",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("reward_sensitivity_raincloud.png", plot = kR_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

# Effort Sensitivity plot
kE_plot <- raincloud_plot(
  dat = combined_params %>%
    mutate(phase = case_when(
      phase == "ovulation" ~ "ovu",
      phase == "perimen"   ~ "men"
    )) %>%
    mutate(phase = factor(phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Effort Sensitivity",
  predictor_var = "phase",
  outcome_var = "kE",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("effort_sensitivity_raincloud.png", plot = kE_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

# Save data for change analysis
saveRDS(combined_params, here("2025_MMM_study", "data", "processed_dat", "main_analysis", "combined_params.RDS"))
