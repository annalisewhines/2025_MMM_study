library(rstanarm)
library(bayestestR)
library(dplyr)
library(lubridate)
library(tidyr)
library(here)
library(patchwork)

# Set working directory
here::i_am("2025_MMM_study/code/analyses/mhq_analyses.R")

source(here("2025_MMM_study", "code", "functions", "helper_fun.R"))
source(here("2025_MMM_study", "code", "functions", "plot_funs.R"))

# Load data
main_mhq_data <- readRDS(here("2025_MMM_study", "data", "processed_dat", "main_analysis", "main_mhq_data.RDS"))

# Sensitivity analysis data
# mhq_data_cycle <- readRDS(here("2025_MMM_study", "data", "processed_dat", "sens_cycle", "mhq_data_cycle.RDS"))
# mhq_data_time  <- readRDS(here("2025_MMM_study", "data", "processed_dat", "sens_time", "mhq_data_time.RDS"))

#  Sensitivity analysis settings
sens_cycle <- FALSE
sens_time  <- FALSE

# Select data based on sensitivity analysis settings
if (sens_cycle) {
  selected_mhq_data <- mhq_data_cycle
} else if (sens_time) {
  selected_mhq_data <- mhq_data_time
} else {
  selected_mhq_data <- main_mhq_data
}

# Merge all the questionnaire data
mhq_analysis <- selected_mhq_data$ovu_shaps %>%
  select(subj_id, total_shaps) %>%
  rename(ovu_shaps = total_shaps) %>%
  left_join(selected_mhq_data$men_shaps %>% select(subj_id, total_shaps) %>% rename(men_shaps = total_shaps), by = "subj_id") %>%
  left_join(selected_mhq_data$ovu_phq %>% select(subj_id, total_phq) %>% rename(ovu_phq = total_phq), by = "subj_id") %>%
  left_join(selected_mhq_data$men_phq %>% select(subj_id, total_phq) %>% rename(men_phq = total_phq), by = "subj_id") %>%
  left_join(selected_mhq_data$ovu_stai %>% select(subj_id, total_stai) %>% rename(ovu_stai = total_stai), by = "subj_id") %>%
  left_join(selected_mhq_data$men_stai %>% select(subj_id, total_stai) %>% rename(men_stai = total_stai), by = "subj_id") %>%
  left_join(selected_mhq_data$ovu_pms %>% select(subj_id, total_pms) %>% rename(ovu_pms = total_pms), by = "subj_id") %>%
  left_join(selected_mhq_data$men_pms %>% select(subj_id, total_pms) %>% rename(men_pms = total_pms), by = "subj_id") %>%
  
# Reshape to long format
  pivot_longer(
    cols = -subj_id,
    names_to = c("cycle_phase", "measure"),
    names_sep = "_",
    values_to = "total"
  ) %>%
  
  # Add group and age
  left_join(selected_mhq_data$screen_demo %>% select(subj_id, age, group), by = "subj_id") %>%
  mutate(age = as.numeric(age)) %>%
  
# Scale total scores to 0-1
  group_by(measure) %>%
  mutate(total_scaled = rescale(total)) %>%
  ungroup()

# Set factor levels
mhq_analysis <- mhq_analysis %>%
  mutate(cycle_phase = factor(cycle_phase, levels = c("ovu", "men")))

# Descriptive statistics
descrip <- bind_rows(
  mhq_analysis %>%
    group_by(measure, cycle_phase) %>%
    summarise(
      group = "All",
      mean = mean(total),
      sd = sd(total),
      range = paste(min(total), "-", max(total)),
      .groups = "drop"
    ),
  mhq_analysis %>%
    group_by(measure, cycle_phase, group) %>%
    summarise(
      mean = mean(total),
      sd = sd(total),
      range = paste(min(total), "-", max(total)),
      .groups = "drop"
    )
)

# SHAPS model
shaps_mlm <- stan_glmer(
  total_scaled ~ cycle_phase + age + group + (1 | subj_id),
  data = mhq_analysis %>% filter(measure == "shaps"),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

shaps_mean <- fixef(shaps_mlm)
shaps_mean
print(bayestestR::hdi(shaps_mlm), digits = 3)


# PHQ model
phq_mlm <- stan_glmer(
  total_scaled ~ cycle_phase + age + group + (1 | subj_id),
  data = mhq_analysis %>% filter(measure == "phq"),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

phq_mean <- fixef(phq_mlm)
phq_mean
print(bayestestR::hdi(phq_mlm), digits = 3)

# STAI model
stai_mlm <- stan_glmer(
  total_scaled ~ cycle_phase + age + group + (1 | subj_id),
  data = mhq_analysis %>% filter(measure == "stai"),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

stai_mean <- fixef(stai_mlm)
stai_mean
print(bayestestR::hdi(stai_mlm), digits = 3)

# PMS model
pms_mlm <- stan_glmer(
  total_scaled ~ cycle_phase + age + group + (1 | subj_id),
  data = mhq_analysis %>% filter(measure == "pms"),
  iter = 4000,
  warmup = 2000,
  chains = 4,
  seed = 123
)

pms_mean <- fixef(pms_mlm)
pms_mean
print(bayestestR::hdi(pms_mlm), digits = 3)

# SHAPS plot
shaps_plot <- raincloud_plot(
  dat = mhq_analysis %>% 
    filter(measure == "shaps") %>%
    mutate(cycle_phase = factor(cycle_phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Anhedonia Score",
  predictor_var = "cycle_phase",
  outcome_var = "total",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("shaps_raincloud.png", plot = shaps_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

# PHQ plot
phq_plot <- raincloud_plot(
  dat = mhq_analysis %>%
    filter(measure == "phq") %>%
    mutate(cycle_phase = factor(cycle_phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Depression Score",
  predictor_var = "cycle_phase",
  outcome_var = "total",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("phq_raincloud.png", plot = phq_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

# STAI plot
stai_plot <- raincloud_plot(
  dat = mhq_analysis %>%
    filter(measure == "stai") %>%
    mutate(cycle_phase = factor(cycle_phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Anxiety Score",
  predictor_var = "cycle_phase",
  outcome_var = "total",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("stai_raincloud.png", plot = stai_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

pms_plot <- raincloud_plot(
  dat = mhq_analysis %>%
    filter(measure == "pms") %>%
    mutate(cycle_phase = factor(cycle_phase, levels = c("ovu", "men"))),
  title = "",
  xlab = "Cycle Phase",
  ylab = "Premenstrual Symptom Score",
  predictor_var = "cycle_phase",
  outcome_var = "total",
  predictor_tick_lab = c("Periovulatory", "Perimenstrual"),
  col = c("#00AB8E", "#CC0077"),
  direction = "vertical",
  include_grouping = FALSE,
  scale_seq = NULL,
  subj_var = "subj_id"
)
ggsave("pms_raincloud.png", plot = pms_plot, 
       width = 5, height = 4, dpi = 800, bg = "white")

# Combine plots
combined_plot <- (pms_plot | phq_plot) / (shaps_plot | stai_plot)

combined_styled <- combined_plot + 
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

ggsave("psychiatric_combined.png", 
       plot = combined_styled,
       width = 10, height = 8, dpi = 800, bg = "white")
