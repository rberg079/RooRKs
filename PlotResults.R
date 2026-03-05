
# load libraries
library(bayesplot)
library(beepr)
library(coda)
library(cowplot)
library(forcats)
library(ggdist)
library(here)
library(lubridate)
library(MCMCvis)
library(nimble)
library(parallel)
library(patchwork)
library(postpack)
library(RColorBrewer)
library(readxl)
library(scales)
library(strex)
library(tidybayes)
library(tidyverse)


## FEMALES vs MALES by scenario ------------------------------------------------

# LOW RK
outF <- readRDS("results/modelF_tObs_aVR_itX_tR_noRecO_wYAFs.rds")
outM <- readRDS("results/modelM_tObs_aVR_itX_tR_noRecO_wYAFs_dexp1.rds")

# # HIGH RK
# outF <- readRDS("results/modelF_tObs_aVR_itX_tR_noRec_wYAFs.rds")
# outM <- readRDS("results/modelM_tObs_aVR_itX_tR_noRec_wYAFs_dexp1.rds")

n.occasions <- 18
n.ageC <- 6

# # CHECK raw data
# # for who was RKed when
# table(age[eh == 2])
# table(colnames(age)[col(age)[eh == 2]], age[eh == 2])

# Posterior means vs year
years <- (1:n.occasions) + 2007
ageCs <- c("0", "1", "2", "3-6", "7-9", "10+")


## Wrangle female data ---------------------------------------------------------

# convert model output into matrix
mcmc.musF <- outF %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

# calculate probability of being roadkilled
phi_cols <- grep("^mu.phi\\[", colnames(mcmc.musF))              # indices of survival columns
R_cols   <- sub("mu.phi", "mu.R", colnames(mcmc.musF)[phi_cols]) # corresponding mu.R columns
R_cols   <- match(R_cols, colnames(mcmc.musF))                   # grab indices

mu_die <- 1 - mcmc.musF[, phi_cols]                              # probability of dying
mu_RK <- mu_die * mcmc.musF[, R_cols]                            # probability of dying of roadkill

colnames(mu_die) <- sub("mu.phi", "mu.die", colnames(mu_die))
colnames(mu_RK)  <- sub("mu.phi", "mu.RK", colnames(mu_RK))

mcmc.musF <- cbind(mcmc.musF, mu_die, mu_RK)                     # add back to mcmc.mus
remove(mu_die, mu_RK)

# pivot long & sort by age
mcmc.musF <- mcmc.musF %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."),
               names_to = "param.full",
               values_to = "value") %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"),
         a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs),
         sex = "Females")

# summarise by age class
summ.musF <- mcmc.musF %>%
  group_by(param, a) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            lcl = quantile(value, 0.025, na.rm = TRUE),
            ucl = quantile(value, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ageC = factor(ageCs[a], levels = ageCs))


## Wrangle male data -----------------------------------------------------------

# convert model output into matrix
mcmc.musM <- outM %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

# calculate probability of being roadkilled
phi_cols <- grep("^mu.phi\\[", colnames(mcmc.musM))              # indices of survival columns
R_cols   <- sub("mu.phi", "mu.R", colnames(mcmc.musM)[phi_cols]) # corresponding mu.R columns
R_cols   <- match(R_cols, colnames(mcmc.musM))                   # grab indices

mu_die <- 1 - mcmc.musM[, phi_cols]                              # probability of dying
mu_RK <- mu_die * mcmc.musM[, R_cols]                            # probability of dying of roadkill

colnames(mu_die) <- sub("mu.phi", "mu.die", colnames(mu_die))
colnames(mu_RK)  <- sub("mu.phi", "mu.RK", colnames(mu_RK))

mcmc.musM <- cbind(mcmc.musM, mu_die, mu_RK)                     # add back to mcmc.mus
remove(mu_die, mu_RK)

# pivot long & sort by age
mcmc.musM <- mcmc.musM %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."),
               names_to = "param.full",
               values_to = "value") %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"),
         a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs),
         sex = "Males")

# summarise by age class
summ.musM <- mcmc.musM %>%
  group_by(param, a) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            lcl = quantile(value, 0.025, na.rm = TRUE),
            ucl = quantile(value, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ageC = factor(ageCs[a], levels = ageCs))

# join females & males
mcmc.mus <- rbind(mcmc.musF, mcmc.musM)


## PLOT ------------------------------------------------------------------------

# Sort colours
mcmc.mus <- mcmc.mus %>%
  mutate(sex_age = paste(sex, ageC, sep = " - "))

blues <- seq_gradient_pal("#08306B", "#9ECAE1")(seq(0, 1, length.out = n.ageC))
pinks <- seq_gradient_pal("#7A3B4B", "#E7B6C2")(seq(0, 1, length.out = n.ageC))
# oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(blues) <- paste("Males", ageCs, sep = " - ")
names(pinks) <- paste("Females", ageCs, sep = " - ")
# names(oranges) <- paste("Females", ageCs, sep = " - ")

colours_fill <- c(blues, pinks)

PHI <- mcmc.mus %>% 
  filter(param == "mu.phi") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = sex_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = sex), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Males" = "#072A5F", "Females" = "#6E3542"), name = NULL) +
  scale_x_continuous(limits = c(0.4, 1.0), breaks = c(0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Survival", y = "Age class") +
  theme_bw() +
  theme(axis.title.y = element_blank(),               # y axis title
        axis.title = element_text(size = 14),         # axis titles
        axis.text  = element_text(size = 12),         # axis tick labels
        legend.title = element_text(size = 14),       # legend title
        legend.text  = element_text(size = 12)); PHI  # legend labels

# ggsave("figures/F&M_PHI_lowRK.jpeg", width = 12.0, height = 20.0, units = c("cm"), dpi = 600)

RK <- mcmc.mus %>% 
  filter(param == "mu.RK") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = sex_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = sex), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6), show.legend = F) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Males" = "#072A5F", "Females" = "#6E3542"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.4), breaks = c(0.0, 0.2, 0.4)) +
  labs(x = "Road mortality", y = "Age class", title = "Low road mortality") +
  theme_bw() +
  theme(title = element_text(size = 16),             # title
        axis.title = element_text(size = 14),        # axis titles
        axis.text  = element_text(size = 12),        # axis tick labels
        legend.title = element_text(size = 14),      # legend title
        legend.text  = element_text(size = 12)); RK  # legend labels

# ggsave("figures/F&M_RK_lowRK.jpeg", width = 12.0, height = 20.0, units = c("cm"), dpi = 600)

# combine PHI & RK figures
RK + PHI + plot_layout(widths = c(0.4, 0.6))

# ggsave("figures/F&M_PHI&RK_lowRK.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)


## Survival & road mortality vs age --------------------------------------------

# survival & roadkill vs age
blues   <- seq_gradient_pal("#08306B", "#9ECAE1")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(blues)   <- paste("Survival", ageCs, sep = " - ")
names(oranges) <- paste("Road mortality", ageCs, sep = " - ")

mcmc.mus %>% 
  filter(sex == "Females",
         param %in% c("mu.phi", "mu.RK")) %>% 
  mutate(param_label = recode(param,
                              "mu.phi" = "Survival",
                              "mu.R"   = "Road mortality"),
         group = paste(param_label, ageC, sep = " - ")) %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = group), adjust = 1, .width = c(0.5, 0.95), alpha = 0.6, colour = NA) +
  stat_pointinterval(aes(colour = param_label), point_interval = median_qi, linewidth = 0.8, size = 1) +
  scale_fill_manual(values = c(blues, oranges), guide = "none") +
  scale_colour_manual(values = c("Survival" = "#08306B", "Road mortality" = "#7F2704"), name = NULL) +
  labs(x = "Probability", y = "Age class") +
  theme_bw()

# ggsave("figures/modF_PHI&RK_lowRK_wVR&X.jpeg", width = 12.0, height = 24.0, units = c("cm"), dpi = 600)


## LOW vs HIGH scenarios by sex ------------------------------------------------

# LOW RK
outF_L <- readRDS("results/modelF_tObs_aVR_itX_tR_noRecO_wYAFs.rds")
outM_L <- readRDS("results/modelM_tObs_aVR_itX_tR_noRecO_wYAFs_dexp1.rds")

# HIGH RK
outF_H <- readRDS("results/modelF_tObs_aVR_itX_tR_noRec_wYAFs.rds")
outM_H <- readRDS("results/modelM_tObs_aVR_itX_tR_noRec_wYAFs_dexp1.rds")

n.occasions <- 18
n.ageC <- 6

# Posterior means vs year
years <- (1:n.occasions) + 2007
ageCs <- c("0", "1", "2", "3-6", "7-9", "10+")

## Plot females ----------------------------------------------------------------

# convert model output into matrix
mcmc.musF_L <- outF_L %>% map(~as.data.frame(as.matrix(.x))) %>% bind_rows()
mcmc.musF_H <- outF_H %>% map(~as.data.frame(as.matrix(.x))) %>% bind_rows()

# calculate probability of being roadkilled
calculateRK <- function(mcmc_df) {
  
  # identify phi & R columns
  phi_cols <- grep("^mu.phi\\[", colnames(mcmc_df))
  R_cols <- sub("mu.phi", "mu.R", colnames(mcmc_df)[phi_cols])
  R_cols <- match(R_cols, colnames(mcmc_df))
  
  # compute probabilities
  mu_die <- 1 - mcmc_df[, phi_cols]
  mu_RK  <- mu_die * mcmc_df[, R_cols]
  
  # name columns & bind to original df
  colnames(mu_die) <- sub("mu.phi", "mu.die", colnames(mu_die))
  colnames(mu_RK)  <- sub("mu.phi", "mu.RK", colnames(mu_RK))
  mcmc_df <- cbind(mcmc_df, mu_die, mu_RK)
  
  return(mcmc_df)
}

mcmc.musF_L <- calculateRK(mcmc.musF_L)
mcmc.musF_H <- calculateRK(mcmc.musF_H)

# pivot long & sort by age
mcmc.musF_L <- mcmc.musF_L %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."), names_to = "param.full", values_to = "value") %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"), a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs), scenario = "Low")

mcmc.musF_H <- mcmc.musF_H %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."), names_to = "param.full", values_to = "value") %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"), a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs), scenario = "High")

# join low & high
mcmc.musF <- rbind(mcmc.musF_L, mcmc.musF_H) %>% 
  mutate(scenario_age = paste(scenario, ageC, sep = " - "))

# plot
greens   <- seq_gradient_pal("#555F45", "#A7B295")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(greens) <- paste("Low", ageCs, sep = " - ")
names(oranges) <- paste("High", ageCs, sep = " - ")

colours_fill <- c(greens, oranges)

PHI.F <- mcmc.musF %>% 
  filter(param == "mu.phi") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = scenario_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = scenario), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Low" = "#3F4733", "High" = "#631E03"), name = NULL) +
  scale_x_continuous(limits = c(0.4, 1.0), breaks = c(0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Survival", y = "Age class") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 14),           # axis titles
        axis.text  = element_text(size = 12),           # axis tick labels
        legend.title = element_text(size = 14),         # legend title
        legend.text  = element_text(size = 12)); PHI.F  # legend labels

# ggsave("figures/F&M_PHI_lowRK.jpeg", width = 12.0, height = 20.0, units = c("cm"), dpi = 600)

RK.F <- mcmc.musF %>% 
  filter(param == "mu.RK") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = scenario_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = scenario), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6), show.legend = F) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Low" = "#3F4733", "High" = "#631E03"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.2), breaks = c(0.0, 0.1, 0.2)) +
  labs(x = "Road mortality", y = "Age class", title = "Females") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),               # title
        axis.title = element_text(size = 14),          # axis titles
        axis.text  = element_text(size = 12),          # axis tick labels
        legend.title = element_text(size = 14),        # legend title
        legend.text  = element_text(size = 12)); RK.F  # legend labels

# ggsave("figures/F&M_RK_lowRK.jpeg", width = 12.0, height = 20.0, units = c("cm"), dpi = 600)

# combine PHI & RK figures
RK.F + PHI.F + plot_layout(widths = c(0.4, 0.6))

ggsave("figures/F_PHI&RK_LOWvsHIGH.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)


## Plot males ------------------------------------------------------------------

# convert model output into matrix
mcmc.musM_L <- outM_L %>% map(~as.data.frame(as.matrix(.x))) %>% bind_rows()
mcmc.musM_H <- outM_H %>% map(~as.data.frame(as.matrix(.x))) %>% bind_rows()

mcmc.musM_L <- calculateRK(mcmc.musM_L)
mcmc.musM_H <- calculateRK(mcmc.musM_H)

# pivot long & sort by age
mcmc.musM_L <- mcmc.musM_L %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."), names_to = "param.full", values_to = "value") %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"), a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs), scenario = "Low")

mcmc.musM_H <- mcmc.musM_H %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."), names_to = "param.full", values_to = "value") %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"), a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs), scenario = "High")

# join low & high
mcmc.musM <- rbind(mcmc.musM_L, mcmc.musM_H) %>% 
  mutate(scenario_age = paste(scenario, ageC, sep = " - "))

# plot
greens   <- seq_gradient_pal("#555F45", "#A7B295")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(greens) <- paste("Low", ageCs, sep = " - ")
names(oranges) <- paste("High", ageCs, sep = " - ")

colours_fill <- c(greens, oranges)

PHI.M <- mcmc.musM %>% 
  filter(param == "mu.phi") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = scenario_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = scenario), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Low" = "#3F4733", "High" = "#631E03"), name = NULL) +
  scale_x_continuous(limits = c(0.4, 1.0), breaks = c(0.4, 0.6, 0.8, 1.0)) +
  labs(x = "Survival", y = "Age class") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 14),           # axis titles
        axis.text  = element_text(size = 12),           # axis tick labels
        legend.title = element_text(size = 14),         # legend title
        legend.text  = element_text(size = 12)); PHI.M  # legend labels

# ggsave("figures/F&M_PHI_lowRK.jpeg", width = 12.0, height = 20.0, units = c("cm"), dpi = 600)

RK.M <- mcmc.musM %>% 
  filter(param == "mu.RK") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = scenario_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = scenario), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6), show.legend = F) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Low" = "#3F4733", "High" = "#631E03"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.4), breaks = c(0.0, 0.2, 0.4)) +
  labs(x = "Road mortality", y = "Age class", title = "Males") +
  theme_bw() +
  theme(title = element_text(size = 16),               # title
        axis.title = element_text(size = 14),          # axis titles
        axis.text  = element_text(size = 12),          # axis tick labels
        legend.title = element_text(size = 14),        # legend title
        legend.text  = element_text(size = 12)); RK.M  # legend labels

# ggsave("figures/F&M_RK_lowRK.jpeg", width = 12.0, height = 20.0, units = c("cm"), dpi = 600)

# combine PHI & RK figures
RK + PHI + plot_layout(widths = c(0.4, 0.6))

# ggsave("figures/M_PHI&RK_LOWvsHIGH.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)

# combine PHI & RK for females & males
(RK.F + PHI.F + plot_layout(widths = c(0.4, 0.6)))/(RK.M + PHI.M + plot_layout(widths = c(0.4, 0.6)))

# ggsave("figures/F&M_PHI&RK_LOWvsHIGH.jpeg", width = 20.0, height = 24.0, units = c("cm"), dpi = 600)

