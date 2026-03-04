library(tidyverse)
library(ggdist)
library(scales)

## Survival vs sex & age -------------------------------------------------------

outF <- readRDS("results/modelF_tObs_aVR_itX_tR_noRecO_wYAFs.rds")
outM <- readRDS("results/modelM_tObs_aVR_itX_tR_noRecO_wYAFs_dexp1.rds")

# Posterior means vs year
years <- (1:n.occasions) + 2007
ageCs <- c("0", "1", "2", "3-6", "7-9", "10+")

# females
mcmc.musF <- outF %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

mcmc.musF <- mcmc.musF %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."),
               names_to = "param.full",
               values_to = "value")

mcmc.musF <- mcmc.musF %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"),
         a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs),
         sex = "Females")

summ.musF <- mcmc.musF %>%
  group_by(param, a) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            lcl = quantile(value, 0.025, na.rm = TRUE),
            ucl = quantile(value, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ageC = factor(ageCs[a], levels = ageCs))

# males
mcmc.musM <- outM %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

mcmc.musM <- mcmc.musM %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."),
               names_to = "param.full",
               values_to = "value")

mcmc.musM <- mcmc.musM %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"),
         a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs),
         sex = "Males")

summ.musM <- mcmc.musM %>%
  group_by(param, a) %>%
  summarise(mean = mean(value, na.rm = TRUE),
            lcl = quantile(value, 0.025, na.rm = TRUE),
            ucl = quantile(value, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(ageC = factor(ageCs[a], levels = ageCs))

# join
mcmc.mus <- rbind(mcmc.musF, mcmc.musM)

# survival vs age
blues   <- seq_gradient_pal("#08306B", "#9ECAE1")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(blues)   <- paste("Males", ageCs, sep = " - ")
names(oranges) <- paste("Females", ageCs, sep = " - ")

mcmc.mus %>% 
  filter(param == "mu.phi") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = sex), .width = c(0.5, 0.95), alpha = 0.6, colour = NA, position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = sex), point_interval = median_qi, size = 1.2, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("Males" = "#9ECAE1", "Females" = "#FDBE85"), name = NULL) +
  scale_colour_manual(values = c("Males" = "#08306B", "Females" = "#7F2704"), name = NULL) +
  labs(x = "Probability", y = "Age class") +
  theme_bw()

# ggsave("figures/modF_PHI&RK_lowRK_wVR&X.jpeg", width = 12.0, height = 24.0, units = c("cm"), dpi = 600)


## Road mortality vs sex & age -------------------------------------------------

mcmc.mus %>% 
  filter(param == "mu.R") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = sex), .width = c(0.5, 0.95), alpha = 0.6, colour = NA, position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = sex), point_interval = median_qi, size = 1.2, position = position_dodge(width = 0.6)) +
  scale_fill_manual(values = c("Males" = "#9ECAE1", "Females" = "#FDBE85"), name = NULL) +
  scale_colour_manual(values = c("Males" = "#08306B", "Females" = "#7F2704"), name = NULL) +
  labs(x = "Probability of death begin\nrelated to road mortality", y = "Age class") +
  theme_bw()

# ggsave("figures/modF_PHI&RK_lowRK_wVR&X.jpeg", width = 12.0, height = 24.0, units = c("cm"), dpi = 600)

## Survival & road mortality vs age --------------------------------------------

mcmc.mus <- out %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

mcmc.mus <- mcmc.mus %>% 
  select(starts_with("mu.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mu."),
               names_to = "param.full",
               values_to = "value")

mcmc.mus <- mcmc.mus %>% 
  mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+"),
         a = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  select(iter, param.full, param, a, value) %>% 
  mutate(ageC = factor(ageCs[a], levels = ageCs))

# summ.mus <- mcmc.mus %>% 
#   group_by(param, a) %>% 
#   summarise(mean = mean(value, na.rm = TRUE),
#             lcl = quantile(value, 0.025, na.rm = TRUE),
#             ucl = quantile(value, 0.975, na.rm = TRUE),
#             .groups = "drop") %>% 
#   mutate(ageC = factor(ageCs[a], levels = ageCs))

# survival & roadkill vs age
blues   <- seq_gradient_pal("#08306B", "#9ECAE1")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(blues)   <- paste("Survival", ageCs, sep = " - ")
names(oranges) <- paste("Road mortality", ageCs, sep = " - ")

mcmc.mus %>% 
  filter(param %in% c("mu.phi", "mu.R")) %>% 
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
  theme_classic()

# ggsave("figures/modF_PHI&RK_lowRK_wVR&X.jpeg", width = 12.0, height = 24.0, units = c("cm"), dpi = 600)

