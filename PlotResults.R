# 4 March 2026
# Round up & plot road mortality results


## Set up ----------------------------------------------------------------------

# load libraries
library(ggdist)
library(patchwork)
library(readxl)
library(scales)
library(tidyverse)

# LOW RK
outFL <- readRDS("results/modelF_tObs_aVR_itX_tR_noRecO_wYAFs.rds")
outML <- readRDS("results/modelM_tObs_aVR_itX_tR_noRecO_wYAFs_dexp1.rds")

# HIGH RK
outFH <- readRDS("results/modelF_tObs_aVR_itX_tR_noRec_wYAFs.rds")
outMH <- readRDS("results/modelM_tObs_aVR_itX_tR_noRec_wYAFs_dexp1.rds")

n.occasions <- 18
n.ageC <- 6

years <- (1:n.occasions) + 2007
ageCs <- c("0", "1", "2", "3-6", "7-9", "10+")

# # CHECK raw data
# # for who was RKed when
# table(age[eh == 2])
# table(colnames(age)[col(age)[eh == 2]], age[eh == 2])


## Wrangle mus -----------------------------------------------------------------

# randomly sample iterations & bind chains
sampleChains <- function(mcmc_list, n = 1000){
  
  map(mcmc_list, function(chain){
    mat <- as.matrix(chain)
    mat[sample(seq_len(nrow(mat)), n), , drop = F]
  }) %>%
    map(as.data.frame) %>% 
    bind_rows()
}

meansFL <- sampleChains(outFL, n = 1000)
meansML <- sampleChains(outML, n = 1000)

meansFH <- sampleChains(outFH, n = 1000)
meansMH <- sampleChains(outMH, n = 1000)

# calculate probability of being roadkilled
calculateRK <- function(mcmc_df){
  
  # identify phi & R columns
  phi_cols <- grep("^mu.phi\\[", colnames(mcmc_df))
  R_cols   <- sub("mu.phi", "mu.R", colnames(mcmc_df)[phi_cols])
  R_cols   <- match(R_cols, colnames(mcmc_df))
  
  # compute probabilities
  mu_die <- 1 - mcmc_df[, phi_cols]
  mu_RK <- mu_die * mcmc_df[, R_cols]
  
  # name columns & bind to original df
  colnames(mu_die) <- sub("mu.phi", "mu.die", colnames(mu_die))
  colnames(mu_RK)  <- sub("mu.phi", "mu.RK", colnames(mu_RK))
  mcmc_df <- cbind(mcmc_df, mu_die, mu_RK)
  
  return(mcmc_df)
}

meansFL <- calculateRK(meansFL)
meansML <- calculateRK(meansML)

meansFH <- calculateRK(meansFH)
meansMH <- calculateRK(meansMH)

# pivot & parse parameters
parseMCMC <- function(mcmc_df, sex_label, scenario_label){
  
  mcmc_df %>% 
    select(starts_with("mu.")) %>% 
    mutate(iter = row_number()) %>% 
    pivot_longer(cols = starts_with("mu."),
                 names_to = "param.full",
                 values_to = "value") %>% 
    mutate(param = str_extract(param.full, "mu\\.[A-Za-z]+")) %>% 
    filter(param != "mu.p") %>% 
    mutate(a = as.numeric(str_extract(param.full, "\\d+")),
           # a = str_extract_all(param.full, "\\d+"),
           sex = sex_label,
           scenario = scenario_label,
           ageC = factor(ageCs[a], levels = ageCs)) %>% 
    select(iter, scenario, sex, param.full, param, ageC, a, value)
}

meansFL <- parseMCMC(meansFL, "Females", "Low")
meansML <- parseMCMC(meansML, "Males", "Low")

meansFH <- parseMCMC(meansFH, "Females", "High")
meansMH <- parseMCMC(meansMH, "Males", "High")

# calculate means & 95% CrIs
results <- function(df){
  
  summary <- df %>%
    group_by(param, a) %>%
    summarise(mean = mean(value, na.rm = T),
              lcl = quantile(value, 0.025, na.rm = T),
              ucl = quantile(value, 0.975, na.rm = T),
              .groups = "drop") %>%
    mutate(ageC = factor(ageCs[a], levels = ageCs))
  
  return(summary)
}

summFL <- results(meansFL)
summML <- results(meansML)

summFH <- results(meansFH)
summMH <- results(meansMH)


## Wrangle means ----------------------------------------------------------------

# randomly sample iterations & bind chains
sampleChains <- function(mcmc_list, n = 1000){
  
  map(mcmc_list, function(chain){
    mat <- as.matrix(chain)
    mat[sample(seq_len(nrow(mat)), n), , drop = F]
  }) %>%
    map(as.data.frame) %>% 
    bind_rows()
}

meansFL <- sampleChains(outFL, n = 1000)
meansML <- sampleChains(outML, n = 1000)

meansFH <- sampleChains(outFH, n = 1000)
meansMH <- sampleChains(outMH, n = 1000)

# calculate probability of being roadkilled
calculateRK <- function(mcmc_df){
  
  # identify phi & R columns
  phi_cols <- grep("^mean.phi\\[", colnames(mcmc_df))
  R_cols   <- sub("mean.phi", "mean.R", colnames(mcmc_df)[phi_cols])
  R_cols   <- match(R_cols, colnames(mcmc_df))
  
  # compute probabilities
  mean_die <- 1 - mcmc_df[, phi_cols]
  mean_RK <- mean_die * mcmc_df[, R_cols]
  
  # name columns & bind to original df
  colnames(mean_die) <- sub("mean.phi", "mean.die", colnames(mean_die))
  colnames(mean_RK)  <- sub("mean.phi", "mean.RK", colnames(mean_RK))
  mcmc_df <- cbind(mcmc_df, mean_die, mean_RK)
  
  return(mcmc_df)
}

meansFL <- calculateRK(meansFL)
meansML <- calculateRK(meansML)

meansFH <- calculateRK(meansFH)
meansMH <- calculateRK(meansMH)

# pivot & parse parameters
parseMCMC <- function(mcmc_df, sex_label, scenario_label){
  
  mcmc_df %>% 
    select(starts_with("mean.")) %>% 
    mutate(iter = row_number()) %>% 
    pivot_longer(cols = starts_with("mean."),
                 names_to = "param.full",
                 values_to = "value") %>% 
    mutate(param = str_extract(param.full, "mean\\.[A-Za-z]+")) %>% 
    filter(param != "mean.p") %>% 
    mutate(index = str_extract_all(param.full, "\\d+"),
           a = map_dbl(index, ~as.numeric(.x[1])),
           t = map_dbl(index, ~as.numeric(.x[2])),
           sex = sex_label,
           scenario = scenario_label,
           ageC = factor(ageCs[a], levels = ageCs)) %>% 
    select(iter, scenario, sex, param.full, param, ageC, a, t, value)
}

meansFL <- parseMCMC(meansFL, "Females", "Low")
meansML <- parseMCMC(meansML, "Males", "Low")

meansFH <- parseMCMC(meansFH, "Females", "High")
meansMH <- parseMCMC(meansMH, "Males", "High")

# calculate means & 95% CrIs
# # FOR AGE EFFECTS ONLY
# results <- function(df){
# 
#   summary <- df %>%
#     group_by(param, a) %>%
#     summarise(mean = mean(value, na.rm = T),
#               lcl = quantile(value, 0.025, na.rm = T),
#               ucl = quantile(value, 0.975, na.rm = T),
#               .groups = "drop") %>%
#     mutate(ageC = factor(ageCs[a], levels = ageCs))
# 
#   return(summary)
# }

# FOR TIME SERIES
results <- function(df){
  
  summary <- df %>%
    group_by(param, a, t) %>%
    summarise(mean = mean(value, na.rm = T),
              lcl = quantile(value, 0.025, na.rm = T),
              ucl = quantile(value, 0.975, na.rm = T),
              .groups = "drop") %>%
    mutate(ageC = factor(ageCs[a], levels = ageCs),
           year = years[t])
  
  return(summary)
}

summFL <- results(meansFL)
summML <- results(meansML)

summFH <- results(meansFH)
summMH <- results(meansMH)


## xMed effect -----------------------------------------------------------------

# randomly sample iterations & bind chains
sampleChains <- function(mcmc_list, n = 1000){
  
  map(mcmc_list, function(chain){
    mat <- as.matrix(chain)
    mat[sample(seq_len(nrow(mat)), n), , drop = F]
  }) %>%
    map(as.data.frame) %>% 
    bind_rows()
}

meansFL <- sampleChains(outFL, n = 1000)
meansML <- sampleChains(outML, n = 1000)

meansFH <- sampleChains(outFH, n = 1000)
meansMH <- sampleChains(outMH, n = 1000)

# calculate probability of being roadkilled
calculateRK <- function(mcmc_df){
  
  # identify phi & R columns
  phi_cols <- grep("^mean.phi\\[", colnames(mcmc_df))
  R_cols   <- sub("mean.phi", "mean.R", colnames(mcmc_df)[phi_cols])
  R_cols   <- match(R_cols, colnames(mcmc_df))
  
  # compute probabilities
  mean_die <- 1 - mcmc_df[, phi_cols]
  mean_RK <- mean_die * mcmc_df[, R_cols]
  
  # name columns & bind to original df
  colnames(mean_die) <- sub("mean.phi", "mean.die", colnames(mean_die))
  colnames(mean_RK)  <- sub("mean.phi", "mean.RK", colnames(mean_RK))
  mcmc_df <- cbind(mcmc_df, mean_die, mean_RK)
  
  return(mcmc_df)
}

meansFL <- calculateRK(meansFL)
meansML <- calculateRK(meansML)

meansFH <- calculateRK(meansFH)
meansMH <- calculateRK(meansMH)

xGrid <- seq(-2.343521, 2.169987, length.out = 100)

predictRK <- function(mcmc_df, xGrid, sex_label, scenario_label) {
  
  betaX_phi <- mcmc_df$betaX.phi
  betaX_R   <- mcmc_df$betaX.R
  preds <- list()
  
  for(a in 1:n.ageC){
    phi_cols <- grep(paste0("^mean\\.phi\\[", a, ", "), colnames(mcmc_df), value = TRUE)
    R_cols   <- grep(paste0("^mean\\.R\\[",   a, ", "), colnames(mcmc_df), value = TRUE)
    
    base_phi <- as.matrix(mcmc_df[, phi_cols])
    base_R   <- as.matrix(mcmc_df[, R_cols])
    
    for(x in xGrid) {
      phi_mat <- plogis(qlogis(base_phi) + betaX_phi * x)
      R_mat   <- plogis(qlogis(base_phi) + betaX_phi * x)
      RK_iter <- rowMeans((1 - phi_mat) * R_mat)
      
      preds[[length(preds) + 1]] <- data.frame(
        ageC = factor(ageCs[a], levels = ageCs),
        xmed = x,
        mean = mean(RK_iter),
        lcl  = quantile(RK_iter, 0.025),
        ucl  = quantile(RK_iter, 0.975),
        scenario = scenario_label,
        sex = sex_label
      )
    }
  }
  return(do.call(rbind, preds))
}

# generate predictions
predFL <- predictRK(meansFL, xGrid, "Females", "Low")
predML <- predictRK(meansML, xGrid, "Males",   "Low")

predFH <- predictRK(meansFH, xGrid, "Females", "High")
predMH <- predictRK(meansMH, xGrid, "Males",   "High")

preds <- bind_rows(predFL, predML, predFH, predMH)
preds$scenario <- factor(preds$scenario, levels = c("Low", "High"))

# plot
xEffect <- preds %>% 
  ggplot(., aes(x = xmed, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  facet_grid(sex ~ scenario) +
  labs(x = "Standardized median easting (X)",
       y = "Road mortality",
       color = "Age class",
       fill = "Age class") +
  scale_y_continuous(limits = c(0, 0.25),
                     breaks = c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); xEffect

ggsave("figures/combXEffect.jpeg", width = 24.0, height = 18.0, units = c("cm"), dpi = 600)

# # separately
# X_FL <- predFL %>% 
#   ggplot(., aes(x = xmed, y = mean)) +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
#   geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
#   labs(title = "Females",
#        x = "Standardized median easting (X)",
#        y = "Road mortality (lowest)",
#        color = "Age class",
#        fill = "Age class") +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = "grey90"),
#         title = element_text(size = 16),
#         axis.title.x = element_blank(),
#         axis.title = element_text(size = 14),
#         axis.text  = element_text(size = 12)); X_FL
# 
# X_ML <- predML %>% 
#   ggplot(., aes(x = xmed, y = mean)) +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
#   geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
#   labs(title = "Males",
#        x = "Standardized median easting (X)",
#        y = "Road mortality (lowest)",
#        color = "Age class",
#        fill = "Age class") +
#   theme_bw() +
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         title = element_text(size = 16),
#         axis.title = element_text(size = 14),
#         axis.text  = element_text(size = 12),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = "grey90")); X_ML
# 
# X_FH <- predFH %>% 
#   ggplot(., aes(x = xmed, y = mean)) +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
#   geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
#   labs(x = "Standardized median easting (X)",
#        y = "Road mortality (highest)",
#        color = "Age class",
#        fill = "Age class") +
#   theme_bw() +
#   theme(axis.title = element_text(size = 14),
#         axis.text  = element_text(size = 12),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = "grey90")); X_FH
# 
# X_MH <- predMH %>% 
#   ggplot(., aes(x = xmed, y = mean)) +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
#   geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
#   labs(x = "Standardized median easting (X)",
#        y = "Road mortality (highest)",
#        color = "Age class",
#        fill = "Age class") +
#   theme_bw() +
#   theme(axis.title.y = element_blank(),
#         axis.title = element_text(size = 14),
#         axis.text  = element_text(size = 12),,
#         legend.title = element_text(size = 14),
#         legend.text  = element_text(size = 12),
#         panel.grid.minor = element_blank(),
#         strip.background = element_rect(fill = "grey90")); X_MH
# 
# library(patchwork)
# (X_FL + X_ML) / (X_FH + X_MH)


## FEMALES vs MALES ------------------------------------------------------------

# join females & males
meansL <- rbind(meansFL, meansML) %>% 
  mutate(sex_age = paste(sex, ageC, sep = " - "))

meansH <- rbind(meansFH, meansMH) %>% 
  mutate(sex_age = paste(sex, ageC, sep = " - "))

# sort colours
blues <- seq_gradient_pal("#08306B", "#9ECAE1")(seq(0, 1, length.out = n.ageC))
pinks <- seq_gradient_pal("#7A3B4B", "#E7B6C2")(seq(0, 1, length.out = n.ageC))

names(blues) <- paste("Males", ageCs, sep = " - ")
names(pinks) <- paste("Females", ageCs, sep = " - ")

colours_fill <- c(blues, pinks)

# plot
PHI_L <- meansL %>% 
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
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); PHI_L

PHI_H <- meansH %>% 
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
  theme(axis.title.y = element_blank(),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); PHI_H

RK_L <- meansL %>% 
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
  theme(title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); RK_L

RK_H <- meansH %>% 
  filter(param == "mu.RK") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = sex_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = sex), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6), show.legend = F) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Males" = "#072A5F", "Females" = "#6E3542"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.4), breaks = c(0.0, 0.2, 0.4)) +
  labs(x = "Road mortality", y = "Age class", title = "High road mortality") +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); RK_H

# combine PHI & RK plots
RK_L + PHI_L + plot_layout(widths = c(0.4, 0.6))
# ggsave("figures/F&M_PHI&RK_lowRK.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)

RK_H + PHI_H + plot_layout(widths = c(0.4, 0.6))
# ggsave("figures/F&M_PHI&RK_lowRK.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)

(RK_L + PHI_L + plot_layout(widths = c(0.4, 0.6)))/(RK_H + PHI_H + plot_layout(widths = c(0.4, 0.6)))
# ggsave("figures/F&M_PHI&RK_LOWvsHIGH.jpeg", width = 20.0, height = 28.0, units = c("cm"), dpi = 600)


## LOW vs HIGH SCENARIOS -------------------------------------------------------

# join low & high scenarios
meansF <- rbind(meansFL, meansFH) %>% 
  mutate(scenario_age = paste(scenario, ageC, sep = " - "))

meansM <- rbind(meansML, meansMH) %>% 
  mutate(scenario_age = paste(scenario, ageC, sep = " - "))

# sort colours
greens  <- seq_gradient_pal("#555F45", "#A7B295")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(greens) <- paste("Low", ageCs, sep = " - ")
names(oranges) <- paste("High", ageCs, sep = " - ")

colours_fill <- c(greens, oranges)

# plot
PHI_F <- meansF %>% 
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
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); PHI_F

PHI_M <- meansM %>% 
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
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); PHI_M

RK_F <- meansF %>% 
  filter(param == "mu.RK") %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = scenario_age), .width = c(0.5, 0.95), alpha = 0.6, colour = NA,
               position = position_dodge(width = 0.6)) +
  stat_pointinterval(aes(colour = scenario), point_interval = median_qi, size = 1.2,
                     position = position_dodge(width = 0.6), show.legend = F) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Low" = "#3F4733", "High" = "#631E03"), name = NULL) +
  scale_x_continuous(limits = c(0, 0.4), breaks = c(0.0, 0.2, 0.4)) +
  labs(x = "Road mortality", y = "Age class", title = "Females") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); RK_F

RK_M <- meansM %>% 
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
  theme(title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); RK_M

# combine
RK_F + PHI_F + plot_layout(widths = c(0.4, 0.6))
# ggsave("figures/F_PHI&RK_LOWvsHIGH.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)

RK_M + PHI_M + plot_layout(widths = c(0.4, 0.6))
# ggsave("figures/F_PHI&RK_LOWvsHIGH.jpeg", width = 24.0, height = 20.0, units = c("cm"), dpi = 600)

(RK_F + PHI_F + plot_layout(widths = c(0.4, 0.6)))/(RK_M + PHI_M + plot_layout(widths = c(0.4, 0.6)))
ggsave("figures/F&M_PHI&RK_LOWvsHIGH.jpeg", width = 20.0, height = 28.0, units = c("cm"), dpi = 600)


## Combined PHI & RK -----------------------------------------------------------

# sort colours
blues   <- seq_gradient_pal("#08306B", "#9ECAE1")(seq(0, 1, length.out = n.ageC))
oranges <- seq_gradient_pal("#7F2704", "#FDBE85")(seq(0, 1, length.out = n.ageC))

names(blues)   <- paste("Survival", ageCs, sep = " - ")
names(oranges) <- paste("Road mortality", ageCs, sep = " - ")

colours_fill <- c(blues, oranges)

# plot
comb_F <- meansF %>% 
  filter(param %in% c("mu.phi", "mu.RK")) %>% 
  mutate(param_label = recode(param,
                              "mu.phi" = "Survival",
                              "mu.RK"   = "Road mortality"),
         group = paste(param_label, ageC, sep = " - ")) %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = group), adjust = 1, .width = c(0.5, 0.95), alpha = 0.6, colour = NA) +
  stat_pointinterval(aes(colour = param_label), point_interval = median_qi, linewidth = 0.8, size = 1) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Survival" = "#06224C", "Road mortality" = "#631E03"), name = NULL) +
  labs(x = "Probability", y = "Age class", title = "Females") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); comb_F

comb_M <- meansM %>% 
  filter(param %in% c("mu.phi", "mu.RK")) %>% 
  mutate(param_label = recode(param,
                              "mu.phi" = "Survival",
                              "mu.RK"   = "Road mortality"),
         group = paste(param_label, ageC, sep = " - ")) %>% 
  ggplot(aes(x = value, y = ageC)) +
  stat_halfeye(aes(fill = group), adjust = 1, .width = c(0.5, 0.95), alpha = 0.6, colour = NA) +
  stat_pointinterval(aes(colour = param_label), point_interval = median_qi, linewidth = 0.8, size = 1) +
  scale_fill_manual(values = colours_fill, guide = "none") +
  scale_colour_manual(values = c("Survival" = "#06224C", "Road mortality" = "#631E03"), name = NULL) +
  labs(x = "Probability", y = "Age class", title = "Males") +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)); comb_M

# combine
comb_F / comb_M
# ggsave("figures/modF_PHI&RK_lowRK_wVR&X.jpeg", width = 12.0, height = 24.0, units = c("cm"), dpi = 600)


## Time series -----------------------------------------------------------------

# roadkill vs year
RKt_FL <- summFL %>% 
  filter(param %in% c("mean.RK")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (lowest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = NA)); RKt_FL

RKt_FH <- summFH %>% 
  filter(param %in% c("mean.RK")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (highest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = NA)); RKt_FH

RKt_ML <- summML %>% 
  filter(param %in% c("mean.RK")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (lowest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = NA)); RKt_ML

RKt_MH <- summMH %>% 
  filter(param %in% c("mean.RK")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (highest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = NA)); RKt_MH

# survival vs year
PHIt_FL <- summFL %>% 
  filter(param %in% c("mean.phi")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Survival", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = NA)); PHIt_FL

PHIt_ML <- summML %>% 
  filter(param %in% c("mean.phi")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Survival", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(title = element_text(size = 16),
        axis.title = element_text(size = 14),
        axis.text  = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12),
        strip.background = element_rect(fill = "grey90", colour = NA)); PHIt_ML

library(patchwork)
RKt_FL / RKt_FH / PHIt_FL
# ggsave("figures/modF_combTimeSeries.jpeg", width = 20.0, height = 28.0, units = c("cm"), dpi = 600)

RKt_ML / RKt_MH / PHIt_ML
# ggsave("figures/modM_combTimeSeries.jpeg", width = 20.0, height = 28.0, units = c("cm"), dpi = 600)


## Original plots --------------------------------------------------------------

# LOW RK
out <- readRDS("results/modelF_tObs_aVR_itX_tR_noRecO_wYAFs.rds")
out <- readRDS("results/modelM_tObs_aVR_itX_tR_noRecO_wYAFs_dexp1.rds")

# HIGH RK
out <- readRDS("results/modelF_tObs_aVR_itX_tR_noRec_wYAFs.rds")
out <- readRDS("results/modelM_tObs_aVR_itX_tR_noRec_wYAFs_dexp1.rds")

# # Posterior predictive checks
# set.seed(1)
# samples.mat <- do.call(rbind, lapply(out, as.matrix))
# s <- samples.mat[sample(nrow(samples.mat), 500),]
# 
# simulate.y <- function(par.row){
#   mean.R   <- par.row["mean.R"]
#   mean.rR  <- par.row["mean.rR"]
#   mean.rO  <- par.row["mean.rO"]
# 
#   tot.dead <- sum(y == 3 | y == 4, na.rm = T)
#   sim.rR <- rbinom(1, tot.dead, mean.R * mean.rR / (mean.R * mean.rR + (1-mean.R) * mean.rO))
#   return(c(sim.rR = sim.rR))
# }
# 
# sim.vals <- t(apply(s, 1, simulate.y))
# hist(sim.vals, main = "Simulated recovered roadkills")
# abline(v = sum(y == 3, na.rm = T), col = "red", lwd = 2)
# # appears to be simulating the right number of roadkills...

# Posterior means vs year
years <- (1:n.occasions) + 2007
ageCs <- c("0", "1", "2", "3-6", "7-9", "10+")

mcmc.means <- out %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

mcmc.means <- mcmc.means %>% 
  select(starts_with("mean.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mean."),
               names_to = "param.full",
               values_to = "value")

mcmc.means <- mcmc.means %>% 
  mutate(param = str_extract(param.full, "mean\\.[A-Za-z]+"),
         # extract all numbers inside brackets
         index = str_extract_all(param.full, "\\d+"),
         index1 = map_dbl(index, ~as.numeric(.x[1])),
         index2 = map_dbl(index, ~ifelse(length(.x) > 1, as.numeric(.x[2]), NA_real_)),
         
         # identify parameter dimensions
         is_time = param %in% c("mean.p"),
         is_both = param %in% c("mean.phi", "mean.R"),
         
         # assign t & a depending on parameter
         t = case_when(is_time ~ index1, is_both ~ index2),
         a = case_when(is_both ~ index1, TRUE ~ NA_real_)) %>% 
  select(iter, param.full, param, a, t, value)

summ.means <- mcmc.means %>% 
  group_by(param, a, t) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            lcl = quantile(value, 0.025, na.rm = TRUE),
            ucl = quantile(value, 0.975, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(year = years[t],
         ageC = factor(ageCs[a], levels = ageCs))

# # observation vs year
# summ.means %>% 
#   filter(param %in% c("mean.p")) %>% 
#   ggplot(., aes(x = year, y = mean, fill = param, colour = param)) +
#   geom_line(linewidth = 1) +
#   geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
#   facet_wrap(~param, scales = "free_y") +
#   labs(x = "Year", y = "Posterior mean (±95% CrI)") +
#   ylim(0, 1) +
#   theme_bw() +
#   theme(legend.position = "none",
#         strip.background = element_rect(fill = "grey90", colour = NA))

# roadkill vs year
lowF <- summ.means %>% 
  filter(param %in% c("mean.R")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (lowest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA)); lowF

# ggsave("figures/modM_RlowRK_wVR&X.jpeg", width = 24.0, height = 12.0, units = c("cm"), dpi = 600)

# roadkill vs year
highF <- summ.means %>% 
  filter(param %in% c("mean.R")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (highest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA)); highF

# ggsave("figures/modM_RlowRK_wVR&X.jpeg", width = 24.0, height = 12.0, units = c("cm"), dpi = 600)

# survival vs year
surv <- summ.means %>% 
  filter(param %in% c("mean.phi")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Survival", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", colour = NA)); surv

# ggsave("figures/PHIlowRK_wVR&X.jpeg", width = 24.0, height = 12.0, units = c("cm"), dpi = 600)

library(patchwork)
lowF / highF / surv
# ggsave("figures/modM_combTimeSeries.jpeg", width = 20.0, height = 28.0, units = c("cm"), dpi = 600)

# # Combined figure
# library(patchwork)
# RK.V + RK.VD + RK.VR # + plot_layout(nrow = 1)
# ggsave("figures/RKvsENV.jpeg", width = 36.0, height = 12.0, units = c("cm"), dpi = 600)
# 
# S.V + S.VD + S.VR # + plot_layout(nrow = 1)
# ggsave("figures/SvsENV.jpeg", width = 36.0, height = 12.0, units = c("cm"), dpi = 600)

# Checks
library(coda)
c1 <- as.mcmc(out$chain1)

MCMCplot(object = out,
         horiz = FALSE,
         rank = TRUE,
         ref_ovl = FALSE)

MCMCtrace(object = out,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE, # add Rhat diagnostics
          n.eff = TRUE, # add eff sample size
          params = c("betaD.phi", "betaV.phi", "betaD.R", "betaV.R",
                     "mu.phi", "mu.R", "mu.p",
                     "sigma.phi", "sigma.R", "sigma.p"))

# Correlation
autocorr.diag(out)
# autocorr.plot(out)
coda::crosscorr.plot(out)

