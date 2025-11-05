# 19 February 2025
# Attempt to add age-specific survival rates,
# time-varying observation & recovery probabilities,
# fixed effect of environment (veg for now)
# & random year effect

library(readxl)
library(tidyverse)
library(lubridate)


## Set up ----------------------------------------------------------------------

# set toggles
females <- TRUE
testRun <- FALSE

# load libraries
library(forcats)
library(nimble)
library(here)
library(coda)
library(beepr)
library(cowplot)
library(patchwork)
library(RColorBrewer)
library(bayesplot)
library(tidybayes)
library(ggdist)
library(postpack)
library(strex)
library(MCMCvis)

# load data
source("PrepDataRK_M.R")
dataRK <- prepDataRK(females = females)
list2env(dataRK, envir = .GlobalEnv)

# # or...
# eh <- read_csv("eh.csv")
# age <- read_csv("age.csv")
# env <- read_csv("env.csv")


## Raw data checks -------------------------------------------------------------

# sex_labels <- c("0" = "Male", "1" = "Female")
# 
# for(s in 0:1){
#   cat("\n", sex_labels[as.character(s)], "\n")
#   
#   # Subset y for that sex
#   y.s <- y[sex == s, ]
#   
#   # 1. Observation frequencies
#   cat("\nObservation frequencies:\n")
#   print(table(as.vector(y.s), useNA = "ifany"))
#   
#   # 2. Roadkill vs. other mortality
#   rk.count <- sum(y.s == 3, na.rm = TRUE)
#   other.count <- sum(y.s == 4, na.rm = TRUE)
#   cat("\nroadkill:", rk.count, "other death:", other.count, "\n")
#   cat("roadkill / (roadkill + other death):",
#       round(rk.count / (rk.count + other.count), 3), "\n")
#   
#   # 3. On-site vs. off-site
#   on.site <- sum(y.s == 1, na.rm = TRUE)
#   off.site <- sum(y.s == 2, na.rm = TRUE)
#   cat("\non-site:", on.site, "off-site:", off.site, "\n")
#   cat("off-site / (off-site + on-site):",
#       round(off.site / (on.site + off.site), 3), "\n")
# }

table.obs <- table(factor(as.vector(y), levels = 1:5, exclude = NULL), useNA = "ifany")
print(table.obs)

rk.count <- sum(y == 3, na.rm = T)
other.count <- sum(y == 4, na.rm = T)
cat("roadkill:", rk.count, "other death:", other.count, "\n")
cat("roadkill / (roadkill + other death):", rk.count / (rk.count + other.count), "/n")

on.site <- sum(y == 1, na.rm = T)
off.site <- sum(y == 2, na.rm = T)
cat("on-site:", on.site, "off-site:", off.site, "/n")
cat("off-site / (off-site + on-site):", off.site / (on.site + off.site), "/n")


## Model -----------------------------------------------------------------------

code <- nimbleCode({
  
  # model missing environment
  for (m in 1:nNoVeg){
    veg[noVeg[m]] ~ dnorm(0, sd = 1)
  }
  
  # model age-class-specific survival probability
  # logit-linear function of time-dependent covariates
  # priors for mu.age defined further in the model 
  for (j in 1:(n.occasions-1)){
    eps.phi[j] ~ dnorm(0, tau.phi) # random year effect
    
    logit(phi.juv[j]) <- logit(mu.juv) + B.veg * veg[j] + eps.phi[j]
    logit(phi.sub[j]) <- logit(mu.sub) + B.veg * veg[j] + eps.phi[j]
    logit(phi.pri[j]) <- logit(mu.pri) + B.veg * veg[j] + eps.phi[j]
    logit(phi.pre[j]) <- logit(mu.pre) + B.veg * veg[j] + eps.phi[j]
    logit(phi.sen[j]) <- logit(mu.sen) + B.veg * veg[j] + eps.phi[j]
    
    mean.phi[1, j] <- phi.juv[j]
    mean.phi[2, j] <- phi.sub[j]
    mean.phi[3, j] <- phi.pri[j]
    mean.phi[4, j] <- phi.pre[j]
    mean.phi[5, j] <- phi.sen[j]
  }
  
  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      
      #### Likelihood ####
      z[i, t] ~ dcat(trans.mat[i, z[i, t-1], 1:n.true.states, t-1]) # z is latent
      y[i, t] ~ dcat(obs.mat[i, z[i, t], 1:n.obs.states, t]) # y is observed
      
    } # t
  } # i
      
  for (i in 1:n.inds){
    for (t in first[i]:(n.occasions-1)){
      
      #### Priors & constraints ####
      
      phi[i,t] <- mean.phi[ageC[age[i,t]], t] # survival
      R[i,t]   <- mean.R # mortality by roadkill
      M[i,t]   <- mean.M # migration (in or out)

      #### Transition matrix ####
      # 1 - alive on-site
      # 2 - alive off-site
      # 3 - dead by roadkill
      # 4 - dead by other
      # 5 - long dead
      
      # ALIVE ON-SITE
      trans.mat[i,1,1,t] <- phi[i,t]*(1-M[i,t])
      trans.mat[i,2,1,t] <- phi[i,t]*M[i,t]
      trans.mat[i,3,1,t] <- 0
      trans.mat[i,4,1,t] <- 0
      trans.mat[i,5,1,t] <- 0
      
      # ALIVE OFF-SITE
      trans.mat[i,1,2,t] <- phi[i,t]*M[i,t]
      trans.mat[i,2,2,t] <- phi[i,t]*(1-M[i,t])
      trans.mat[i,3,2,t] <- 0
      trans.mat[i,4,2,t] <- 0
      trans.mat[i,5,2,t] <- 0
      
      # DEAD BY ROADKILL
      trans.mat[i,1,3,t] <- (1-phi[i,t])*R[i,t]
      trans.mat[i,2,3,t] <- (1-phi[i,t])*R[i,t]
      trans.mat[i,3,3,t] <- 0
      trans.mat[i,4,3,t] <- 0
      trans.mat[i,5,3,t] <- 0
      
      # DEAD BY OTHER
      trans.mat[i,1,4,t] <- (1-phi[i,t])*(1-R[i,t])
      trans.mat[i,2,4,t] <- (1-phi[i,t])*(1-R[i,t])
      trans.mat[i,3,4,t] <- 0
      trans.mat[i,4,4,t] <- 0
      trans.mat[i,5,4,t] <- 0
      
      # LONG DEAD
      trans.mat[i,1,5,t] <- 0
      trans.mat[i,2,5,t] <- 0
      trans.mat[i,3,5,t] <- 1
      trans.mat[i,4,5,t] <- 1
      trans.mat[i,5,5,t] <- 1
      
    } # t
  } # i
      
  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      
      Pi[i,t]  <- mean.Pi[t]  # probability of observation, on-site
      Po[i,t]  <- mean.Po[t]  # probability of observation, off-site
      rR[i,t]  <- mean.rR[t]  # probability of recovery, roadkill
      rO[i,t]  <- mean.rO[t]  # probability of recovery, natural death
      
      #### Observation matrix ####
      # 1 - seen on-site
      # 2 - seen off-site
      # 3 - recovered roadkill
      # 4 - recovered other
      # 5 - undetected
      
      # SEEN ON-SITE
      obs.mat[i,1,1,t] <- Pi[i,t]
      obs.mat[i,2,1,t] <- 0
      obs.mat[i,3,1,t] <- 0
      obs.mat[i,4,1,t] <- 0
      obs.mat[i,5,1,t] <- 0
      
      # SEEN OFF-SITE
      obs.mat[i,1,2,t] <- 0
      obs.mat[i,2,2,t] <- Po[i,t]
      obs.mat[i,3,2,t] <- 0
      obs.mat[i,4,2,t] <- 0
      obs.mat[i,5,2,t] <- 0
      
      # RECOVERED ROADKILL
      obs.mat[i,1,3,t] <- 0
      obs.mat[i,2,3,t] <- 0
      obs.mat[i,3,3,t] <- rR[i,t]
      obs.mat[i,4,3,t] <- 0
      obs.mat[i,5,3,t] <- 0
      
      # RECOVERED OTHER
      obs.mat[i,1,4,t] <- 0
      obs.mat[i,2,4,t] <- 0
      obs.mat[i,3,4,t] <- 0
      obs.mat[i,4,4,t] <- rO[i,t]
      obs.mat[i,5,4,t] <- 0
      
      # UNDETECTED
      obs.mat[i,1,5,t] <- 1-Pi[i,t]
      obs.mat[i,2,5,t] <- 1-Po[i,t]
      obs.mat[i,3,5,t] <- 1-rR[i,t]
      obs.mat[i,4,5,t] <- 1-rO[i,t]
      obs.mat[i,5,5,t] <- 1
      
      # END priors & constraints
      
    } # t
  } # i
  
  # # quick simulations
  # hist(rbeta(1000, 8, 2))
  # hist(rbeta(1000, 4, 4))
  # hist(rbeta(1000, 2, 8))
  # hist(rbeta(1000, 1, 8))
  # hist(rbeta(1000, 1, 1))
  
  mu.juv   ~ dbeta(4, 4)
  mu.sub   ~ dbeta(8, 2)
  mu.pri   ~ dbeta(8, 2)
  mu.pre   ~ dbeta(8, 2)
  mu.sen   ~ dbeta(4, 4)
  
  mean.R   ~ dbeta(1, 8) 
  mean.M   ~ dbeta(1, 8) 
  
  for (t in 1:n.occasions){
    mean.Pi[t]  ~ dbeta(8, 2)
    mean.Po[t]  ~ dbeta(4, 4)
    mean.rR[t]  ~ dbeta(4, 4) # model is not sensitive to mean.rR & mean.rO priors
    mean.rO[t]  ~ dbeta(4, 4) # even fixing at 0.5 does not affect other params
  }
  
  B.veg ~ dnorm(0, 1)
  sigma.phi ~ dunif(0, 4)
  tau.phi <- 1 / (sigma.phi * sigma.phi)
  
}) # nimbleCode


## Initial values --------------------------------------------------------------

# data is coded as OBSERVATION states
# providing model with initial values for deterministic TRANSITIONS

# TRANSITION STATES
# 1 - alive on-site
# 2 - alive off-site
# 3 - dead by roadkill
# 4 - dead by other
# 5 - long dead

z_inits <- y
z_dat <- y
z_inits[y == 999] <- NA
z_dat[y == 999] <- NA

# seen on-site -> alive on-site (1 -> 1)
z_inits[y == 1] <- NA
z_dat[y == 1] <- 1

# seen off-site -> alive off-site (2 -> 2)
z_inits[y == 2] <- NA
z_dat[y == 2] <- 2

# recovered roadkill -> dead by roadkill (3 -> 3)
z_inits[y == 3] <- NA
z_dat[y == 3] <- 3

# recovered other -> dead by other (4 -> 4)
z_inits[y == 4] <- NA
z_dat[y == 4] <- 4

# undetected -> alive (5 -> 1)? TO DISCUSS
z_inits[y == 5] <- 1
z_dat[y == 5] <- NA

# undetected after mortality observed -> long dead (5 -> 5)
n.inds <- nrow(y)
for (i in 1:n.inds) {
  tmp <- z_dat[i, ]
  if (any(tmp == 3, na.rm = T) | any(tmp == 4, na.rm = T)) {
    indexD <- which(tmp == 3 | tmp == 4) # index of 3 or 4 (dead)
    indexLD <- (indexD + 1):length(tmp) # indices of 5 after 3|4 -> 5 (long dead)
    tmp[indexLD] <- 5
  }
  z_dat[i, ] <- tmp
}

z_inits[z_dat == 5] <- NA

# # SPECIAL CASES: fill in intermediate non-detections
# # IDs where applicable in females: 57, 88, 492, 832
# z_dat[37, 5:7] <- 2
# z_dat[55, 7:9] <- 2
# z_dat[177, 13:14] <- 2
# z_dat[242, 14] <- 2


## Assemble --------------------------------------------------------------------

# Inits
inits <- list(
  z         = z_inits,
  mu.juv    = rbeta(1, 4, 4),
  mu.sub    = rbeta(1, 8, 2),
  mu.pri    = rbeta(1, 8, 2),
  mu.pre    = rbeta(1, 8, 2),
  mu.sen    = rbeta(1, 4, 4),
  mean.R    = rbeta(1, 1, 8),
  mean.M    = rbeta(1, 1, 8),
  mean.Pi   = rbeta(n.occasions, 8, 2),
  mean.Po   = rbeta(n.occasions, 4, 4),
  mean.rR   = rbeta(n.occasions, 4, 4),
  mean.rO   = rbeta(n.occasions, 4, 4),
  B.veg     = 0,
  eps.phi   = rep(0, n.occasions-1),
  sigma.phi = runif(1, 0.5, 1.5)
  )

# Data
y[y == 999] <- NA
dat <- list(y = y, 
            z = z_dat, 
            age = age,
            ageC = ageC,
            veg = veg)

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("mu.juv", "mu.sub", "mu.pri", "mu.pre", "mu.sen",
            "mean.R", "mean.M", "mean.Pi", "mean.Po", "mean.rR", "mean.rO",
            "B.veg", "sigma.phi", "veg")

# Constants
const <- list(
  n.inds = n.inds,
  n.occasions = n.occasions,
  n.true.states = n.true.states,
  n.obs.states = n.obs.states,
  first = first,
  noVeg = noVeg,
  nNoVeg = nNoVeg)

# # Check that z[, first] is known for all inds...
# for (ii in 1:n.inds) {
#   print(z_dat[ii, first[ii]])
# }


## Run model -------------------------------------------------------------------

# MCMC settings
if(testRun){
  nb <- 0          # burn-in
  ni <- 10         # total iterations
  nt <- 1          # thin
  nc <- 3          # chains
}else{
  nb <- 5000       # burn-in
  ni <- 5000 + nb  # total iterations
  nt <- 1          # thin
  nc <- 3          # chains
}

# Compile configuration & build
Rmodel <- nimbleModel(code = code, constants = const, data = dat,
                      check = F, calculate = F, inits = inits)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt, useConjugacy = F)
Rmcmc <- buildMCMC(conf) # take a look at sampler

Cmodel <- compileNimble(Rmodel, showCompilerOutput = F)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
beep(sound = 2)

# Run MCMC
t.start <- Sys.time()
# sink("Results/errorreport.txt") # for debugging
# sink("errorreport.txt") # for debugging
out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits = inits,
               setSeed = F, progressBar = T, samplesAsCodaMCMC = T)
t.end <- Sys.time()
# sink() # closing txt file
(runTime <- t.end - t.start)
beep(sound = 2)




## Plots -----------------------------------------------------------------------

# Save model, outputs, & diagnostics
MCMCdiag(out, dir = "./Results", save_object = TRUE,
         obj_name = "modelM_varObs_veg.rds", file_name = "modelM_varObs_veg_summary.txt")
# mod <- readRDS(here("Results", "model1.rds")
# mod_summary <- read.lines(here("Results", "model_summary.txt"))

# Calculate numerical summaries
model.summary <- MCMCsummary(object = out, round = 3)

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

mcmc.df <- as.tibble(as.matrix(out)) %>% 
  select(starts_with("mean.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mean."),
               names_to = "param.full",
               values_to = "value")

mcmc.df <- mcmc.df %>% 
  mutate(param = str_extract(param.full, "mean\\.[A-Za-z]+"),
         t = as.numeric(str_extract(param.full, "\\d+"))) %>% 
  filter(!is.na(t))

summaries <- mcmc.df %>% 
  group_by(param, t) %>% 
  summarise(mean = mean(value, na.rm = T),
            lcl = quantile(value, 0.025, na.rm = T),
            ucl = quantile(value, 0.975, na.rm = T),
            .groups = "drop") %>% 
  mutate(year = years[t])

ggplot(summaries, aes(x = year, y = mean, fill = param, colour = param)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
  facet_wrap(~param, scales = "free_y") +
  labs(x = "Year", y = "Posterior mean (±95% CrI)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90", colour = NA))

library(coda)
c1 <- as.mcmc(out$chain1)

# Caterpillar plot posterior distribution of phi
# Posterior median, 50% & 95% credible interval
MCMCplot(object = out,
         horiz = FALSE,
         rank = TRUE,
         ref_ovl = FALSE)

# Check convergence of chains
MCMCtrace(object = out,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE, # add Rhat diagnostics
          n.eff = TRUE, # add eff sample size
          params = c("mu.juv", "mu.sub", "mu.pri", "mu.pre", "mu.sen",
                     "mean.R", "mean.M", "mean.Pi", "mean.Po", "mean.rR", "mean.rO"))

# Correlation plots
autocorr.diag(out)
# autocorr.plot(out)
coda::crosscorr.plot(out)

# Scatterplot of variables
color_scheme_set("purple")
mcmc_scatter(out, pars = c("mu.pri", "mean.R"))

