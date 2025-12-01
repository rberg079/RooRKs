# 28 October 2025
# Attempt to add age-specific environmental effects
# & perhaps age-specific random year effects


## Set up ----------------------------------------------------------------------

# set toggles
females <- TRUE
testRun <- FALSE
parallelRun <- TRUE

# name outputs
out.model <- "modelF_tObs_aVeg_atMR_obs.rds"
out.sum <- "modelF_tObs_aVeg_atMR_obs_sum.txt"

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
library(strex)
library(tidybayes)
library(tidyverse)

# load data
source("PrepDataRK_M.R")
dataRK <- prepDataRK(females = females)
list2env(dataRK, envir = .GlobalEnv)

# # or...
# eh <- read_csv("eh.csv")
# age <- read_csv("age.csv")
# env <- read_csv("env.csv")

# # or... fetch simulated data
# source("SimDataRK.R")
# dataRK <- simulateDataRK(mu.age = c(0.4, 0.6, 0.6, 0.6, 0.4),
#                          B.veg = c(1, 0.8, 0.6, 0.6, 1),
#                          sigma.phi = 1,
#                          mean.R = 0.1,
#                          mean.M = 0.1,
#                          mean.Pi = NULL,
#                          mean.Po = NULL,
#                          mean.rR = NULL,
#                          mean.rO = NULL,
#                          seed = 123)
# list2env(dataRK, envir = .GlobalEnv)


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

# table.obs <- table(factor(as.vector(y), levels = 1:5, exclude = NULL), useNA = "ifany")
# print(table.obs)
# 
# rk.count <- sum(y == 3, na.rm = T)
# other.count <- sum(y == 4, na.rm = T)
# cat("roadkill:", rk.count, "other death:", other.count, "\n")
# cat("roadkill / (roadkill + other death):", rk.count / (rk.count + other.count), "/n")
# 
# on.site <- sum(y == 1, na.rm = T)
# off.site <- sum(y == 2, na.rm = T)
# cat("on-site:", on.site, "off-site:", off.site, "/n")
# cat("off-site / (off-site + on-site):", off.site / (on.site + off.site), "/n")


## Model -----------------------------------------------------------------------

myCode <- nimbleCode({
  
  ## MISSING VALUES
  ## ---------------------------------------------------------------------------
  
  for (m in 1:nNoVeg){
    veg[noVeg[m]] ~ dnorm(0, sd = 1)
    # dens[noDens[m]] ~ dnorm(0, sd = 1)
  } # m
  
  # win[noWin] ~ dnorm(0, sd = 1)
  
  
  ## SURVIVAL & MOVEMENT MODELS
  ## ---------------------------------------------------------------------------
  
  for (a in 1:n.ageC){
    for (t in 1:(n.occasions-1)){
      
      # random year effect
      eps.phi[a, t] ~ dnorm(0, tau.phi)
      eps.M[a, t]   ~ dnorm(0, tau.M)
      eps.R[a, t]   ~ dnorm(0, tau.R)
      
      # logit-linear function of covariates
      logit(mean.phi[a, t]) <- logit(mu.phi[a]) + B.veg[a] * veg[t] + eps.phi[a, t]
      logit(mean.M[a, t])   <- logit(mu.M[a]) + eps.M[a, t]
      logit(mean.R[a, t])   <- logit(mu.R[a]) + eps.R[a, t]
      
    } # t
  } # a
  
  
  ## OBSERVATION & RECOVERY MODELS
  ## ---------------------------------------------------------------------------
  
  # for (t in 1:(n.occasions-1)){
  #   eps.rR[t] ~ dnorm(0, tau.rR)
  #   eps.rO[t] ~ dnorm(0, tau.rO)
  #   
  #   logit(mean.rR[t]) <- logit(mu.rR) + B.obsR * obs[t] + eps.rR[t]
  #   logit(mean.rO[t]) <- logit(mu.rO) + B.obsO * obs[t] + eps.rO[t]
  # } # t
  
  for (t in 1:n.occasions){
    eps.rR[t] ~ dnorm(0, tau.rR)
    eps.rO[t] ~ dnorm(0, tau.rO)
    
    logit(mean.rR[t]) <- logit(mu.rR) + B.obsR * obs[t] + eps.rR[t]
    logit(mean.rO[t]) <- logit(mu.rO) + B.obsO * obs[t] + eps.rO[t]
  } # t
  
  
  ## LIKELIHOOD
  ## ---------------------------------------------------------------------------

  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      
      # likelihood
      z[i, t] ~ dcat(trans.mat[i, z[i, t-1], 1:n.true.states, t-1]) # z is latent
      y[i, t] ~ dcat(obs.mat[i, z[i, t], 1:n.obs.states, t]) # y is observed
      
    } # t
  } # i
  
  
  ## TRANSITION MODEL
  ## ---------------------------------------------------------------------------
  
  for (i in 1:n.inds){
    for (t in first[i]:(n.occasions-1)){
      
      phi[i,t] <- mean.phi[ageC[age[i,t]], t] # survival
      R[i,t]   <- mean.R[ageC[age[i,t]], t]   # roadkill
      M[i,t]   <- mean.M[ageC[age[i,t]], t]   # migration
      
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
  
  
  ## OBSERVATION MODEL
  ## ---------------------------------------------------------------------------
  
  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      
      Pi[i,t] <- mean.Pi[t]  # observation on-site
      Po[i,t] <- mean.Po[t]  # observation off-site
      rR[i,t] <- mean.rR[t]  # recovery of roadkill
      rO[i,t] <- mean.rO[t]  # recovery of other death
      
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
      
    } # t
  } # i
  
  
  ## PRIORS
  ## ---------------------------------------------------------------------------
  
  # # quick simulations
  # hist(rbeta(1000, 8, 2))
  # hist(rbeta(1000, 4, 4))
  # hist(rbeta(1000, 2, 8))
  # hist(rbeta(1000, 1, 8))
  # hist(rbeta(1000, 1, 1))
  
  # mu.phi[1] ~ dbeta(4, 4)
  # mu.phi[2] ~ dbeta(8, 2)
  # mu.phi[3] ~ dbeta(8, 2)
  # mu.phi[4] ~ dbeta(8, 2)
  # mu.phi[5] ~ dbeta(4, 4)
  
  for (a in 1:n.ageC){
    mu.phi[a] ~ dbeta(4, 4)
    mu.M[a]   ~ dbeta(1, 8)
    mu.R[a]   ~ dbeta(1, 8)
    B.veg[a]  ~ dnorm(0, 1)
  } # a
  
  for (t in 1:n.occasions){
    mean.Pi[t] ~ dbeta(8, 1)
    mean.Po[t] ~ dbeta(4, 4)
  } # t
  
  mu.rR  ~ dbeta(4, 4)
  mu.rO  ~ dbeta(4, 4)
  B.obsR ~ dnorm(0, 1)
  B.obsO ~ dnorm(0, 1)
  
  sigma.phi ~ dunif(0, 4)
  sigma.M   ~ dunif(0, 4)
  sigma.R   ~ dunif(0, 4)
  sigma.rR  ~ dunif(0, 4)
  sigma.rO  ~ dunif(0, 4)
  
  tau.phi <- 1 / (sigma.phi * sigma.phi)
  tau.M   <- 1 / (sigma.M * sigma.M)
  tau.R   <- 1 / (sigma.R * sigma.R)
  tau.rR  <- 1 / (sigma.rR * sigma.rR)
  tau.rO  <- 1 / (sigma.rO * sigma.rO)
  
}) # nimbleCode


## Initial values --------------------------------------------------------------

# data is coded as OBSERVATION states
# initial values are provided for deterministic TRANSITIONS

# TRANSITION STATES
# 1 - alive on-site
# 2 - alive off-site
# 3 - dead by roadkill
# 4 - dead by other
# 5 - long dead

prepZs <- function(y){
  z_inits <- y
  z_dat   <- y
  
  z_inits[y == 999] <- NA
  z_dat[y == 999]   <- NA
  
  # Observed states -> deterministic transitions
  z_inits[y == 1] <- NA; z_dat[y == 1] <- 1  # alive on-site
  z_inits[y == 2] <- NA; z_dat[y == 2] <- 2  # alive off-site
  z_inits[y == 3] <- NA; z_dat[y == 3] <- 3  # dead by roadkill
  z_inits[y == 4] <- NA; z_dat[y == 4] <- 4  # dead by other cause
  z_inits[y == 5] <- 1 ; z_dat[y == 5] <- NA # undetected (alive)
  
  # Undetected after mortality -> long dead (5)
  n.inds <- nrow(y)
  for (i in 1:n.inds) {
    tmp <- z_dat[i, ]
    if (any(tmp == 3, na.rm = T) | any(tmp == 4, na.rm = T)) {
      indexD <- min(which(tmp == 3 | tmp == 4)) # index of 3 or 4 (dead)
      if(indexD < length(tmp)){ # safety in case roos die in last occasion
        indexLD <- (indexD + 1):length(tmp) # 5 after 3|4 -> 5 (long dead)
        tmp[indexLD] <- 5
      }
    }
    z_dat[i, ] <- tmp
  }
  z_inits[z_dat == 5] <- NA
  
  return(list(z_inits = z_inits,
              z_dat = z_dat))
}

ZZs <- prepZs(y)
z_inits <- ZZs$z_inits
z_dat <- ZZs$z_dat


## Assemble --------------------------------------------------------------------

# Inits
myInits <- list(
  z         = z_inits,
  mu.phi    = rbeta(n.ageC, 4, 4),
  mu.R      = rbeta(n.ageC, 1, 8),
  mu.M      = rbeta(n.ageC, 1, 8),
  mu.rR     = rbeta(1, 4, 4),
  mu.rO     = rbeta(1, 4, 4),
  mean.Pi   = rbeta(n.occasions, 8, 1),
  mean.Po   = rbeta(n.occasions, 4, 4),
  B.veg     = rnorm(n.ageC, 0, 0.1),
  B.obsR    = rnorm(1, 0, 0.1),
  B.obsO    = rnorm(1, 0, 0.1),
  eps.phi   = matrix(rnorm(n.ageC * (n.occasions-1), 0, 0.1), nrow = n.ageC, ncol = n.occasions-1),
  eps.M     = matrix(rnorm(n.ageC * (n.occasions-1), 0, 0.1), nrow = n.ageC, ncol = n.occasions-1),
  eps.R     = matrix(rnorm(n.ageC * (n.occasions-1), 0, 0.1), nrow = n.ageC, ncol = n.occasions-1),
  eps.rR    = rnorm(n.occasions, 0, 0.1),
  eps.rO    = rnorm(n.occasions, 0, 0.1),
  sigma.phi = runif(1, 0.5, 1.5),
  sigma.M   = runif(1, 0.5, 1.5),
  sigma.R   = runif(1, 0.5, 1.5),
  sigma.rR  = runif(1, 0.5, 1.5),
  sigma.rO  = runif(1, 0.5, 1.5)
)

# Data
y[y == 999] <- NA
myData <- list(y = y, 
               z = z_dat, 
               age = age,
               ageC = ageC,
               obs = obs,
               veg = veg)
               # dens = dens
               # win = win

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("B.veg", "B.obsR", "B.obsO",
            "mean.phi", "mean.M", "mean.R",
            "mu.phi", "mu.M", "mu.R", "mu.rR", "mu.rO",
            "mean.Pi", "mean.Po", "mean.rR", "mean.rO",
            "sigma.phi", "sigma.M", "sigma.R", "sigma.rR", "sigma.rO",
            "veg")

# Constants
myConst <- list(n.inds = n.inds,
                n.ageC = n.ageC,
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

# MCMC settings
if(testRun){
  nburn   <- 0            # burn-in
  niter   <- 10           # iterations
  nthin   <- 1            # thinning
  nchains <- 3            # chains
}else{
  nburn   <- 5000         # burn-in
  niter   <- 5000 + nburn # iterations
  nthin   <- 1            # thinning
  nchains <- 3            # chains
}


## Run model -------------------------------------------------------------------

if(parallelRun){
  # run one chain inside cluster
  runChain <- function(chainID, code, data, const, inits, params,
                       nburn, niter, nthin, seed){
    
    library(nimble)
    set.seed(seed)
    
    inits <- myInits
    model <- nimbleModel(code = myCode,
                         data = myData,
                         constants = myConst,
                         inits = myInits,
                         calculate = F,
                         check = F)
    
    conf <- configureMCMC(model, monitors = params, useConjugacy = F)
    mcmc <- buildMCMC(conf) # take a look at sampler
    Cmodel <- compileNimble(model, showCompilerOutput = F)
    Cmcmc <- compileNimble(mcmc, project = Cmodel)
    
    samples <- runMCMC(Cmcmc,
                       nburnin = nburn,
                       niter = niter,
                       thin = nthin,
                       inits = myInits,
                       setSeed = F,
                       progressBar = T,
                       samplesAsCodaMCMC = T)
    return(samples)
  }
  
  # create a cluster & export everything needed to each worker
  cl <- makeCluster(nchains)
  clusterExport(cl, varlist = c("myCode", "myData", "myConst", "myInits",
                                "params", "nburn", "niter", "nthin",
                                "runChain"))
}

if(parallelRun){
  # run chains in parallel
  start <- Sys.time()
  out <- parLapply(cl, 1:nchains, function(i){
    runChain(i,
             code = myCode,
             data = myData,
             const = myConst,
             inits = myInits,
             params = params,
             nburn = nburn,
             niter = niter,
             nthin = nthin,
             seed = i)})
  dur <- Sys.time() - start; dur
  stopCluster(cl)
  beep(2)
}else{
  # run chains sequentially
  start <- Sys.time()
  out <- nimbleMCMC(code = myCode,
                    data = myData,
                    constants = myConst,
                    inits = myInits,
                    monitors = params,
                    nburnin = nburn,
                    niter = niter,
                    thin = nthin,
                    nchains = nchains,
                    samplesAsCodaMCMC = T,
                    setSeed = 1:3)
  dur <- Sys.time() - start; dur
  beep(2)
}

# Save model, outputs, & diagnostics
MCMCdiag(out,
         dir = "./Results",
         save_object = T,
         obj_name = out.model,
         file_name = out.sum)

# # TEMP
# # which columns contain NAs
# library(coda)
# post <- as.mcmc(do.call(rbind, out))
# which(apply(post, 2, function(x) any(is.na(x) | is.nan(x))))


## Plots -----------------------------------------------------------------------

out.model <- "modelF_tObs_aVeg_atMR.rds"
out <- readRDS(paste0("results/", out.model))
model.summary <- MCMCsummary(object = out, round = 3)
model.summary

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
ageCs <- c("YAF", "1-2", "3-6", "7-9", "10+")

mcmc.df <- out %>% 
  map(~as.data.frame(as.matrix(.x))) %>% 
  bind_rows()

mcmc.df <- mcmc.df %>% 
  select(starts_with("mean.")) %>% 
  mutate(iter = row_number()) %>% 
  pivot_longer(cols = starts_with("mean."),
               names_to = "param.full",
               values_to = "value")

mcmc.df <- mcmc.df %>% 
  mutate(param = str_extract(param.full, "mean\\.[A-Za-z]+"),
         # extract all numbers inside brackets
         index = str_extract_all(param.full, "\\d+"),
         index1 = map_dbl(index, ~as.numeric(.x[1])),
         index2 = map_dbl(index, ~ifelse(length(.x) > 1, as.numeric(.x[2]), NA_real_)),
         
         # identify parameter dimensions
         is_time = param %in% c("mean.Pi", "mean.Po", "mean.rO", "mean.rR"),
         is_both = param %in% c("mean.M", "mean.R"),
         
         # assign t & a depending on parameter
         t = case_when(is_time ~ index1, is_both ~ index2),
         a = case_when(is_both ~ index1, TRUE ~ NA_real_)) %>% 
  select(iter, param.full, param, a, t, value)

summaries <- mcmc.df %>% 
  group_by(param, a, t) %>% 
  summarise(mean = mean(value, na.rm = TRUE),
            lcl = quantile(value, 0.025, na.rm = TRUE),
            ucl = quantile(value, 0.975, na.rm = TRUE),
            .groups = "drop") %>% 
  mutate(year = years[t],
         ageC = factor(ageCs[a], levels = ageCs))

# observation/recovery
summaries %>% 
  filter(param %in% c("mean.Pi","mean.Po","mean.rO","mean.rR")) %>% 
  ggplot(., aes(x = year, y = mean, fill = param, colour = param)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
  facet_wrap(~param, scales = "free_y") +
  labs(x = "Year", y = "Posterior mean (±95% CrI)") +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90", colour = NA))

# roadkill/migration
summaries %>% 
  filter(param %in% c("mean.M","mean.R")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  facet_wrap(~param, scales = "free_y") +
  labs(x = "Year", y = "Posterior mean (±95% CrI)", colour = "Age class", fill = "Age class") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", colour = NA))

# rs7.age <- pred.df %>%
#   ggplot(., aes(x = .data[[var]], y = rs.mean)) +
#   geom_ribbon(aes(ymin = rs.lower, ymax = rs.upper, fill = prs, group = prs), alpha = 0.2, show.legend = F) +
#   geom_line(aes(colour = prs, group = prs), linewidth = 1, show.legend = F) +
#   geom_jitter(aes(x = .data[[var]], y = y, colour = prs),
#               inherit.aes = FALSE,
#               data = obs.data %>% filter(!is.na(prs)),
#               height = 0.02,
#               alpha = 0.4,
#               show.legend = F) +
#   scale_colour_manual(values = prs.colours) +
#   scale_fill_manual(values = prs.colours) +
#   # xlim(this.range) +
#   scale_x_continuous(breaks = c(2,6,10,14,18)) +
#   theme_bw() +
#   theme(axis.title = element_text(size = 12),
#         axis.title.x = element_blank()) +
#   labs(x = this.name,
#        y = "Probability of producing a LPY",
#        fill = "Previous\nreproductive\nsuccess",
#        colour = "Previous\nreproductive\nsuccess"); rs7.age

# Checks
library(coda)
c1 <- as.mcmc(out$chain1)

MCMCplot(object = out,
         horiz = FALSE,
         rank = TRUE,
         ref_ovl = FALSE)

MCMCtrace(object = out,
          pdf = TRUE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE, # add Rhat diagnostics
          n.eff = TRUE, # add eff sample size
          params = c("B.veg", # "B.obsR", "B.obsO",
                     "mean.phi", "mean.M", "mean.R",
                     "mu.phi", "mu.M", "mu.R", "mu.rR", "mu.rO",
                     "mean.Pi", "mean.Po", "mean.rR", "mean.rO",
                     "sigma.phi", "sigma.M", "sigma.R",
                     "sigma.rR", "sigma.rO"))

# Correlation
autocorr.diag(out)
# autocorr.plot(out)
coda::crosscorr.plot(out)

