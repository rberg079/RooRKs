# 29 January 2026
# Road mortality analyses


## Set up ----------------------------------------------------------------------

# set toggles
females <- TRUE
testRun <- FALSE
parallelRun <- TRUE

# name outputs
out.model <- "modelF_tObs_tR_Zinits3.rds"
out.sum <- "modelF_tObs_tR_Zinits3_sum.txt"

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
source("PrepDataRK.R")
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
#                          mean.p = NULL,
#                          mean.rR = NULL,
#                          mean.rO = NULL,
#                          seed = 123)
# list2env(dataRK, envir = .GlobalEnv)


## Raw data checks -------------------------------------------------------------

# table(y)
# 
# rk.count <- sum(y == 3, na.rm = T)
# other.count <- sum(y == 4, na.rm = T)
# cat("roadkill:", rk.count, "other death:", other.count)
# cat("roadkill / (roadkill + other death):", rk.count / (rk.count + other.count))
# 
# on.site <- sum(y == 1, na.rm = T)
# off.site <- sum(y == 2, na.rm = T)
# cat("on-site:", on.site, "off-site:", off.site)
# cat("off-site / (off-site + on-site):", off.site / (on.site + off.site))


## Model -----------------------------------------------------------------------

myCode <- nimbleCode({
  
  ## SURVIVAL & MOVEMENT MODELS
  ## ---------------------------------------------------------------------------
  
  for (a in 1:n.ageC){
    for (t in 1:(n.occasions-1)){
      
      # random year effect
      eps.phi[a, t] ~ dnorm(0, tau.phi)
      eps.R[t]      ~ dnorm(0, tau.R)
      
      # logit-linear functions
      logit(mean.phi[a, t]) <- logit(mu.phi[a]) + eps.phi[a, t]
      logit(mean.R[a, t])   <- logit(mu.R[a]) + eps.R[t]
      
    } # t
  } # a
  
  
  ## OBSERVATION & RECOVERY MODELS
  ## ---------------------------------------------------------------------------
  
  for (t in 1:n.occasions){
    
    # logit-linear functions
    logit(mean.p[t]) <- logit(mu.p)
    logit(mean.rR[t]) <- logit(mu.rR)
    logit(mean.rO[t]) <- logit(mu.rO)
    
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
      
      #### Transition matrix ####
      # 1 - alive
      # 2 - dead by roadkill
      # 3 - dead by other
      # 4 - long dead
      
      # ALIVE
      trans.mat[i,1,1,t] <- phi[i,t]
      trans.mat[i,2,1,t] <- 0
      trans.mat[i,3,1,t] <- 0
      trans.mat[i,4,1,t] <- 0
      
      # DEAD BY ROADKILL
      trans.mat[i,1,2,t] <- (1-phi[i,t])*R[i,t]
      trans.mat[i,2,2,t] <- 0
      trans.mat[i,3,2,t] <- 0
      trans.mat[i,4,2,t] <- 0
      
      # DEAD BY OTHER
      trans.mat[i,1,3,t] <- (1-phi[i,t])*(1-R[i,t])
      trans.mat[i,2,3,t] <- 0
      trans.mat[i,3,3,t] <- 0
      trans.mat[i,4,3,t] <- 0
      
      # LONG DEAD
      trans.mat[i,1,4,t] <- 0
      trans.mat[i,2,4,t] <- 1
      trans.mat[i,3,4,t] <- 1
      trans.mat[i,4,4,t] <- 1
      
    } # t
  } # i
  
  
  ## OBSERVATION MODEL
  ## ---------------------------------------------------------------------------
  
  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      
      p[i,t]  <- mean.p[t]   # observation
      rR[i,t] <- mean.rR[t]  # recovery
      rO[i,t] <- mean.rO[t]  # recovery
      
      #### Observation matrix ####
      # 1 - seen on-site
      # 2 - recovered roadkill
      # 3 - recovered other
      # 4 - undetected
      
      # SEEN
      obs.mat[i,1,1,t] <- p[i,t]
      obs.mat[i,2,1,t] <- 0
      obs.mat[i,3,1,t] <- 0
      obs.mat[i,4,1,t] <- 0
      
      # RECOVERED ROADKILL
      obs.mat[i,1,2,t] <- 0
      obs.mat[i,2,2,t] <- rR[i,t]
      obs.mat[i,3,2,t] <- 0
      obs.mat[i,4,2,t] <- 0
      
      # RECOVERED OTHER
      obs.mat[i,1,3,t] <- 0
      obs.mat[i,2,3,t] <- 0
      obs.mat[i,3,3,t] <- rO[i,t]
      obs.mat[i,4,3,t] <- 0
      
      # UNDETECTED
      obs.mat[i,1,4,t] <- 1-p[i,t]
      obs.mat[i,2,4,t] <- 1-rR[i,t]
      obs.mat[i,3,4,t] <- 1-rO[i,t]
      obs.mat[i,4,4,t] <- 1
      
    } # t
  } # i
  
  
  ## PRIORS
  ## ---------------------------------------------------------------------------
  
  # # quick simulations
  # hist(rbeta(1000, 8, 1))
  # hist(rbeta(1000, 8, 2))
  # hist(rbeta(1000, 8, 4))
  # hist(rbeta(1000, 4, 4))
  # hist(rbeta(1000, 2, 8))
  # hist(rbeta(1000, 1, 8))
  
  for (a in 1:n.ageC){
    mu.R[a] ~ dbeta(2, 8)
  } # a
  
  # informative priors on survival
  # based on CJS models in Ecology paper
  mu.phi[1] ~ dbeta(8, 2)   # 1 year-old subadults
  mu.phi[2] ~ dbeta(8, 2)   # 2 year-old subadults
  mu.phi[3] ~ dbeta(8, 2)   # prime-aged adults
  mu.phi[4] ~ dbeta(8, 2)   # pre-senescent
  mu.phi[5] ~ dbeta(8, 4)   # senescent
  
  # Pi known to be extremely high &
  # to vary little from Ecology paper
  mu.p ~ dbeta(20, 2)
  mu.rR ~ dbeta(4, 4)
  mu.rO ~ dbeta(4, 4)
  
  sigma.phi ~ dexp(10)
  sigma.R   ~ dexp(10)
  
  tau.phi <- 1 / (sigma.phi * sigma.phi)
  tau.R   <- 1 / (sigma.R * sigma.R)
  
}) # nimbleCode


## Initial values --------------------------------------------------------------

# data is coded as OBSERVATION states
# initial values are provided for deterministic TRANSITIONS

# LATENT STATES
# 1 - alive
# 2 - dead by roadkill
# 3 - dead by other
# 4 - long dead

# OBSERVATION CODES
# 1 - seen alive
# 2 - recovered roadkill
# 3 - recovered other
# 4 - undetected

# # version 1
# prepZs <- function(y){
#   z_inits <- y
#   z_dat   <- y
# 
#   z_inits[y == 999] <- NA
#   z_dat[y == 999]   <- NA
# 
#   # Observed states -> deterministic transitions
#   z_inits[y == 1] <- NA; z_dat[y == 1] <- 1  # alive
#   z_inits[y == 2] <- NA; z_dat[y == 2] <- 2  # dead by roadkill
#   z_inits[y == 3] <- NA; z_dat[y == 3] <- 3  # dead by other cause
#   z_inits[y == 4] <- 1 ; z_dat[y == 4] <- NA # undetected
# 
#   # Undetected after mortality -> long dead (4)
#   n.inds <- nrow(y)
#   for (i in 1:n.inds){
#     tmp <- z_dat[i, ]
#     if (any(tmp == 2, na.rm = T) | any(tmp == 3, na.rm = T)) {
#       indexD <- min(which(tmp == 2 | tmp == 3)) # index of 2 or 3 (dead)
#       if(indexD < length(tmp)){ # safety in case roos die in last occasion
#         indexLD <- (indexD + 1):length(tmp) # 4 after 2|3 -> 4 (long dead)
#         tmp[indexLD] <- 4
#       }
#     }
#     z_dat[i, ] <- tmp
#   }
#   z_inits[z_dat == 4] <- NA
# 
#   return(list(z_inits = z_inits,
#               z_dat = z_dat))
# }
# 
# ZZs <- prepZs(y)
# z_inits <- ZZs$z_inits
# z_dat <- ZZs$z_dat

# # version 2
# # fills in occasions between first & last, when known alive
# # proporses random time of death sometime after disappearance
# prepZs <- function(y){
#   n.inds <- nrow(y)
#   n.occ <- ncol(y)
# 
#   z_inits <- matrix(NA, n.inds, n.occ)
#   z_dat   <- matrix(NA, n.inds, n.occ)
#   
#   for(i in 1:n.inds){
#     # find first and last observation
#     obs_idx <- which(y[i,] != 999 & y[i,] != 4) # 4 is undetected
#     
#     if(length(obs_idx) > 0){
#       f <- min(obs_idx)
#       l <- max(obs_idx)
#       
#       # DATA: known alive from first to last seen
#       # Note: z_dat is whatever y is at l (last)
#       z_dat[i, f:l] <- y[i, f:l]
#       
#       # INITS: fill the gap between last seen & end
#       if(l < n.occ){
#         # if last was a death (2 or 3), fill with 4
#         if(y[i, l] %in% c(2, 3)){
#           z_dat[i, (l + 1):n.occ] <- 4 # known long dead
#         }else{
#           # if last was alive (1), propose a time of death
#           # pick a random time of death to help mixing
#           maybe_die <- sample((l + 1):(n.occ + 1), 1)
#           if(maybe_die <= n.occ){
#             z_inits[i, (l + 1):maybe_die] <- 1 # alive until death
#             # pick a random cause for inits
#             z_inits[i, maybe_die] <- sample(2:3, 1)
#             if(maybe_die < n.occ){
#               z_inits[i, (maybe_die + 1):n.occ] <- 4
#             }
#           }else{
#             z_inits[i, (l + 1):n.occ] <- 1 # stays alive
#           }
#         }
#       }
#     }
#   }
#   return(list(z_inits = z_inits,
#               z_dat = z_dat))
# }
# 
# ZZs <- prepZs(y)
# z_inits <- ZZs$z_inits
# z_dat <- ZZs$z_dat
# 
# # # safeties
# # z_dat[z_dat == 999] <- NA
# # z_dat[z_dat == 4] <- NA # undetected is NOT a known state
# # # add deterministic post-death states calculated in prepZs
# # z_dat[!is.na(ZZs$z_dat) & ZZs$z_dat == 4] <- 4

# version 3
# chatGPT's suggested middle ground
prepZs <- function(y){
  
  n.inds <- nrow(y)
  n.occ  <- ncol(y)
  
  z_dat   <- matrix(NA, n.inds, n.occ)  # DATA: logically fixed states
  z_inits <- matrix(NA, n.inds, n.occ)  # INITS: diffuse guesses only
  
  for(i in 1:n.inds){
    
    # 1. map known latent states
    z_dat[i, y[i, ] == 1] <- 1  # seen alive
    z_dat[i, y[i, ] == 2] <- 2  # recovered roadkill
    z_dat[i, y[i, ] == 3] <- 3  # recovered other
    # undetected (4) and 999 remain NA
    
    # 2. enforce long-dead after known death
    death_idx <- which(z_dat[i, ] %in% c(2, 3))
    if(length(death_idx) > 0){
      d <- min(death_idx)
      if(d < n.occ){
        z_dat[i, (d + 1):n.occ] <- 4
      }
    }
    
    # 3. diffuse initial values for unknowns
    for(t in 1:n.occ){
      if(is.na(z_dat[i, t])){
        if(length(death_idx) == 0 || t < min(death_idx)){ # before any known death
          z_inits[i, t] <- sample(c(1, 2, 3), 1, prob = c(0.2, 0.4, 0.4)) # sample state
        }else{
          z_inits[i, t] <- 4 # long dead after known death
        }
      }
    }
  }
  
  return(list(
    z_dat   = z_dat,
    z_inits = z_inits
  ))
}


ZZs <- prepZs(y)
z_dat   <- ZZs$z_dat
z_inits <- ZZs$z_inits


## Assemble --------------------------------------------------------------------

# Inits
myInits <- list(
  z         = z_inits,
  mu.phi    = rbeta(n.ageC, 8, 4),
  mu.R      = rbeta(n.ageC, 2, 8),
  mu.p      = rbeta(1, 8, 1),
  mu.rR     = rbeta(1, 4, 4),
  mu.rO     = rbeta(1, 4, 4),
  eps.phi   = matrix(rnorm(n.ageC * (n.occasions-1), 0, 0.1), nrow = n.ageC, ncol = n.occasions-1),
  eps.R     = rnorm(n.occasions, 0, 0.1),
  sigma.phi = rexp(1, 10),
  sigma.R   = rexp(1, 10)
)

# Data
y[y == 999] <- NA
myData <- list(y = y, 
               z = z_dat, 
               age = age,
               ageC = ageC)

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("mu.phi", "mu.R", "mu.p", "mu.rR", "mu.rO",
            "mean.phi", "mean.R", "mean.p", "mean.rR", "mean.rO",
            "sigma.phi", "sigma.R")

# Constants
myConst <- list(n.inds = n.inds,
                n.ageC = n.ageC,
                n.occasions = n.occasions,
                n.true.states = n.true.states,
                n.obs.states = n.obs.states,
                first = first)
# noVeg = noVeg
# nNoVeg = nNoVeg

# # Check that z[, first] is known for all inds...
# for (ii in 1:n.inds) {
#   print(z_dat[ii, first[ii]])
# }

# MCMC settings
if(testRun){
  nburn   <- 0             # burn-in
  niter   <- 10            # iterations
  nthin   <- 1             # thinning
  nchains <- 3             # chains
}else{
  nburn   <- 10000         # burn-in
  niter   <- 10000 + nburn # iterations
  nthin   <- 1             # thinning
  nchains <- 3             # chains
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

# # IF NEEDED:
# # which columns contain NAs?
# library(coda)
# post <- as.mcmc(do.call(rbind, out))
# which(apply(post, 2, function(x) any(is.na(x) | is.nan(x))))


## Plots -----------------------------------------------------------------------

# out.model <- "modelF_tObs_tR_newZinits.rds"
# out <- readRDS(paste0("results/", out.model))
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
ageCs <- c("1", "2", "3-6", "7-9", "10+")

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
         is_time = param %in% c("mean.p", "mean.rR", "mean.rO"),
         is_both = param %in% c("mean.phi", "mean.R"),
         
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
  filter(param %in% c("mean.p","mean.rR","mean.rO")) %>% 
  ggplot(., aes(x = year, y = mean, fill = param, colour = param)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl), alpha = 0.2, colour = NA) +
  facet_wrap(~param, scales = "free_y") +
  labs(x = "Year", y = "Posterior mean (±95% CrI)") +
  ylim(0, 1) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey90", colour = NA))

# roadkill/migration
summaries %>% 
  filter(param %in% c("mean.R")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  # facet_wrap(~param, scales = "free_y") +
  labs(x = "Year", y = "Posterior mean (±95% CrI)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", colour = NA))

# survival
summaries %>% 
  filter(param %in% c("mean.phi")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1) +
  # facet_wrap(~param, scales = "free_y") +
  labs(x = "Year", y = "Posterior mean (±95% CrI)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", colour = NA))

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
          params = c("mu.phi", "mu.R", "mu.p", "mu.rR", "mu.rO",
                     # "mean.phi", "mean.R", "mean.p", "mean.r",
                     "sigma.phi", "sigma.R"))

# Correlation
autocorr.diag(out)
# autocorr.plot(out)
coda::crosscorr.plot(out)


## Compare model outputs -------------------------------------------------------

source('compareModels.R')
CompareModels(postPaths = c(
  "results/modelF_tObs_atMR_oner.rds",
  "results/modelF_tObs_atMR_noPiRE.rds",
  "results/modelF_tObs_atMR_noREs.rds",
  "results/modelF_tObs_atMR_tREsMR.rds",
  "results/modelF_tObs_atMR_tREsPhi.rds"
),
modelNames = c(
  "modF_oner",
  "modF_noPiRE",
  "modF_noREs",
  "modF_tREsMR",
  "modF_tREsPhi"
),
plotFolder = c("figures/11.noAgedREs"),
returnSumData = TRUE)

