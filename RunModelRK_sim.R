# 29 January 2026
# Road mortality analyses


## Set up ----------------------------------------------------------------------

# set toggles
females <- FALSE
testRun <- FALSE
parallelRun <- TRUE

# name outputs
out.model <- "modelF_tObs_aVR_itX_tR_noRec_4ageCs.rds"
out.sum <- "modelF_tObs_aVR_itX_tR_noRec_4ageCs_sum.txt"

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

# load data
source("PrepDataRK_sim.R")
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
# rk.count <- sum(y == 2, na.rm = T)
# other.count <- sum(y == 3, na.rm = T)
# cat("roadkill:", rk.count, "other death:", other.count)
# cat("roadkill / (roadkill + other death):", rk.count / (rk.count + other.count))


## Model -----------------------------------------------------------------------

myCode <- nimbleCode({
  
  ## MISSING VALUES
  ## ---------------------------------------------------------------------------
  
  for (m in 1:nNoVegR){
    vegr[noVegR[m]] ~ dnorm(0, 1)
  } # m

  for (i in 1:n.inds) {
    for (t in 1:n.occasions) {
      xmed[i, t] ~ dnorm(0, sd = 1)
    }
  }
  
  
  ## SURVIVAL & ROAD MORTALITY MODELS
  ## ---------------------------------------------------------------------------
  
  # random year effect
  for (t in 1:(n.occasions-1)){
    eps.R[t] ~ dnorm(0, tau.R)
  }
  
  for (a in 1:n.ageC){
    for (t in 1:(n.occasions-1)){
      
      # random year effect
      eps.phi[a, t] ~ dnorm(0, tau.phi)
      
      # logit-linear functions
      logit(mean.phi[a, t]) <- logit(mu.phi[a]) +
        betaVR.phi[a] * vegr[t] +
        eps.phi[a, t]
      
      logit(mean.R[a, t]) <- logit(mu.R[a]) + 
        betaVR.R[a] * vegr[t] +
        eps.R[t]
      
    } # t
  } # a
  
  
  ## OBSERVATION MODEL
  ## ---------------------------------------------------------------------------
  
  for (t in 2:n.occasions){
    
    # random year effect
    eps.p[t] ~ dnorm(0, tau.p)
    
    # logit-linear functions
    logit(mean.p[t]) <- logit(mu.p) + eps.p[t]
    
  } # t
  
  mean.p[1] <- mu.p
  eps.p[1]  <- 0
  
  
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
      
      logit(phi[i,t]) <- logit(mean.phi[ageC[age[i,t]], t]) + betaX.phi * xmed[i,t] # survival
      logit(R[i,t])   <- logit(mean.R[ageC[age[i,t]], t]) + betaX.R * xmed[i,t]     # roadkill
      
      #### Transition matrix ####
      # 1 - alive
      # 2 - new roadkill
      # 3 - new other death
      # 4 - long dead (absorbing)
      
      # ALIVE
      trans.mat[i,1,1,t] <- phi[i,t]
      trans.mat[i,1,2,t] <- (1-phi[i,t]) * R[i,t]     # died by road
      trans.mat[i,1,3,t] <- (1-phi[i,t]) * (1-R[i,t]) # died by other
      trans.mat[i,1,4,t] <- 0
      
      # NEW ROADKILL
      trans.mat[i,2,1,t] <- 0
      trans.mat[i,2,2,t] <- 0
      trans.mat[i,2,3,t] <- 0
      trans.mat[i,2,4,t] <- 1
      
      # NEW OTHER DEATH
      trans.mat[i,3,1,t] <- 0
      trans.mat[i,3,2,t] <- 0
      trans.mat[i,3,3,t] <- 0
      trans.mat[i,3,4,t] <- 1
      
      # LONG DEAD (ABSORBING)
      trans.mat[i,4,1,t] <- 0
      trans.mat[i,4,2,t] <- 0
      trans.mat[i,4,3,t] <- 0
      trans.mat[i,4,4,t] <- 1
      
    } # t
  } # i
  
  
  ## OBSERVATION MODEL
  ## ---------------------------------------------------------------------------
  
  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      
      p[i,t]  <- mean.p[t] # observation
      
      # #### Observation matrix ####
      # # LOW ROADKILL SCENARIO:
      # # 1 - seen
      # # 2 - recovered roadkill
      # # 3 - undetected
      # 
      # # ALIVE
      # obs.mat[i,1,1,t] <- p[i,t]
      # obs.mat[i,1,2,t] <- 0
      # obs.mat[i,1,3,t] <- 1-p[i,t]
      # 
      # # NEW ROADKILL
      # obs.mat[i,2,1,t] <- 0
      # obs.mat[i,2,2,t] <- 1
      # obs.mat[i,2,3,t] <- 0
      # 
      # # NEW OTHER DEATH
      # obs.mat[i,3,1,t] <- 0
      # obs.mat[i,3,2,t] <- 0
      # obs.mat[i,3,3,t] <- 1
      # 
      # # LONG DEAD (ABSORBING)
      # obs.mat[i,4,1,t] <- 0
      # obs.mat[i,4,2,t] <- 0
      # obs.mat[i,4,3,t] <- 1
      # # assuming other deaths are not recovered
      # # but rather disappear & remain undetected
      
      #### Observation matrix ####
      # HIGH ROADKILL SCENARIO:
      # 1 - seen
      # 2 - undetected

      # ALIVE
      obs.mat[i,1,1,t] <- p[i,t]
      obs.mat[i,1,2,t] <- 1-p[i,t]

      # NEW ROADKILL
      obs.mat[i,2,1,t] <- 0
      obs.mat[i,2,2,t] <- 1

      # NEW OTHER DEATH
      obs.mat[i,3,1,t] <- 0
      obs.mat[i,3,2,t] <- 1

      # LONG DEAD (ABSORBING)
      obs.mat[i,4,1,t] <- 0
      obs.mat[i,4,2,t] <- 1
      # assuming no one is recovered
      # all rather remain undetected
      
    } # t
  } # i
  
  
  ## PRIORS
  ## ---------------------------------------------------------------------------
  
  # # quick simulations
  # # according to CJS models in Ecology paper:
  # hist(rnorm(1000, 0.89, 0.029), xlim = c(0,1)) # phi age 1-2
  # hist(rnorm(1000, 0.95, 0.014), xlim = c(0,1)) # phi age 3-6
  # hist(rnorm(1000, 0.90, 0.024), xlim = c(0,1)) # phi age 7-9
  # hist(rnorm(1000, 0.70, 0.092), xlim = c(0,1)) # phi age 10+
  # 
  # hist(rbeta(1000, 12, 12), xlim = c(0,1)) # phi age 0
  # hist(rbeta(1000, 20, 4), xlim = c(0,1))  # phi age 1-2
  # hist(rbeta(1000, 20, 2), xlim = c(0,1))  # phi age 3-6
  # hist(rbeta(1000, 20, 4), xlim = c(0,1))  # phi age 7-9
  # hist(rbeta(1000, 20, 12), xlim = c(0,1)) # phi age 10+
  # 
  # hist(rbeta(1000, 12, 16), xlim = c(0,1)) # phi age 0
  # hist(rbeta(1000, 20, 6), xlim = c(0,1))  # phi age 1-2
  # hist(rbeta(1000, 20, 4), xlim = c(0,1))  # phi age 3-6
  # hist(rbeta(1000, 20, 8), xlim = c(0,1))  # phi age 7-9
  # hist(rbeta(1000, 20, 14), xlim = c(0,1)) # phi age 10+
  
  for (a in 1:n.ageC){
    mu.R[a] ~ dbeta(2, 8)
    betaVR.phi[a] ~ dnorm(0, 1)
    betaVR.R[a]   ~ dnorm(0, 1)
  } # a
  
  betaX.phi ~ dnorm(0, 1)
  betaX.R   ~ dnorm(0, 1)
  
  # informative priors on survival
  # based on CJS models in Ecology paper
  # females:
  mu.phi[1] ~ dbeta(20, 4)  # subadults
  mu.phi[2] ~ dbeta(20, 2)  # prime-aged adults
  mu.phi[3] ~ dbeta(20, 4)  # pre-senescent
  mu.phi[4] ~ dbeta(20, 12) # senescent
  
  # # males:
  # mu.phi[1] ~ dbeta(20, 6)  # subadults
  # mu.phi[2] ~ dbeta(20, 4)  # prime-aged adults
  # mu.phi[3] ~ dbeta(20, 8)  # pre-senescent
  # mu.phi[4] ~ dbeta(20, 14) # senescent
  
  # Pi known to be extremely high &
  # to vary little from Ecology paper
  mu.p  ~ dbeta(40, 2) # females
  # mu.p  ~ dbeta(20, 6) # males
  
  sigma.phi ~ dexp(1)  # was 10, 1 helped w convergence
  sigma.R   ~ dexp(1)  # was 10, 1 helped w convergence
  sigma.p   ~ dexp(10) # was 10, 1 helped w convergence
  
  tau.phi <- 1 / (sigma.phi * sigma.phi)
  tau.R   <- 1 / (sigma.R * sigma.R)
  tau.p   <- 1 / (sigma.p * sigma.p)
  
}) # nimbleCode


## Initial values --------------------------------------------------------------

# data is coded as OBSERVATION states
# initial values are provided for deterministic TRANSITIONS

# LATENT STATES
# 1 - alive
# 2 - new roadkill
# 3 - new other death
# 4 - long dead (absorbing)

# OBSERVATION CODES (LOW RK)
# 1 - seen alive
# 2 - recovered roadkill
# 3 - undetected

# OBSERVATION CODES (HIGH RK)
# 1 - seen alive
# 2 - undetected

y[y == 999] <- NA

prepZs <- function(y, first, last, id){
  
  zData  <- y
  zInits <- y
  
  # map observed events to known latent states
  zInits[y == 1] <- NA; zData[y == 1] <- 1  # alive
  zInits[y == 2] <- NA; zData[y == 2] <- 2  # roadkill
  zInits[y == 3] <- NA; zData[y == 3] <- 3  # other death
  zInits[y == 4] <- NA; zData[y == 4] <- NA # undetected
  
  for(i in 1:nrow(y)){
    f <- first[i]
    l <- last[i]
    
    fate <- id$fate[i]
    gone <- id$gone[i]
    
    # undetected between f & l
    # must be alive but missed
    if(f != l){
      zInits[i, f:(l-1)] <- 1
    }
    
    if(!is.na(fate)){
      # disappeared roos w known fates
      # could only be undetected once gone
      if(gone <= ncol(y)) zData[i, gone:ncol(y)] <- 4 
    }else{
      # disappeared roos w unknown fates
      # likely died then went undetected
      if(l < ncol(y)){
        zInits[i, l + 1] <- 3             # first new death
        if((l + 2) <= ncol(y)){
          zInits[i, (l + 2):ncol(y)] <- 4 # then long dead
        }
      }
      
      # # ...or were alive but missed
      # if(l < ncol(y)){
      #   zInits[i, (l + 1):ncol(y)] <- 1
      # }
    }
  }
  
  zInits[!is.na(zData)] <- NA
  
  return(list(zInits = zInits,
              zData  = zData))
}

ZZs <- prepZs(y, first, last, id)
zInits <- ZZs$zInits
zData <- ZZs$zData

# HIGH RK
y[y == 2] <- 2
y[y == 3] <- 2
y[y == 4] <- 2

# # LOW RK
# y[y == 4] <- 3


## Assemble --------------------------------------------------------------------

# Inits
myInits <- list(
  z          = zInits,
  mu.phi     = rbeta(n.ageC, 4, 2),
  mu.R       = rbeta(n.ageC, 2, 8),
  mu.p       = rbeta(1, 20, 4),
  betaVR.phi = rnorm(n.ageC, 0, 0.1),
  betaX.phi  = rnorm(1, 0, 0.1),
  betaVR.R   = rnorm(n.ageC, 0, 0.1),
  betaX.R    = rnorm(1, 0, 0.1),
  eps.phi    = matrix(rnorm(n.ageC * (n.occasions-1), 0, 0.1), nrow = n.ageC, ncol = n.occasions-1),
  eps.R      = rnorm(n.occasions-1, 0, 0.1),
  eps.p      = rnorm(n.occasions, 0, 0.1),
  sigma.phi  = rexp(1, 10),
  sigma.R    = rexp(1, 10),
  sigma.p    = rexp(1, 10)
)

# Data
myData <- list(y = y, 
               z = zData, 
               age = age,
               ageC = ageC,
               vegr = vegr,
               xmed = xmed)

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c(
  "betaVR.phi", "betaVR.R",
  "betaX.phi", "betaX.R",
  "mu.phi", "mu.R", "mu.p",
  "mean.phi", "mean.R", "mean.p",
  "eps.phi", "eps.R", "eps.p",
  "sigma.phi", "sigma.R", "sigma.p",
  "vegr"
)

# Constants
myConst <- list(n.inds = n.inds,
                n.ageC = n.ageC,
                n.occasions = n.occasions,
                n.true.states = n.true.states,
                n.obs.states = n.obs.states,
                first = first,
                noVegR = noVegR,
                nNoVegR = nNoVegR,
                noX = noX,
                nNoX = nNoX)

# # Check that z[, first] is known for all inds...
# for (ii in 1:n.inds) {
#   print(zData[ii, first[ii]])
# }

# MCMC settings
if(testRun){
  nburn   <- 0             # burn-in
  niter   <- 10            # iterations
  nthin   <- 1             # thinning
  nchains <- 3             # chains
}else{
  nburn   <- 70000         # burn-in
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


## Plot results ----------------------------------------------------------------

# out.model <- "modelM_tObs_aVD_itX_tR_noRecO_wYAFs.rds"
# out <- readRDS(paste0("results/", out.model))

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
summ.means %>% 
  filter(param %in% c("mean.R")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Road mortality (lowest)", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        strip.background = element_rect(fill = "grey90", colour = NA))

# survival vs year
summ.means %>% 
  filter(param %in% c("mean.phi")) %>% 
  ggplot(., aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lcl, ymax = ucl, fill = ageC, group = ageC), alpha = 0.2, show.legend = F) +
  geom_line(aes(colour = ageC, group = ageC), linewidth = 1, show.legend = F) +
  scale_x_continuous(breaks = c(2008, 2012, 2016, 2020, 2024)) +
  labs(x = "Year", y = "Survival", colour = "Age class", fill = "Age class") +
  ylim(0, 1) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90", colour = NA))


## Compare model outputs -------------------------------------------------------

# source('compareModels.R')
# CompareModels(postPaths = c(
#   "results/modelF_tObs_aV_tR_noRecO_wYAFs.rds",
#   "results/modelF_tObs_aV&D_tR_noRecO_wYAFs.rds",
#   "results/modelF_tObs_aVR_tR_noRecO_wYAFs.rds"
# ),
# modelNames = c(
#   "modF_V",
#   "modF_V&D",
#   "modF_VR"
# ),
# plotFolder = c("figures/15.envCovars"),
# returnSumData = TRUE,
# nAgeC = 6)

