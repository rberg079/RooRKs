# 5 Septembre 2024
# Attempt to run same model
# with dDHMMo funciton instead

library(readxl)
library(tidyverse)
library(lubridate)


## Set up ----------------------------------------------------------------------

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
library(nimbleEcology)

# load data
source("PrepDataRK_M.R")

# # or...
# eh <- read_csv("eh.csv")
# env <- read_csv("env.csv")

y <- eh[,-1] %>% as.matrix() %>% unname()

n.occasions <- 16
n.obs.states <- 5
n.true.states <- 5

# missing veg


## Model -----------------------------------------------------------------------

code <- nimbleCode({
  
  
  # # model age-class-specific survival probability
  # # logit-linear function of time-dependent covariates
  # # priors for mu.age defined further in the model 
  # for (j in 1:(n.occasions-1)){
  #   logit(phi.juv[j]) <- mu.juv # + B.veg[1]*veg[j] + eps.phi[j]
  #   logit(phi.sub[j]) <- mu.sub # + B.veg[2]*veg[j] + eps.phi[j]
  #   logit(phi.pri[j]) <- mu.pri # + B.veg[3]*veg[j] + eps.phi[j]
  #   logit(phi.pre[j]) <- mu.pre # + B.veg[4]*veg[j] + eps.phi[j]
  #   logit(phi.sen[j]) <- mu.sen # + B.veg[5]*veg[j] + eps.phi[j]
  #   mean.phi[1, j] <- phi.juv[j]
  #   mean.phi[2, j] <- phi.sub[j]
  #   mean.phi[3, j] <- phi.pri[j]
  #   mean.phi[4, j] <- phi.pre[j]
  #   mean.phi[5, j] <- phi.sen[j]
  # }
  
  for (i in 1:n.inds){
    
    y[i, first[i]:n.occasions] ~ dDHMMo(init = inits[i, 1:5],                                                                  # initial state probabilities, provide this as data (1s and 0s, assuming first state known)
                                        probObs = obs.mat[i, 1:n.true.states, 1:n.obs.states, first[i]:n.occasions],           # time dependent 3d array (index by i) [nstates, nevents, t]; assume known at first and provide as data
                                        probTrans = trans.mat[i, 1:n.true.states, 1:n.true.states, (first[i]+1):n.occasions],  # time dependent 3d array (index by i) [nstates, nstates, t]
                                        len = n.occasions - first[i] + 1,                                                      # length of observations
                                        checkRowSums = 0) 
    
    for (t in (first[i]+1):n.occasions){
      
      #### Likelihood ####
      # z[i, t] ~ dcat(trans.mat[i, t, z[i, t-1], 1:n.true.states]) # z is latent
      # y[i, t] ~ dcat(obs.mat[i, t, z[i, t], 1:n.obs.states]) # y is observed
      
      
      #### Priors & constraints ####
      # beta prior preferable to uniform because survival is very high
      # explore beta shape parameters using rbeta(10000, x, y)
      
      # Survival
      logit(phi[i,t]) <- mean.phi[ageC[age[i,t-1]], t-1]  # + B.veg[ageC[age[i,t]]]*veg[t] + eps.phi[t] or eps.phi[t,ageC[age[i,t]]] ?
      
      
      # Fixed
      R[i,t]   <- mean.R   # mortality by roadkill
      M[i,t]   <- mean.M   # migration (in or out)
      
      Pi[i,t]  <- mean.Pi  # probability of observation, on-site
      Po[i,t]  <- mean.Po  # probability of observation, off-site
      rR[i,t]  <- mean.rR  # probability of recovery, roadkill
      rO[i,t]  <- mean.rO  # probability of recovery, natural death
      
      
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
  
  # mu.juv   ~ dbeta(8, 2)
  # mu.sub   ~ dbeta(8, 2)
  # mu.pri   ~ dbeta(8, 2)
  # mu.pre   ~ dbeta(8, 2)
  # mu.sen   ~ dbeta(8, 2)
  
  for (a in 1:n.age){
    for (t in (first[i]+1):n.occasions){
      mean.phi[a,t] ~ dbeta(8, 2)
    }
  }
  
  mean.R   ~ dbeta(1, 4)
  mean.M   ~ dbeta(1, 4)
  
  mean.Pi  ~ dbeta(8, 2)
  mean.Po  ~ dbeta(8, 2)
  mean.rR  ~ dbeta(2, 2)
  mean.rO  ~ dbeta(2, 2)
  
  
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

# z_inits <- y
# z_dat <- y
# z_inits[y == 999] <- NA
# z_dat[y == 999] <- NA
# 
# # seen on-site -> alive on-site (1 -> 1)
# z_inits[y == 1] <- NA
# z_dat[y == 1] <- 1
# 
# # seen off-site -> alive off-site (2 -> 2)
# z_inits[y == 2] <- NA
# z_dat[y == 2] <- 2
# 
# # recovered roadkill -> dead by roadkill (3 -> 3)
# z_inits[y == 3] <- NA
# z_dat[y == 3] <- 3
# 
# # recovered other -> dead by other (4 -> 4)
# z_inits[y == 4] <- NA
# z_dat[y == 4] <- 4
# 
# # undetected -> alive (5 -> 1)? TO DISCUSS
# z_inits[y == 5] <- 1
# z_dat[y == 5] <- NA
# 
# # undetected after mortality observed -> long dead (5 -> 5)
# for (i in 1:n.inds) {
#   tmp <- z_dat[i, ]
#   if (any(tmp == 3, na.rm = T) | any(tmp == 4, na.rm = T)) { # if undetected after mortality observed
#     # index of 3 or 4 (dead)
#     indexD <- which(tmp == 3 | tmp == 4)
#     # indices of 5s after 3|4 (dead) -> replace with 5 (LD)
#     indexLD <- (indexD + 1):length(tmp) # indexLD <- (indexD + 1)[indexD + 1 <= length(tmp)]
#     tmp[indexLD] <- 5
#   }
#   z_dat[i, ] <- tmp
# }
# 
# z_inits[z_dat == 5] <- NA


## Assemble --------------------------------------------------------------------
# Inits
# inits <- list(
#   z        = z_inits,
#   # mu.juv   = rbeta(1, 8, 2),
#   # mu.sub   = rbeta(1, 8, 2),
#   # mu.pri   = rbeta(1, 8, 2),
#   # mu.pre   = rbeta(1, 8, 2),
#   # mu.sen   = rbeta(1, 8, 2),
#   # mean.phi = rbeta(n.age, 8, 2),
#   mean.R   = rbeta(1, 1, 4),
#   mean.M   = rbeta(1, 1, 4),
#   mean.Pi  = rbeta(1, 8, 2),
#   mean.Po  = rbeta(1, 8, 2),
#   mean.rR  = rbeta(1, 2, 2),
#   mean.rO  = rbeta(1, 2, 2))

# figure out how to initialize for dDHMMo!!

# Data
y[y == 999] <- NA
dat <- list(y = y, age = age, ageC = ageC) # z = z_dat

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("mean.phi",
            # "mu.juv", "mu.sub", "mu.pri", "mu.pre", "mu.sen",
            "mean.R", "mean.M", "mean.Pi", "mean.Po", "mean.rR", "mean.rO")

# Constants
const <- list(
  n.inds = n.inds,
  n.occasions = n.occasions,
  n.true.states = n.true.states,
  n.obs.states = n.obs.states,
  first = first)


## Run model -------------------------------------------------------------------
# MCMC settings
nb <- 5000       # burn-in
ni <- 5000 + nb  # total iterations
nt <- 1          # thin
nc <- 3          # chains

# Compile configuration & build
Rmodel <- nimbleModel(code = code, constants = const, data = dat,
                      check = F, calculate = F, inits = inits)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt) # slow
Rmcmc <- buildMCMC(conf)

Cmodel <- compileNimble(Rmodel, showCompilerOutput = F)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
beep(sound = 2)

# Run MCMC
t.start <- Sys.time()
sink("Results/errorreport.txt") # for debugging
out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits = inits,
               setSeed = F, progressBar = T, samplesAsCodaMCMC = T)
t.end <- Sys.time()
sink() # closing txt file
(runTime <- t.end - t.start)
beep(sound = 2)

