# 2 February 2025
# Attempt to add age-specific survival rates
# as well as time-varying covariates (veg for now)

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

# load data
source("PrepDataRK_M.R")

# # or...
# eh <- read_csv("eh.csv")
# env <- read_csv("env.csv")

y <- eh[,-1] %>% as.matrix() %>% unname()

n.occasions <- ncol(y) # 16
n.obs.states <- 5
n.true.states <- 5

# missing veg

# AEB edit
unkAge <- which(apply(age, 1, function(x){all(is.na(x))}))

y <- y[-unkAge, ]
age <- age[-unkAge, ]
first <- apply(y, 1, get.first)
last  <- apply(y, 1, get.last)

# # quick simulations
# hist(plogis(rnorm(1000, 1, 0.75)))
# 
# hist(rbeta(1000, .8*10, (1-.8)*1))
# hist(rbeta(1000, 8, 2))


## Model -----------------------------------------------------------------------

code <- nimbleCode({
  
  # model age-class-specific survival probability
  # logit-linear function of time-dependent covariates
  # priors for mu.age defined further in the model 
  for (j in 1:(n.occasions-1)){
    # logit(phi.juv[j]) <- mu.juv # + B.veg[1]*veg[j] + eps.phi[j]
    # logit(phi.sub[j]) <- mu.sub # + B.veg[2]*veg[j] + eps.phi[j]
    # logit(phi.pri[j]) <- mu.pri # + B.veg[3]*veg[j] + eps.phi[j]
    # logit(phi.pre[j]) <- mu.pre # + B.veg[4]*veg[j] + eps.phi[j]
    # logit(phi.sen[j]) <- mu.sen # + B.veg[5]*veg[j] + eps.phi[j]
    
    # AEB note - changed here to test
    phi.juv[j] <- mu.juv # + B.veg[1]*veg[j] + eps.phi[j]
    phi.sub[j] <- mu.sub # + B.veg[2]*veg[j] + eps.phi[j]
    phi.pri[j] <- mu.pri # + B.veg[3]*veg[j] + eps.phi[j]
    phi.pre[j] <- mu.pre # + B.veg[4]*veg[j] + eps.phi[j]
    phi.sen[j] <- mu.sen # + B.veg[5]*veg[j] + eps.phi[j]
    
    mean.phi[1, j] <- phi.juv[j]
    mean.phi[2, j] <- phi.sub[j]
    mean.phi[3, j] <- phi.pri[j]
    mean.phi[4, j] <- phi.pre[j]
    mean.phi[5, j] <- phi.sen[j]
  }
  
  for (i in 1:n.inds){
    for (t in (first[i]+1):n.occasions){
      # TO DISCUSS: should it not be (t in (first[i]):(n.occasions-1))?
      # TO DISCUSS: should n.occasions not be last[i]?
      
      #### Likelihood ####
      z[i, t] ~ dcat(trans.mat[i, z[i, t-1], 1:n.true.states, t]) # z is latent
      y[i, t] ~ dcat(obs.mat[i, z[i, t], 1:n.obs.states, t]) # y is observed
      
      
      #### Priors & constraints ####
      # beta prior preferable to uniform because survival is very high
      # explore beta shape parameters using rbeta(10000, x, y)
      
      # Survival
      # AEB note - do you want index t on left and t-1 on right? SOLVED
      phi[i,t] <- mean.phi[ageC[age[i,t]], t-1]
      
      # Fixed
      R[i,t]   <- mean.R   # mortality by roadkill
      M[i,t]   <- mean.M   # migration (in or out)
      
      # AEB - not sure why this is being estimated with a posterior predictive sampler??
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
  
  # mu.juv   ~ dbeta(8, 2) # AEB note - these priors would be on the probability scale
  # mu.sub   ~ dbeta(8, 2) # AEB note - think you want the real scale if you want to add covariates
  # mu.pri   ~ dbeta(8, 2) 
  # mu.pre   ~ dbeta(8, 2)
  # mu.sen   ~ dbeta(8, 2)
  
  mu.juv   ~ dunif(-2, 2)
  mu.sub   ~ dunif(-2, 2)
  mu.pri   ~ dunif(-2, 2) 
  mu.pre   ~ dunif(-2, 2)
  mu.sen   ~ dunif(-2, 2)
  
  mean.R   ~ dbeta(2, 8) 
  mean.M   ~ dbeta(2, 8) 
  
  mean.Pi  ~ dbeta(8, 2)
  mean.Po  ~ dbeta(8, 2) # AEB note - do you expect your prob of observation off-site to be as high as on-site?
  mean.rR  ~ dbeta(2, 2) # AEB note - looks good, although you might want to add some effort metric if this varies by season
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
n.inds <- nrow(y) # AEB note - added here
for (i in 1:n.inds) {
  tmp <- z_dat[i, ]
  if (any(tmp == 3, na.rm = T) | any(tmp == 4, na.rm = T)) { # if undetected after mortality observed
    # index of 3 or 4 (dead)
    indexD <- which(tmp == 3 | tmp == 4)
    # indices of 5s after 3|4 (dead) -> replace with 5 (LD)
    indexLD <- (indexD + 1):length(tmp) # indexLD <- (indexD + 1)[indexD + 1 <= length(tmp)]
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
  z        = z_inits,
  mu.juv   = runif(1, -2, 2),
  mu.sub   = runif(1, -2, 2),
  mu.pri   = runif(1, -2, 2),
  mu.pre   = runif(1, -2, 2),
  mu.sen   = runif(1, -2, 2),
  mean.R   = rbeta(1, 2, 8),
  mean.M   = rbeta(1, 2, 8),
  mean.Pi  = rbeta(1, 8, 2),
  mean.Po  = rbeta(1, 8, 2),
  mean.rR  = rbeta(1, 2, 2),
  mean.rO  = rbeta(1, 2, 2))

# Data
y[y == 999] <- NA
dat <- list(y = y, 
            z = z_dat, 
            age = age, ageC = ageC)

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("mu.juv", "mu.sub", "mu.pri", "mu.pre", "mu.sen",
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
         obj_name = "model1.rds", file_name = "model_summary.txt")
# mod <- readRDS(here("Results", "model1.rds")
# mod_summary <- read.lines(here("Results", "model_summary.txt"))

# Calculate numerical summaries
model_summary <- MCMCsummary(object = out, round = 3)

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

ex <- MCMCchains(out)

# Correlation plots
autocorr.diag(out)
autocorr.plot(out)
coda::crosscorr.plot(out)

posterior <- as.array(out)

# Scatterplot of variables
color_scheme_set("purple")
mcmc_scatter(out, pars = c("mu.pri", "mean.R"))
# mcmc_pairs(out, diag_fun = "dens", off_diag_args = list(size = 0.5))
