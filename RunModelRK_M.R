# 28 August 2024
# Original known-fate model
# with an extra state for alive but off the study area

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

y <- eh[, -1] %>% as.matrix() %>% unname()
# n_occasions <- 16
n_obs_states <- 5
n_true_states <- 5


## Model -----------------------------------------------------------------------

code <- nimbleCode({
  #### Priors & constraints ####
  # TODO more informative priors
  # beta prior preferable to uniform because survival is very high
  # explore beta shape parameters using rbeta(10000, x, y)
  
  # Survival probability
  phi ~ dbeta(8, 2)  # roughly equivalent to dunif(0, 1)
  R   ~ dbeta(1, 4)  # mortality by roadkill
  M   ~ dbeta(1, 4)  # migration (in or out)
  
  # Observations
  Pi  ~ dbeta(8, 2)  # probability of observation, on-site
  Po  ~ dbeta(1, 1)  # probability of observation, off-site
  rR  ~ dbeta(2, 2)  # probability of recovery, roadkill
  rO  ~ dbeta(2, 2)  # probability of recovery, natural death
  # rO <- rR
  
  #### Transition matrix ####
  # 1 - alive on-site
  # 2 - alive off-site
  # 3 - dead by roadkill
  # 4 - dead by other
  # 5 - long dead
  
  # ALIVE ON-SITE
  trans.mat[1, 1] <- phi*(1-M)
  trans.mat[2, 1] <- phi*M
  trans.mat[3, 1] <- 0
  trans.mat[4, 1] <- 0
  trans.mat[5, 1] <- 0
  
  # ALIVE OFF-SITE
  trans.mat[1, 2] <- phi*M
  trans.mat[2, 2] <- phi*(1-M)
  trans.mat[3, 2] <- 0
  trans.mat[4, 2] <- 0
  trans.mat[5, 2] <- 0
  
  # DEAD BY ROADKILL
  trans.mat[1, 3] <- (1-phi)*R
  trans.mat[2, 3] <- (1-phi)*R
  trans.mat[3, 3] <- 0
  trans.mat[4, 3] <- 0
  trans.mat[5, 3] <- 0
  
  # DEAD BY OTHER
  trans.mat[1, 4] <- (1-phi)*(1-R)
  trans.mat[2, 4] <- (1-phi)*(1-R)
  trans.mat[3, 4] <- 0
  trans.mat[4, 4] <- 0
  trans.mat[5, 4] <- 0
  
  # LONG DEAD
  trans.mat[1, 5] <- 0
  trans.mat[2, 5] <- 0
  trans.mat[3, 5] <- 1
  trans.mat[4, 5] <- 1
  trans.mat[5, 5] <- 1
  
  #### Observation matrix ####
  # 1 - seen on-site
  # 2 - seen off-site
  # 3 - recovered roadkill
  # 4 - recovered other
  # 5 - undetected
  
  # SEEN ON-SITE
  obs.mat[1, 1] <- Pi
  obs.mat[2, 1] <- 0
  obs.mat[3, 1] <- 0
  obs.mat[4, 1] <- 0
  obs.mat[5, 1] <- 0
  
  # SEEN OFF-SITE
  obs.mat[1, 2] <- 0
  obs.mat[2, 2] <- Po
  obs.mat[3, 2] <- 0
  obs.mat[4, 2] <- 0
  obs.mat[5, 2] <- 0
  
  # RECOVERED ROADKILL
  obs.mat[1, 3] <- 0
  obs.mat[2, 3] <- 0
  obs.mat[3, 3] <- rR
  obs.mat[4, 3] <- 0
  obs.mat[5, 3] <- 0
  
  # RECOVERED OTHER
  obs.mat[1, 4] <- 0
  obs.mat[2, 4] <- 0
  obs.mat[3, 4] <- 0
  obs.mat[4, 4] <- rO
  obs.mat[5, 4] <- 0
  
  # UNDETECTED
  obs.mat[1, 5] <- 1-Pi
  obs.mat[2, 5] <- 1-Po
  obs.mat[3, 5] <- 1-rR
  obs.mat[4, 5] <- 1-rO
  obs.mat[5, 5] <- 1
  
  # END priors & constraints
  
  #### Likelihood ####
  for (i in 1:n_inds) {
    for (t in (first[i] + 1):n_occasions) { # conditioning on first capture
      z[i, t] ~ dcat(trans.mat[z[i, t-1], 1:n_true_states]) # z is latent
      y[i, t] ~ dcat(obs.mat[z[i, t], 1:n_obs_states]) # y is observed
    } # t
  } # i
  
  # END likelihood
})


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
for (i in 1:n_inds) {
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

# # Laurie (to test loop)
# # undetected after mortality observed -> long dead (4 -> 4)
# for (i in 1:n_inds) {
#   i = 55
#   tmp <- z_dat[i, ]
#   tmp
#   print(i)
#   if (any(tmp == 2, na.rm = T) | any(tmp == 3, na.rm = T)) { # if undetected after mortality observed
#     # index of 2 or 3 (dead)
#     indexD <- which(tmp == 2 | tmp == 3)
#     # indices of 4s after 2|3 (dead) -> replace with 4 (LD)
#     indexLD <- (indexD + 1)[indexD + 1 <= length(tmp)]  # indexLD <- (indexD + 1):length(tmp)
#     tmp[indexLD] <- 4
#   }
#   tmp
#   z_dat[i, ] <- tmp
# }
# 
# z_inits[z_dat == 4] <- NA


## Assemble --------------------------------------------------------------------
# Inits
inits <- list(
  z   = z_inits,
  phi = rbeta(1, 1, 1),
  R   = rbeta(1, 1, 1),
  M   = rbeta(1, 1, 1),
  Pi  = rbeta(1, 1, 1),
  Po  = rbeta(1, 1, 1))

# Data
y[y == 999] <- NA
dat <- list(y = y, z = z_dat)

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("phi", "R", "M", "Pi", "Po", "rR", "rO")

# Constants
const <- list(
  n_inds = n_inds,
  n_occasions = dim(dat$y)[2],
  n_true_states = n_true_states,
  n_obs_states = n_obs_states,
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
          params = c("phi", "R", "M"))

ex <- MCMCchains(out)

# Correlation plots
autocorr.diag(out)
autocorr.plot(out)
coda::crosscorr.plot(out)

posterior <- as.array(out)

# Scatterplot of variables
color_scheme_set("purple")
phi_R <- mcmc_scatter(out, pars = c("phi", "R"))
out_pairsplot <- mcmc_pairs(out, diag_fun = "dens", off_diag_args = list(size = 0.5))

