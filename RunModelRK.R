# 9 September 2023

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
source("PrepDataRK.R")

y <- eh[, -1] %>% as.matrix() %>% unname()
n_obs_states <- 4
n_true_states <- 4


## Model -----------------------------------------------------------------------

code <- nimbleCode({
  #### Priors & constraints ####
  # TODO more informative priors
  # beta prior preferable to uniform because survival is very high
  # explore beta shape parameters using rbeta(10000, x, y)
  
  # Survival probability
  phi ~ dbeta(1, 1) # roughly equivalent to dunif(0, 1)
  mR  ~ dbeta(1, 1) # mortality by roadkill
  
  # Observations
  p  ~ dbeta(1, 1)  # probability of observation
  rR ~ dbeta(1, 1)  # probability of recovery, roadkill
  rO ~ dbeta(1, 1)  # probability of recovery, natural death
  # rO <- rR        # probability of recovery, natural death
  
  #### Transition matrix ####
  # 1 - alive
  # 2 - dead by roadkill
  # 3 - dead by other
  # 4 - long dead
  
  # ALIVE
  trans.mat[1, 1] <- phi
  trans.mat[2, 1] <- 0                                                  
  trans.mat[3, 1] <- 0                                                 
  trans.mat[4, 1] <- 0
  
  # DEAD BY ROADKILL
  trans.mat[1, 2] <- (1 - phi) * mR
  trans.mat[2, 2] <- 0
  trans.mat[3, 2] <- 0
  trans.mat[4, 2] <- 0
  
  # DEAD BY OTHER
  trans.mat[1, 3] <- (1 - phi) * (1 - mR)
  trans.mat[2, 3] <- 0
  trans.mat[3, 3] <- 0
  trans.mat[4, 3] <- 0
  
  # LONG DEAD
  trans.mat[1, 4] <- 0
  trans.mat[2, 4] <- 1
  trans.mat[3, 4] <- 1
  trans.mat[4, 4] <- 1
  
  #### Observation matrix ####
  # 1 - seen
  # 2 - recovered roadkill
  # 3 - recovered other
  # 4 - undetected
  
  # SEEN
  obs.mat[1, 1] <- p
  obs.mat[2, 1] <- 0
  obs.mat[3, 1] <- 0
  obs.mat[4, 1] <- 0
  
  # RECOVERED ROADKILL
  obs.mat[1, 2] <- 0
  obs.mat[2, 2] <- rR
  obs.mat[3, 2] <- 0
  obs.mat[4, 2] <- 0
  
  # RECOVERED ROADKILL
  obs.mat[1, 3] <- 0
  obs.mat[2, 3] <- 0
  obs.mat[3, 3] <- rO
  obs.mat[4, 3] <- 0
  
  # UNDETECTED
  obs.mat[1, 4] <- 1 - p
  obs.mat[2, 4] <- 1 - rR
  obs.mat[3, 4] <- 1 - rO
  obs.mat[4, 4] <- 1
  
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
# 1 - alive
# 2 - dead by roadkill
# 3 - dead by other
# 4 - long dead

z_inits <- y
z_dat <- y
z_inits[y == 999] <- NA
z_dat[y == 999] <- NA

# seen -> alive (1 -> 1)
z_inits[y == 1] <- NA
z_dat[y == 1] <- 1

# recovered roadkill -> dead by roadkill (2 -> 2)
z_inits[y == 2] <- NA
z_dat[y == 2] <- 2

# recovered other -> dead by other (3 -> 3)
z_inits[y == 3] <- NA
z_dat[y == 3] <- 3

# undetected -> alive (4 -> 1)? TO DISCUSS
z_inits[y == 4] <- 1
z_dat[y == 4] <- NA

# undetected after mortality observed -> long dead (4 -> 4)
for (i in 1:n_inds) {
  tmp <- z_dat[i, ]
  if (any(tmp == 2, na.rm = T) | any(tmp == 3, na.rm = T)) { # if undetected after mortality observed
    # index of 2 or 3 (dead)
    indexD <- which(tmp == 2 | tmp == 3)
    # indices of 4s after 2|3 (dead) -> replace with 4 (LD)
    indexLD <- (indexD + 1)[indexD + 1 <= length(tmp)]  # indexLD <- (indexD + 1):length(tmp)
    tmp[indexLD] <- 4
  }
  z_dat[i, ] <- tmp
}

z_inits[z_dat == 4] <- NA


## Assemble --------------------------------------------------------------------
# Inits
inits <- list(
  z = z_inits,
  phi = rbeta(1, 1, 1),
  mR = rbeta(1, 1, 1),
  p = rbeta(1, 1, 1))

# Data
y[y == 999] <- NA
dat <- list(y = y, z = z_dat)

# Parameters to monitor
# best practice is to only include things that are directly sampled (i.e. have a prior)
# anything derived can be done post-hoc, unless you want the model to give annual survival
# when debugging, could add trans.mat & obs.mat, or even z, etc.

params <- c("phi", "mR", "p", "rR", "rO")

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
beep(sound = 1)

# Run MCMC
t.start <- Sys.time()
sink("Results/errorreport.txt") # for debugging
out <- runMCMC(Cmcmc, niter = ni, nburnin = nb, nchains = nc, inits = inits,
               setSeed = F, progressBar = T, samplesAsCodaMCMC = T)
t.end <- Sys.time()
sink() # closing txt file
(runTime <- t.end - t.start)
beep(sound = 1)


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
          params = c("phi", "mR"))

ex <- MCMCchains(out)

# Correlation plots
autocorr.diag(out)
autocorr.plot(out)
coda::crosscorr.plot(out)

posterior <- as.array(out)

# Scatterplot of variables
color_scheme_set("purple")
phi_mR <- mcmc_scatter(out, pars = c("phi", "mR"))
out_pairsplot <- mcmc_pairs(out, diag_fun = "dens", off_diag_args = list(size = 0.5))

