############################################################################## #
# Bayesian known-fate model for deer GPS collar 
# WEEKLY models of survival
# Author: Rachel Hickcox, Abby Bratt (Proteus)
# Client: DRNSW
# Year: 2023
############################################################################## #

#### LOAD LIBRARIES ####
install_and_load_packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) { install.packages(new.pkg, dependencies = TRUE)}
  lapply(pkg, library, character.only = TRUE)
}
install_and_load_packages(c("tidyverse", # for data processing; includes ggplot
                            "forcats", # specific to factors
                            "nimble", # to run the model
                            "here", # easy file paths across users
                            "coda", # to compute diagnostics
                            "beepr", # make computer beep when code is done running
                            "cowplot", # arrange multiple plots
                            "patchwork", # arrange multiple plots
                            "RColorBrewer", # color palettes 
                            "bayesplot", # plotting results from Bayesian MCMC using ggplot
                            "tidybayes", # plotting results from Bayesian MCMC using ggplot
                            "ggdist", # plotting results from Bayesian MCMC using ggplot
                            "postpack", # handy functions for post-processing of Bayesian MCMC results
                            "strex", # helpful string functions
                            "MCMCvis")) # making beautiful plots

#### SOURCE FILES #####
source(here("Code", "01_ProcessData.R"))

#### DATA INPUT #####
# Weekly model of survival
y <- eh_weekly[, -1] %>% as.matrix() %>% unname() 
n_obs_states <- 5
n_true_states <- 8

############################################################################# #
### WEEKLY MODEL ##############################################################
############################################################################# #

code <- nimbleCode({
  ##### Priors and constraints ####
  # TODO more informative priors. survival probably very high (i.e., over 0.9)
  # beta prior preferable to uniform for this reason
  # can explore beta shape parameters using rbeta(10000, x, y)
  
  # SURVIVAL PROBABILITY
  phi ~ dbeta(1, 1) # roughly equivalent to dunif(0,1)
  mort.h ~ dbeta(1, 1) # human cause
  
  # OBSERVATIONS
  active ~ dbeta(1, 1)
  p ~ dbeta(1, 1) # detection prob alive + active
  q ~ dbeta(1, 1) # detection prob collar failure
  r_O ~ dbeta(1,1) # detection prob dead other
  r_H <- r_O # detection prob dead humans = dead other
  
  ##### Transition matrix ####
  # 1 - alive and active collar
  # 2 - alive and inactive collar, observed
  # 3 - alive and inactive collar, undetected
  # 4 - newly caused natural mortality + active collar
  # 5 - newly caused natural mortality + inactive collar
  # 6 - newly caused human mortality + active collar
  # 7 - newly caused human mortality + inactive collar
  # 8 - long dead
  
  # ALIVE AND ACTIVE (AA)
  trans.mat[1, 1] <- phi*active                
  trans.mat[2, 1] <- 0                                                    
  trans.mat[3, 1] <- 0                                                 
  trans.mat[4, 1] <- 0 
  trans.mat[5, 1] <- 0
  trans.mat[6, 1] <- 0
  trans.mat[7, 1] <- 0
  trans.mat[8, 1] <- 0
  
  # ALIVE AND INACTIVE, OBSERVED (AI_O)
  trans.mat[1, 2] <- phi * (1 - active)
  trans.mat[2, 2] <- 0 
  trans.mat[3, 2] <- 0 
  trans.mat[4, 2] <- 0 
  trans.mat[5, 2] <- 0 
  trans.mat[6, 2] <- 0 
  trans.mat[7, 2] <- 0 
  trans.mat[8, 2] <- 0 
  
  # ALIVE AND INACTIVE, UNDETECTED (AI_U)
  trans.mat[1, 3] <- 0
  trans.mat[2, 3] <- phi 
  trans.mat[3, 3] <- phi 
  trans.mat[4, 3] <- 0 
  trans.mat[5, 3] <- 0 
  trans.mat[6, 3] <- 0 
  trans.mat[7, 3] <- 0 
  trans.mat[8, 3] <- 0 
  
  # NEWLY DEAD-OTHER AND ACTIVE (DO_A)
  trans.mat[1, 4] <- (1-phi)*(1-mort.h)*active
  trans.mat[2, 4] <- 0 
  trans.mat[3, 4] <- 0 
  trans.mat[4, 4] <- 0 
  trans.mat[5, 4] <- 0 
  trans.mat[6, 4] <- 0 
  trans.mat[7, 4] <- 0 
  trans.mat[8, 4] <- 0 
  
  # NEWLY DEAD-OTHER AND INACTIVE (DO_I)
  trans.mat[1, 5] <- (1-phi)*(1-active)*(1-mort.h)
  trans.mat[2, 5] <- (1-phi)*(1-mort.h)
  trans.mat[3, 5] <- (1-phi)*(1-mort.h)
  trans.mat[4, 5] <- 0 
  trans.mat[5, 5] <- 0 
  trans.mat[6, 5] <- 0 
  trans.mat[7, 5] <- 0 
  trans.mat[8, 5] <- 0 
  
  # NEWLY DEAD-HUMAN AND ACTIVE (DH_A)
  trans.mat[1, 6] <- (1-phi)*mort.h*active
  trans.mat[2, 6] <- 0
  trans.mat[3, 6] <- 0
  trans.mat[4, 6] <- 0
  trans.mat[5, 6] <- 0
  trans.mat[6, 6] <- 0
  trans.mat[7, 6] <- 0
  trans.mat[8, 6] <- 0
  
  # NEWLY DEAD-HUMAN AND INACTIVE (DH_I)
  trans.mat[1, 7] <- (1-phi)*mort.h*(1-active)
  trans.mat[2, 7] <- (1-phi)*mort.h
  trans.mat[3, 7] <- (1-phi)*mort.h
  trans.mat[4, 7] <- 0
  trans.mat[5, 7] <- 0
  trans.mat[6, 7] <- 0
  trans.mat[7, 7] <- 0
  trans.mat[8, 7] <- 0
  
  # LONG DEAD (LD)
  trans.mat[1, 8] <- 0
  trans.mat[2, 8] <- 0
  trans.mat[3, 8] <- 0
  trans.mat[4, 8] <- 1
  trans.mat[5, 8] <- 1
  trans.mat[6, 8] <- 1
  trans.mat[7, 8] <- 1
  trans.mat[8, 8] <- 1
  
  ##### Observation matrix ####
  # 1 - alive and active collar
  # 2 - failed collar, observed
  # 3 - newly dead other mortality + active collar
  # 4 - newly dead human mortality + active collar
  # 5 - undetected
  
  # ALIVE AND ACTIVE
  obs.mat[1, 1] <- p              
  obs.mat[2, 1] <- 0
  obs.mat[3, 1] <- 0                                             
  obs.mat[4, 1] <- 0
  obs.mat[5, 1] <- 0
  obs.mat[6, 1] <- 0
  obs.mat[7, 1] <- 0
  obs.mat[8, 1] <- 0
  
  # Fail
  obs.mat[1, 2] <- 0
  obs.mat[2, 2] <- q
  obs.mat[3, 2] <- 0
  obs.mat[4, 2] <- 0
  obs.mat[5, 2] <- 0
  obs.mat[6, 2] <- 0
  obs.mat[7, 2] <- 0
  obs.mat[8, 2] <- 0
  
  # NEWLY DEAD-OTHER 
  obs.mat[1, 3] <- 0
  obs.mat[2, 3] <- 0 
  obs.mat[3, 3] <- 0 
  obs.mat[4, 3] <- r_O
  obs.mat[5, 3] <- 0 
  obs.mat[6, 3] <- 0 
  obs.mat[7, 3] <- 0 
  obs.mat[8, 3] <- 0 
  
  # NEWLY DEAD-HUMAN 
  obs.mat[1, 4] <- 0
  obs.mat[2, 4] <- 0  
  obs.mat[3, 4] <- 0  
  obs.mat[4, 4] <- 0
  obs.mat[5, 4] <- 0
  obs.mat[6, 4] <- r_H
  obs.mat[7, 4] <- 0
  obs.mat[8, 4] <- 0
  
  # NOT DETECTED
  obs.mat[1, 5] <- 1-p
  obs.mat[2, 5] <- 1-q 
  obs.mat[3, 5] <- 1 
  obs.mat[4, 5] <- 1-r_O
  obs.mat[5, 5] <- 1
  obs.mat[6, 5] <- 1-r_H
  obs.mat[7, 5] <- 1
  obs.mat[8, 5] <- 1
  
  # END priors and constraints
  
  ##### Likelihood #####
  
  for (i in 1:n_ind) {
    for (t in (first[i] + 1):n_occasions) { # note conditioning on first capture here
      z[i, t] ~ dcat(trans.mat[z[i, t-1], 1:n_true_states]) # z is latent (true state); categorical distribution
      y[i, t] ~ dcat(obs.mat[z[i, t], 1:n_obs_states]) # y is observed data; categorical distribution
    } # t last
  } # i
  
  # END likelihood
})

#### INIITAL VALUES ####
# Data is coded as OBSERVATION states. Providing model with initial values for deterministic TRANSITION states

# TRANSITION STATES
# 1 - alive and active collar
# 2 - alive and inactive collar, observed
# 3 - alive and inactive collar, undetected
# 4 - newly caused natural mortality + active collar
# 5 - newly caused natural mortality + inactive collar
# 6 - newly caused human mortality + active collar
# 7 - newly caused human mortality + inactive collar
# 8 - long dead

##### Defining initial states #####
z_inits <- y
z_dat <- y
z_inits[y == 999] <- NA
z_dat[y == 999] <- NA

# active collar -> alive (1 -> 1)
z_inits[y == 1] <- NA
z_dat[y == 1] <- 1
# collar failed -> alive and observed inactive collar (2 -> 2)
z_inits[y == 2] <- NA
z_dat[y == 2] <- 2
# newly caused natural mortality -> newly caused natural mortality (3 -> 4)
z_inits[y == 3] <- NA
z_dat[y == 3] <- 4
# newly caused human mortality -> newly caused human mortality (4 -> 6)
z_inits[y == 4] <- NA
z_dat[y == 4] <- 6
# undetected -> alive and unobserved inactive collar (5 -> 1)
z_inits[y == 5] <- 1
z_dat[y == 5] <- NA

# If long inactive after mortality signal observed -> long dead (5 -> 8)
for (i in 1:n_inds) {
  tmp <- z_dat[i, ]
  if (any(tmp == 4, na.rm = TRUE) | any(tmp == 6, na.rm = TRUE)) { # if long inactive after mortality signal observed
    # Index of 3 of 4 (dead)
    indexD <- which(tmp == 4 | tmp == 6)
    # Indices of 5's after 4|6 (dead) -> replace with 8 (LD)
    indexLD <- (indexD + 1):length(tmp) 
    tmp[indexLD] <- 8
  }
  z_dat[i,] <- tmp
}

z_inits[z_dat == 8] <- NA

# If collar failure (confirmed inactive) -> all values after are inactive but undetected (2 -> 3)
for (i in 1:n_inds) {
  tmp <- z_inits[i,]
  if (any(y[i,] == 2, na.rm = TRUE)) {  
    # Index of 2 inactive confirmed
    indexD <- which(y[i,] == 2)
    # Indices after 2 inactive -> replace with 3
    indexLD <- (indexD + 1):length(tmp) 
    tmp[indexLD] <- 3
  } 
  z_inits[i,] <- tmp
}

z_dat[z_inits == 3] <- NA

##### Inits #####
inits <- list(
  z = z_inits,
  phi = rbeta(1, 1, 1), 
  mort.h = rbeta(1, 1, 1), 
  active = rbeta(1, 1, 1),
  p = rbeta(1, 1, 1), # detection prob alive + active
  q = rbeta(1, 1, 1)) # detection prob collar failure

#### DATA ####
y[y == 999] <- NA
dat <- list(
  y = y,
  z = z_dat) 

#### PARAMETERS TO MONITOR ####
# put anything that you want to track in here
# the more parameters, the slower the model
# best practice to only include things that are directly sampled (i.e., you have a prior for them)
# anything derived can be done post-hoc, unless you want the model to spit out annual survival probability
# when debugging model, could add trans.mat and obs.mat, or even z, etc. 
params <- c(
  "phi", 
  "mort.h",
  "active",
  "p",  "q", 
  "r_O", "r_H")

#### CONSTANTS ####
const <- list(
  n_ind = n_inds, 
  n_occasions = dim(dat$y)[2], 
  n_true_states = n_true_states, 
  n_obs_states = n_obs_states, 
  first = first_weekly)# TODO sub this out for whatever time-scale you care about

#### RUN WEEKLY MODEL ####
##### MCMC settings ####
nb <- 5000 #burn-in # 5000
ni <- 5000 + nb #total iterations # 5000 + nb
nt <- 1  #thin
nc <- 3  #chains

##### Compile configuration and build ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      check = FALSE, calculate = FALSE, inits = inits)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt) # this step is slow
Rmcmc <- buildMCMC(conf)  

Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
beep(sound = 1)

##### Run MCMC ####
t.start <- Sys.time()
sink("Results/errorreport.txt") # for debugging, writing errors to txt file
out <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
               setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
t.end <- Sys.time()
sink() # closing txt file
(runTime <- t.end - t.start)
beep(sound = 1)

#### WEEKLY MODEL PLOTS ####
# Saves model and model, outputs, and diagnostics
MCMCdiag(out, dir = "./Results", save_object = TRUE, obj_name = "final_model.rds", file_name = "model_summary.txt")
# mod <- readRDS(here("Results", "final_model.rds")
# mod_summary <- read.lines(here("Results", "model_summary.txt"))

# Calculate numerical summaries
model_summary <- MCMCsummary(object = out, round = 3)

library(coda)
c1 <- as.mcmc(out$chain1)

# Caterpillar plot posterior distribution of phi
# Point = the posterior median,
# Thick line = 50% credible interval an
# Thin line = 95% credible interval
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
          params = c("phi", "mort.h"))

ex <- MCMCchains(out)

# Correlation plots
autocorr.diag(out) 
autocorr.plot(out)
coda::crosscorr.plot(out)

posterior <- as.array(out)

# Scatterplot of variables
color_scheme_set("purple")
phi_mort <- mcmc_scatter(out, pars = c("phi", "mort.h"))
phi_active <- mcmc_scatter(out, pars = c("phi", "active"))
mort_active <- mcmc_scatter(out, pars = c("mort.h", "active"))
p_q <- mcmc_scatter(out, pars = c("p", "q"))

out_pairsplot <- mcmc_pairs(out, diag_fun = "dens", off_diag_args = list(size = 0.5))

############################################################################# #
### WEEKLY COVARIATE MODEL  ###################################################
############################################################################# #

code <- nimbleCode({
  ##### Priors and constraints ####
  # SURVIVAL PROBABILITY
  mu.phi[1] ~ dunif(0,1) 
  mu.phi[2] ~ dunif(0,1) 
  mu.mort[1] ~ dunif(0,1)
  mu.mort[2] ~ dunif(0,1)
  
  # OBSERVATIONS
  active ~ dbeta(1, 1)
  mu.p[1] ~ dunif(0,1)
  mu.p[2] ~ dunif(0,1)
  mu.q[1] ~ dunif(0,1)
  mu.q[2] ~ dunif(0,1)
  mu.ro[1] ~ dunif(0,1)
  mu.ro[2] ~ dunif(0,1)
  mu.rh[1] <- mu.ro[1]
  mu.rh[2] <- mu.ro[2]
  
  ##### Transition matrix ####
  # 1 - alive and active collar
  # 2 - alive and inactive collar, observed
  # 3 - alive and inactive collar, undetected
  # 4 - newly caused natural mortality + active collar
  # 5 - newly caused natural mortality + inactive collar
  # 6 - newly caused human mortality + active collar
  # 7 - newly caused human mortality + inactive collar
  # 8 - long dead
  
  # Indexing
  # Including indexed values for phi, mort.h, p, q, r_O, r_H 
  # inside the loop for each individual 
  # priors for all mu. values are above 
  
  for (i in 1:n_ind) {
    phi[i] <- mu.phi[sex[i]]
    mort.h[i] <- mu.mort[sex[i]]
    p[i] <- mu.p[mfr[i]]  
    q[i] <- mu.q[mfr[i]]  
    r_O[i] <- mu.ro[mfr[i]]  
    r_H[i] <- r_O[i]
    
    # ALIVE AND ACTIVE (AA)
    trans.mat[1, 1, i] <- phi[i]*active             
    trans.mat[2, 1, i] <- 0                                                    
    trans.mat[3, 1, i] <- 0                                                 
    trans.mat[4, 1, i] <- 0 
    trans.mat[5, 1, i] <- 0
    trans.mat[6, 1, i] <- 0
    trans.mat[7, 1, i] <- 0
    trans.mat[8, 1, i] <- 0
    
    # ALIVE AND INACTIVE, OBSERVED (AI_O)
    trans.mat[1, 2, i] <- phi[i] * (1 - active)
    trans.mat[2, 2, i] <- 0 
    trans.mat[3, 2, i] <- 0 
    trans.mat[4, 2, i] <- 0 
    trans.mat[5, 2, i] <- 0 
    trans.mat[6, 2, i] <- 0 
    trans.mat[7, 2, i] <- 0 
    trans.mat[8, 2, i] <- 0 
    
    # ALIVE AND INACTIVE, UNDETECTED (AI_U)
    trans.mat[1, 3, i] <- 0
    trans.mat[2, 3, i] <- phi[i]
    trans.mat[3, 3, i] <- phi[i]
    trans.mat[4, 3, i] <- 0 
    trans.mat[5, 3, i] <- 0 
    trans.mat[6, 3, i] <- 0 
    trans.mat[7, 3, i] <- 0 
    trans.mat[8, 3, i] <- 0 
    
    # NEWLY DEAD-OTHER AND ACTIVE (DO_A)
    trans.mat[1, 4, i] <- (1-phi[i])*(1-mort.h[i])*active
    trans.mat[2, 4, i] <- 0 
    trans.mat[3, 4, i] <- 0 
    trans.mat[4, 4, i] <- 0 
    trans.mat[5, 4, i] <- 0 
    trans.mat[6, 4, i] <- 0 
    trans.mat[7, 4, i] <- 0 
    trans.mat[8, 4, i] <- 0 
    
    # NEWLY DEAD-OTHER AND INACTIVE (DO_I)
    trans.mat[1, 5, i] <- (1-phi[i])*(1-active)*(1-mort.h[i])
    trans.mat[2, 5, i] <- (1-phi[i])*(1-mort.h[i])
    trans.mat[3, 5, i] <- (1-phi[i])*(1-mort.h[i])
    trans.mat[4, 5, i] <- 0 
    trans.mat[5, 5, i] <- 0 
    trans.mat[6, 5, i] <- 0 
    trans.mat[7, 5, i] <- 0 
    trans.mat[8, 5, i] <- 0 
    
    # NEWLY DEAD-HUMAN AND ACTIVE (DH_A)
    trans.mat[1, 6, i] <- (1-phi[i])*mort.h[i]*active
    trans.mat[2, 6, i] <- 0
    trans.mat[3, 6, i] <- 0
    trans.mat[4, 6, i] <- 0
    trans.mat[5, 6, i] <- 0
    trans.mat[6, 6, i] <- 0
    trans.mat[7, 6, i] <- 0
    trans.mat[8, 6, i] <- 0
    
    # NEWLY DEAD-HUMAN AND INACTIVE (DH_I)
    trans.mat[1, 7, i] <- (1-phi[i])*mort.h[i]*(1-active)
    trans.mat[2, 7, i] <- (1-phi[i])*mort.h[i]
    trans.mat[3, 7, i] <- (1-phi[i])*mort.h[i]
    trans.mat[4, 7, i] <- 0
    trans.mat[5, 7, i] <- 0
    trans.mat[6, 7, i] <- 0
    trans.mat[7, 7, i] <- 0
    trans.mat[8, 7, i] <- 0
    
    # LONG DEAD (LD)
    trans.mat[1, 8, i] <- 0
    trans.mat[2, 8, i] <- 0
    trans.mat[3, 8, i] <- 0
    trans.mat[4, 8, i] <- 1
    trans.mat[5, 8, i] <- 1
    trans.mat[6, 8, i] <- 1
    trans.mat[7, 8, i] <- 1
    trans.mat[8, 8, i] <- 1
    
    ##### Observation matrix ####
    # 1 - alive and active collar
    # 2 - failed collar, observed
    # 3 - newly dead other mortality + active collar
    # 4 - newly dead human mortality + active collar
    # 5 - undetected
    
    # ALIVE AND ACTIVE
    obs.mat[1, 1, i] <- p[i]
    obs.mat[2, 1, i] <- 0
    obs.mat[3, 1, i] <- 0                                             
    obs.mat[4, 1, i] <- 0
    obs.mat[5, 1, i] <- 0
    obs.mat[6, 1, i] <- 0
    obs.mat[7, 1, i] <- 0
    obs.mat[8, 1, i] <- 0
    
    # Fail
    obs.mat[1, 2, i] <- 0
    obs.mat[2, 2, i] <- q[i]
    obs.mat[3, 2, i] <- 0
    obs.mat[4, 2, i] <- 0
    obs.mat[5, 2, i] <- 0
    obs.mat[6, 2, i] <- 0
    obs.mat[7, 2, i] <- 0
    obs.mat[8, 2, i] <- 0
    
    # NEWLY DEAD-OTHER 
    obs.mat[1, 3, i] <- 0
    obs.mat[2, 3, i] <- 0 
    obs.mat[3, 3, i] <- 0 
    obs.mat[4, 3, i] <- r_O[i]
    obs.mat[5, 3, i] <- 0 
    obs.mat[6, 3, i] <- 0 
    obs.mat[7, 3, i] <- 0 
    obs.mat[8, 3, i] <- 0 
    
    # NEWLY DEAD-HUMAN 
    obs.mat[1, 4, i] <- 0
    obs.mat[2, 4, i] <- 0  
    obs.mat[3, 4, i] <- 0  
    obs.mat[4, 4, i] <- 0
    obs.mat[5, 4, i] <- 0
    obs.mat[6, 4, i] <- r_H[i]
    obs.mat[7, 4, i] <- 0
    obs.mat[8, 4, i] <- 0
    
    # NOT DETECTED
    obs.mat[1, 5, i] <- 1-p[i]
    obs.mat[2, 5, i] <- 1-q[i]
    obs.mat[3, 5, i] <- 1 
    obs.mat[4, 5, i] <- 1-r_O[i]
    obs.mat[5, 5, i] <- 1
    obs.mat[6, 5, i] <- 1-r_H[i]
    obs.mat[7, 5, i] <- 1
    obs.mat[8, 5, i] <- 1
    
    # END priors and constraints
    
    ##### Likelihood #####
    
    for (t in (first[i] + 1):n_occasions) { # conditioning on first capture 
      z[i, t] ~ dcat(trans.mat[z[i, t-1], 1:n_true_states, i]) # z is latent (true state); categorical distribution
      y[i, t] ~ dcat(obs.mat[z[i, t], 1:n_obs_states, i]) # y is observed data; categorical distribution
    } # t last
  } # i
  
  # END likelihood
})

#### INIITAL VALUES ####
# Data is coded as OBSERVATION states. Providing model with initial values for deterministic TRANSITION states

# TRANSITION STATES
# 1 - alive and active collar
# 2 - alive and inactive collar, observed
# 3 - alive and inactive collar, undetected
# 4 - newly caused natural mortality + active collar
# 5 - newly caused natural mortality + inactive collar
# 6 - newly caused human mortality + active collar
# 7 - newly caused human mortality + inactive collar
# 8 - long dead

##### Defining initial states #####
z_inits <- y
z_dat <- y
z_inits[y == 999] <- NA
z_dat[y == 999] <- NA

# active collar -> alive (1 -> 1)
z_inits[y == 1] <- NA
z_dat[y == 1] <- 1
# collar failed -> alive and observed inactive collar (2 -> 2)
z_inits[y == 2] <- NA
z_dat[y == 2] <- 2
# newly caused natural mortality -> newly caused natural mortality (3 -> 4)
z_inits[y == 3] <- NA
z_dat[y == 3] <- 4
# newly caused human mortality -> newly caused human mortality (4 -> 6)
z_inits[y == 4] <- NA
z_dat[y == 4] <- 6
# undetected -> alive and unobserved inactive collar (5 -> 1)
z_inits[y == 5] <- 1
z_dat[y == 5] <- NA

# If long inactive after mortality signal observed -> long dead (5 -> 8)
for (i in 1:n_inds) {
  tmp <- z_dat[i, ]
  if (any(tmp == 4, na.rm = TRUE) | any(tmp == 6, na.rm = TRUE)) { # if long inactive after mortality signal observed
    # Index of 3 of 4 (dead)
    indexD <- which(tmp == 4 | tmp == 6)
    # Indices of 5's after 4|6 (dead) -> replace with 8 (LD)
    indexLD <- (indexD + 1):length(tmp) 
    tmp[indexLD] <- 8
  }
  z_dat[i,] <- tmp
}

z_inits[z_dat == 8] <- NA

# If collar failure (confirmed inactive) -> all values after are inactive but undetected (2 -> 3)
for (i in 1:n_inds) {
  tmp <- z_inits[i,]
  if (any(y[i,] == 2, na.rm = TRUE)) {  
    # Index of 2 inactive confirmed
    indexD <- which(y[i,] == 2)
    # Indices after 2 inactive -> replace with 3
    indexLD <- (indexD + 1):length(tmp) 
    tmp[indexLD] <- 3
  } 
  z_inits[i,] <- tmp
}

z_dat[z_inits == 3] <- NA

##### Inits #####
inits <- list(
  z = z_inits,
  mu.phi = rbeta(2, 1, 1), 
  mu.mort = rbeta(2, 1, 1), 
  active = rbeta(1, 1, 1),
  mu.p = rbeta(2, 1, 1), # detection prob alive + active
  mu.q = rbeta(2, 1, 1)) # detection prob collar failure

#### DATA ####
y[y == 999] <- NA
dat <- list(
  sex = X$sex %>% as_factor() %>% as.numeric(),
  mfr = X$mfr %>% as_factor() %>% as.numeric(),
  y = y,
  z = z_dat) 

#### PARAMETERS TO MONITOR ####
# put anything that you want to track in here
# the more parameters, the slower the model
# best practice to only include things that are directly sampled (i.e., you have a prior for them)
# anything derived can be done post-hoc, unless you want the model to spit out annual survival probability
# when debugging model, could add trans.mat and obs.mat, or even z, etc. 
params <- c(
  "mu.phi", 
  "mu.mort",
  "active",
  "mu.p",  
  "mu.q", 
  "mu.ro")

#### CONSTANTS ####
const <- list(
  n_ind = n_inds, 
  n_occasions = dim(dat$y)[2], 
  n_true_states = n_true_states, 
  n_obs_states = n_obs_states, 
  first = first_weekly, # TODO sub this out for whatever time-scale you care about 
  n_sex = 2,
  n_mfr = 2
)

#### RUN WEEKLY MODEL ####
##### MCMC settings ####
nb <- 5000 #burn-in # 5000
ni <- 5000 + nb #total iterations # 5000 + nb
nt <- 1  #thin
nc <- 3  #chains

##### Compile configure and build ####
Rmodel <- nimbleModel(code = code, constants = const, data = dat, 
                      check = FALSE, calculate = FALSE, inits = inits)
conf <- configureMCMC(Rmodel, monitors = params, thin = nt) # this step is slow
Rmcmc <- buildMCMC(conf)  

Cmodel <- compileNimble(Rmodel, showCompilerOutput = FALSE)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
beep(sound = 1)

##### Run  MCMC ####
t.start <- Sys.time()
sink("Results/errorreport.txt") # for debugging, writing errors to txt file
out_cov <- runMCMC(Cmcmc, niter = ni , nburnin = nb , nchains = nc, inits = inits,
                   setSeed = FALSE, progressBar = TRUE, samplesAsCodaMCMC = TRUE)  
t.end <- Sys.time()
sink() # closing txt file
(runTime <- t.end - t.start)
beep(sound = 1)

#### WEEKLY COVARIATE MODEL PLOTS ####
# Saves model and model, outputs, and diagnostics
MCMCdiag(out_cov, dir = "./Results", save_object = TRUE, obj_name = "final_model_covariates.rds", file_name = "model_summary_covariates.txt")
# mod_cov <- readRDS(here("Results", "final_model_covariates.rds"))
# mod_summary_cov <- read_table(here("Results", "model_summary_covariates.txt"), 
#                               skip = 11)

# Calculate numerical summaries
model_summary_cov <- MCMCsummary(object = out_cov, round = 3)

# Caterpillar plot posterior distribution of phi
# Point = the posterior median,
# Thick line = 50% credible interval an
# Thin line = 95% credible interval
MCMCplot(object = out_cov,
         horiz = FALSE,
         rank = TRUE,
         ref_ovl = FALSE)

# Check convergence of chains
MCMCtrace(object = out_cov,
          pdf = FALSE, # no export to PDF
          ind = TRUE, # separate density lines per chain
          Rhat = TRUE, # add Rhat diagnostics
          n.eff = TRUE, # add eff sample size
          params = c("mu.phi", "mu.mort"))

ex <- MCMCchains(out_cov)

# Some post-hoc tests and plots
autocorr.diag(out_cov) 
autocorr.plot(out_cov)

posterior <- as.array(out_cov)

# scatterplot of variables
color_scheme_set("pink")
phi_mort_cov <- mcmc_scatter(out_cov, pars = c("mu.phi", "mu.mort"))
phi_active_cov <- mcmc_scatter(out_cov, pars = c("mu.phi", "active"))
mort_active_cov <- mcmc_scatter(out_cov, pars = c("mu.mort", "active"))
p_q_cov <- mcmc_scatter(out_cov, pars = c("mu.p", "mu.q"))
out_pairsplot_cov <- mcmc_pairs(out_cov, diag_fun = "dens", off_diag_args = list(size = 0.5))
