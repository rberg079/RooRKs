#' Simulate data for RK model
#'
#' @param n.age integer. Number of ages. n.age = 20 by default.
#' @param n.ageC integer. Number of age classes. n.ageC = 5 by default.
#' @param n.occasions integer. Number of years. n.occasions = 20 by default.
#' @param first vector of length n.inds. Year when each roo was first caught. first = rep(1:n.occasions, each = n.age) by default.
#' @param mu.age vector of length n.ageC. Mean survival probability of each age class. mu.age = c(0.7, 0.85, 0.9, 0.9, 0.8) by default.
#' @param B.veg vector of length n.ageC. Age-specific slope of vegetation effect. B.veg = c(0.6, 0.4, 0.2, 0.2, 0.4) by default.
#' @param mean.veg integer. Mean value of centered & scaled vegetation. mean.veg = 0 by default.
#' @param sd.veg integer. Standard deviation of vegetation effect. sd.veg = 1 by default.
#' @param sigma.phi integer. Standard deviation of random year effect on survival. sigma.phi = 0.2 by default.
#' @param mean.R integer. Mean probability of a mortality being road-related. mean.R = 0.05 by default.
#' @param mean.M integer. Mean probability of migration, in either direction. mean.M = 0.05 by default.
#' @param mean.Pi vector of length n.occasions. Probability of observation in the study area. mean.Pi ~ N(qlogis(0.6), 0.2) by default.
#' @param mean.Po vector of length n.occasions. Probability of observation outside the study area. mean.Po ~ N(qlogis(0.4), 0.2) by default.
#' @param mean.rR vector of length n.occasions. Probability of recovery of roadkilled roos. mean.rR ~ N(qlogis(0.6), 0.2) by default.
#' @param mean.rO vector of length n.occasions. Probability of recovery of other dead roos. mean.rO ~ N(qlogis(0.4), 0.2) by default.
#' @param seed integer. Seed. seed = 42 by default.
#'
#' @returns a list containing all simulated data needed to run the RK model.
#' @export
#'
#' @examples

simulateDataRK <- function(n.age = 20,
                           n.ageC = 5,
                           n.occasions = 20,
                           first = NULL,
                           mu.age = c(0.7, 0.85, 0.9, 0.9, 0.8),
                           B.veg = c(0.6, 0.4, 0.2, 0.2, 0.4),
                           mean.veg = 0,
                           sd.veg = 1,
                           sigma.phi = 0.2,
                           mean.R = 0.05,
                           mean.M = 0.05,
                           mean.Pi = NULL,
                           mean.Po = NULL,
                           mean.rR = NULL,
                           mean.rO = NULL,
                           seed = 42){
  
  # # for testing purposes
  # library(dplyr)
  # library(purrr)
  # library(tidyr)
  # 
  # n.age = 20
  # n.ageC = 5
  # n.occasions = 20
  # n.inds = n.age * (n.occasions-1)
  # first = rep(1:(n.occasions-1), each = n.age)
  # mu.age = c(0.7, 0.85, 0.9, 0.9, 0.8)
  # B.veg = c(0.6, 0.4, 0.2, 0.2, 0.4)
  # mean.veg = 0
  # sd.veg = 1
  # sigma.phi = 0.2
  # mean.R = 0.05
  # mean.M = 0.05
  # mean.Pi = NULL
  # mean.Po = NULL
  # mean.rR = NULL
  # mean.rO = NULL
  # seed = 42
  
  set.seed(seed)
  n.inds = n.age * (n.occasions-1)
  
  if(is.null(mean.Pi)) mean.Pi <- plogis(rnorm(n.occasions, qlogis(0.6), 0.2))
  if(is.null(mean.Po)) mean.Po <- plogis(rnorm(n.occasions, qlogis(0.4), 0.2))
  if(is.null(mean.rR)) mean.rR <- plogis(rnorm(n.occasions, qlogis(0.6), 0.2))
  if(is.null(mean.rO)) mean.rO <- plogis(rnorm(n.occasions, qlogis(0.4), 0.2))
  
  if(is.null(first)) first <- rep(1:(n.occasions-1), each = n.age)
  
  # simulate standard normal veg
  veg <- rnorm(n.occasions - 1, mean = mean.veg, sd = sd.veg)
  
  
  ## Survival function ---------------------------------------------------------
  
  # convert mu.age to logit
  logit <- function(p) log(p / (1-p))
  invlogit <- function(x) 1 / (1 + exp(-x))
  mu.logit.age <- logit(mu.age)
  
  # eps.phi random year effects
  eps.phi <- matrix(rnorm(n.ageC * (n.occasions - 1), 0, sigma.phi),
                    nrow = n.ageC, ncol = n.occasions - 1)
  
  # mean.phi per age class & time
  mean.phi <- matrix(NA, nrow = n.ageC, ncol = n.occasions - 1)
  for(a in 1:n.ageC){
    for(t in 1:(n.occasions-1)){
      mean.phi[a,t] <- invlogit(mu.logit.age[a] + B.veg[a] * veg[t] + eps.phi[a,t])
    }
  }
  
  
  ## Age matrix ----------------------------------------------------------------
  
  # prep age matrix
  ageC <- c(0, rep(1,2), rep(2,4), rep(3,3), rep(4,10))+1 # age classes
  age <- matrix(NA_integer_, nrow = n.inds, ncol = n.occasions)
  
  # set initial age
  firstAge <- sort(rpois(n.age, 2), decreasing = T)
  firstAge <- firstAge/sum(firstAge) # as probs summing to 1
  firstAge <- sample(1:n.age, n.inds, replace = TRUE, prob = firstAge)
  # hist(firstAge)
  
  # fill age matrix
  for(i in 1:n.inds){
    for(t in first[i]:n.occasions){          # from first capture
      age.it <- firstAge[i] + (t - first[i]) # age each year
      age.it <- pmin(age.it, n.age)          # cap at max age
      age[i,t] <- age.it                     # fill into age.mat
    }
  }
  
  
  ## Transitions & observations ------------------------------------------------
  
  # transition from t -> t + 1
  # function to draw z given previous z & phi at t
  sample.z <- function(z.prev, phi, R = mean.R, M = mean.M){
    p <- numeric(5)
    if(z.prev == 1){ 
      # if z.prev is 1,
      # take transition probabilities for state 1
      p[1] <- phi * (1 - M)
      p[2] <- phi * M
      p[3] <- (1-phi) * R
      p[4] <- (1-phi) * (1-R)
      p[5] <- 0
    } else if(z.prev == 2){ 
      # if z.prev is 2,
      # take transition probabilities for state 2
      p[1] <- phi * M
      p[2] <- phi * (1 - M)
      p[3] <- (1-phi) * R
      p[4] <- (1-phi) * (1-R)
      p[5] <- 0
    } else if(z.prev %in% c(3,4,5)){ 
      # if z.prev is 3, 4, or 5,
      # only possible transition is to state 5
      p <- c(0,0,0,0,1)
    }
    # # numerical safety
    # if(sum(p) <= 0){
    #   p <- p + 1e-12
    # }
    p <- p / sum(p)
    sample(1:5, size = 1, prob = p)
  }
  
  # initialize z & y matrices
  z <- matrix(NA_integer_, nrow = n.inds, ncol = n.occasions)
  y <- matrix(NA_integer_, nrow = n.inds, ncol = n.occasions)
  
  # simulate initial z & y
  # individuals are alive & observed at first[i]
  for(i in 1:n.inds){
    z[i, first[i]] <- 1
    y[i, first[i]] <- 1
  }
  
  # simulate z & y forward with sample.z function
  for(i in 1:n.inds){
    if(first[i] + 1 > n.occasions) next
    for(t in (first[i]+1):n.occasions){
      
      # transition from t-1 -> t
      age.prev <- ageC[age[i, t-1]]
      phi.prev <- mean.phi[age.prev, t-1]
      
      # draw z at t
      z.next <- sample.z(z.prev = z[i, t-1],
                         phi = phi.prev,
                         R = mean.R,
                         M = mean.M)
      z[i,t] <- z.next
      
      # draw y at t
      Pi <- mean.Pi[t]
      Po <- mean.Po[t]
      rR <- mean.rR[t]
      rO <- mean.rO[t]
      
      if(z.next == 1){
        probs.obs <- c(Pi, 0,  0,  0,  1-Pi)
      }else if(z.next == 2){
        probs.obs <- c(0,  Po, 0,  0,  1-Po)
      }else if(z.next == 3){
        probs.obs <- c(0,  0,  rR, 0,  1-rR)
      }else if(z.next == 4){
        probs.obs <- c(0,  0,  0,  rO, 1-rO)
      }else if(z.next == 5){
        probs.obs <- c(0,  0,  0,  0,  1)
      }
      y[i,t] <- sample(1:5, 1, prob = probs.obs)
    }
  }
  
  
  ## Assemble list -------------------------------------------------------------
  
  # code missing values as 999
  y[is.na(y)] <- 999
  
  n.true.states <- 5
  n.obs.states <- 5
  
  out <- list(
    n.age = n.age,
    n.ageC = n.ageC,
    n.inds = n.inds,
    n.occasions = n.occasions,
    n.true.states = n.true.states,
    n.obs.states = n.obs.states,
    
    y = y,
    z = z,
    ageC = ageC,
    age = age,
    veg = veg,
    first = first,
    mean.phi = mean.phi,
    eps.phi = eps.phi,
    mean.Pi = mean.Pi,
    mean.Po = mean.Po,
    mean.rR = mean.rR,
    mean.rO = mean.rO,
    mean.R = mean.R,
    mean.M = mean.M,
    params = list(B.veg = B.veg,
                  mu.age = mu.age,
                  sigma.phi = sigma.phi)
  )
  return(out)
}

