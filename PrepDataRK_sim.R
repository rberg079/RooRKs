#' Prepare data for road mortality model
#'
#' @param females logical. If TRUE, data is prepped for female kangaroos only. females = TRUE by default.
#'
#' @returns a list containing all data & constants needed for known-fate road mortality models.
#' @export
#'
#' @examples

prepDataRK <- function(females = T){
  
  # # for testing purposes
  # females = FALSE
  
  library(readxl)
  library(tidyverse)
  library(lubridate)
  
  
  ## Load & clean up data ------------------------------------------------------
  
  yafs <- read_excel("data/RSmain_Jan26.xlsx")
  surv <- read_excel("data/PromSurvivalNov25_RB.xlsx", sheet = "YEARLY SURV")
  obs <- read_excel("data/PromObs_2008-2024.xlsx")
  env <- read_csv("data/Env_Mar25.csv")
  
  # modify column names to read as survived to s20XX
  surv <- surv %>%
    select(1, 7, 9:45) %>%
    rename("2009" = "08to09", "2010" = "09to10", "2011" = "10to11",
           "2012" = "11to12", "2013" = "12to13", "2014" = "13to14",
           "2015" = "14to15", "2016" = "15to16", "2017" = "16to17",
           "2018" = "17to18", "2019" = "18to19", "2020" = "19to20",
           "2021" = "20to21", "2022" = "21to22", "2023" = "22to23",
           "2024" = "23to24", "2025" = "24to25", Dead = "Found dead") %>% 
    mutate(Sex = Sex-1)
  
  # select sex
  if(females){
    surv <- surv %>%
      filter(Sex == 1) # females
  }else{
    surv <- surv %>%
      filter(Sex == 0) # males
  }
  
  
  ## YAF survival data ---------------------------------------------------------
  
  # filter out the too young/old & the firstborn "twins"
  # & select YAFs that made it to their first October
  yafs <- yafs %>% 
    mutate(Age = as.numeric(Age)) %>% 
    filter(Exclude == 0,
           between(Age, 3, 20) | is.na(Age),
           PYid != 308 & PYid != 340 & PYid != 672 & PYid != 885 & PYid != 900 &
           PYid != 891 & PYid != 912 & PYid != 1023 & PYid != 1106,
           SurvOct1 == 1) %>% 
    select(Year, PYsex, PYid, SurvOct2) %>% 
    rename(Sex = PYsex, ID = PYid) %>% 
    mutate(Sex = Sex-1)
  
  # select sex
  if(females){
    yafs <- yafs %>%
      filter(Sex == 1) %>%
      select(-Sex) # females
  }else{
    yafs <- yafs %>%
      filter(Sex == 0) %>%
      select(-Sex) # males
  }
  
  # align survival to second October
  # to the corresponding year
  yafs <- yafs %>% 
    rename(yafs = SurvOct2) %>% 
    mutate(Year = Year + 1,
           yafs = as.numeric(yafs),
           yafs = ifelse(yafs == 2, NA, yafs)) %>% 
    filter(!is.na(yafs))
  
  
  ## Encounter history ---------------------------------------------------------
  
  # round up survival data
  eh <- surv %>%
    select(1,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39) %>%
    mutate("2008" = NA) %>%
    select(1,19,2:18) %>%
    pivot_longer(-ID, names_to = "Year", values_to = "eh") %>%
    mutate(Year = as.numeric(Year))
  
  # # check for yafs that were recovered dead
  # # first pull all yafs not already in eh
  # check <- yafs %>%
  #   distinct(ID) %>%
  #   anti_join(eh %>% distinct(ID),
  #             by = "ID")
  
  # then manually check from fine-scale survival data
  # whether each yaf in check was recovered or disappeared
  deadYAFs <- c(223, 278, 285, 299, 330, 362, 382, 449, 452, 465, 467, 473, 479,
                489, 502, 527, 530, 532, 539, 569, 570, 573, 575, 582, 583, 589,
                590, 592, 593, 599, 605, 606, 615, 622, 628, 632, 634, 642, 681,
                693, 703, 717, 723, 724, 748, 758, 759, 774, 819, 850,
                1041, 1051, 1053, 1159, 1261, 1270, 1277, 1293, 1321)
  
  # join yafs data
  eh <- eh %>% 
    mutate(Year = as.numeric(Year)) %>% 
    full_join(yafs, by = c("ID", "Year")) %>% 
    mutate(surv = coalesce(yafs, eh)) %>% 
    select(ID, Year, surv)
  
  # complete ID*year pairs & add 1s
  # at the start of each survival interval
  eh <- eh %>% 
    complete(ID, Year = 2008:2025) %>%
    arrange(ID, Year) %>%
    group_by(ID) %>%
    mutate(surv = ifelse(is.na(surv) & !is.na(lead(surv)), 1, surv)) %>%
    ungroup() %>%
    pivot_wider(names_from = Year, values_from = surv)
  
  id <- eh %>% select(1)
  eh <- eh %>% mutate(ID = 999)
  
  # change values before first obs to 999
  for(i in 1:nrow(eh)){
    eh[i, 1:(min(which(eh[i,] == 1)) -1)] <- 999
  }
  
  eh <- eh %>% select(2:19) # fix later
  
  # OBSERVATION STATES
  # 1 - seen
  # 2 - recovered roadkill
  # 3 - undetected
  
  # gather ID & which were found dead
  id <- left_join(id, surv) %>% 
    select(ID, Dead) %>%
    mutate(Dead = ifelse(!is.na(Dead), 1, 0),
           Dead = ifelse(ID %in% deadYAFs, 1, Dead))
  
  eh <- cbind(id, eh)
  
  # create a Fate column
  # Fate = 2 when found dead of human cause
  # Fate = 3 when found dead of natural cause
  eh <- eh %>%
    rowwise() %>%
    mutate(Fate = ifelse(any(c_across(3:20) == 2), 2, NA),
           Fate = ifelse(is.na(Fate) & Dead == 1, 3, Fate),
           Fate = ifelse(ID == 1153, 2, Fate))
  
  # SPECIAL CASES: change Dead to 1 for handful of RKs
  # that were not recovered according to survival file
  eh <- eh %>%
    mutate(Dead = ifelse(!is.na(Fate), 1, Dead),
           `2022` = ifelse(ID == 1293, 0, `2022`))
  # ID 1293 should have been in the survival file
  # & therefore gotten a 0 in first year unobserved
  
  id <- id %>%
    select(ID) %>%
    left_join(., eh) %>%
    select(ID, Dead)
  
  # now that cause of death is stored in the Fate column
  # replace observations of 0, 2, & 3 with observation state 4
  # currently 2s are at first October when inds were not alive to be seen
  eh <- eh %>%
    mutate_at(3:20, ~replace(., . == 0, 4)) %>% # 0: new natural death
    mutate_at(3:20, ~replace(., . == 2, 4)) %>% # 2: new roadkill
    mutate_at(3:20, ~replace(., . == 3, 1)) %>% # 3: seen off site
    mutate_at(3:20, ~replace(., . == 4, 4))     # 4: undetected emigrant
         
  # create column Gone with year where each ID is first unobserved
  # Gone is NA for individuals still alive & observed in 2025
  suppressWarnings(
    eh <- eh %>%
      rowwise() %>%
      mutate(Gone = max(which(c_across(3:20) == 4)),
             Gone = ifelse(Gone == "-Inf" | ID == 834 | ID == 1125, NA, Gone)) %>%
      ungroup()) # suppose to give -Inf warnings
  
  # replace last observation with Fate
  for(i in 1:nrow(eh)) {
    fate <- eh$Fate
    gone <- eh$Gone
    
    if(!is.na(gone[i]) && !is.na(fate[i])) {
      eh[i, gone[i]+1] <- fate[i]
    }
  }
  
  # select relevant columns
  # replace remaining NAs with obs 4
  id <- cbind(id, fate, gone)
  remove(i, fate, gone)
  
  eh <- eh %>%
    select(-c(Dead, Fate, Gone)) %>%
    mutate(across(2:19, ~replace(., is.na(.), 4)))
  
  # eh <- eh %>% select(1:19)
  # write_csv(eh, "eh.csv")
  # remove(yafs, deadYAFs)
  
  
  ## Individual data -----------------------------------------------------------
  
  # # observer days data
  # obs <- obs %>%
  #   select(Date, Month, Day, Year, Time, ID, X, Y, Observer) %>% 
  #   mutate(Time = format(as.POSIXct(Time), format = "%H:%M")) %>% 
  #   distinct(Year, Date, Observer) %>% 
  #   group_by(Year) %>% 
  #   summarise(n = n()) %>% 
  #   ungroup()
  # 
  # # TEMP; until I figure out Emily's obs*days
  # obs <- rbind(obs, c(2021, 10)) %>% 
  #   arrange(Year)
  
  # age data
  age <- surv %>%
    select(1,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38) %>%
    mutate_at(vars(-ID), ~ifelse(.== "A", NA, .)) %>%
    mutate_at(vars(-ID), ~as.numeric(.)) %>%
    mutate("2025" = as.numeric(NA)) %>%
    rename("2008" = "Age08", "2009" = "Age09", "2010" = "Age10",
           "2011" = "Age11", "2012" = "Age12", "2013" = "Age13",
           "2014" = "Age14", "2015" = "Age15", "2016" = "Age16",
           "2017" = "Age17", "2018" = "Age18", "2019" = "Age19",
           "2020" = "Age20", "2021" = "Age21", "2022" = "Age22",
           "2023" = "Age23", "2024" = "Age24")
  
  yafs <- yafs %>%
    pivot_wider(names_from = Year, values_from = yafs) %>%
    mutate_at(vars(-ID), ~ifelse(!is.na(.), 1, NA)) %>%
    anti_join(age %>% select(ID), by = "ID") %>%
    mutate("2008" = NA) %>%
    select(ID, "2008", everything())

  age <- bind_rows(age, yafs) %>% arrange(ID)
  tmp <- age %>% select(2:19)
  
  for (j in 1:nrow(tmp)) {
    # forward fill
    for (i in 2:ncol(tmp)) {
      if (is.na(tmp[j, i])) {
        tmp[j, i] <- tmp[j, i - 1] + 1
      }
    }
    
    # backward fill
    for (i in (ncol(tmp) - 1):1) {
      if (is.na(tmp[j, i])) {
        tmp[j, i] <- tmp[j, i + 1] - 1
      }
    }
  }
  
  tmp[tmp < 0] <- NA
  age <- tmp
  remove(tmp, i, j, surv, yafs, deadYAFs)
  # write_csv(age, "ageF.csv")
  
  age <- as.matrix(age)
  
  # age classes
  # when including YAFs
  age <- age + 1
  ageC <- c(0, 1, 2, rep(3,4), rep(4,3), rep(5,30)) + 1 
  
  # when excluding YAFs
  # ageC <- c(1, 2, rep(3,4), rep(4,3), rep(5,31))
  
  
  ## X coordinate data ---------------------------------------------------------
  
  obs <- obs %>% 
    select(Date, Year, Month, Day, Time, ID, X) %>% 
    mutate(ttime = format(as.POSIXct(Time), format = "%H:%M"),
           X = as.numeric(X)) %>% 
    select(-Time) %>% 
    rename(Time = ttime) %>% 
    filter(X < 40000, X > 32000,  # within the realm of reality
           !is.na(ID), !is.na(X), # only known roos w/ coordinates
           Month >= 7)            # during the main field season
  
  # limit to IDs seen at least 10 times during the year
  # calculate from those 10+ observations a median X
  obs <- obs %>%
    group_by(Year, ID) %>% 
    mutate(DaysObs = n_distinct(Date)) %>% 
    ungroup() %>% 
    group_by(ID, Year) %>% 
    mutate(xmed = median(X, na.rm = T)) %>% 
    ungroup() %>% 
    arrange(ID, Year) %>% 
    filter(DaysObs >= 10)
  
  # hist(obs$X)
  # hist(obs$xmed)
  
  # set up matrix
  obs <- obs %>% 
    arrange(Year) %>% 
    distinct(Year, ID, xmed) %>% 
    complete(ID, Year = 2008:2025) %>%
    pivot_wider(names_from = Year, values_from = xmed) %>% 
    filter(ID %in% id$ID)
  
  # impute missing values
  # for IDs with some xmeds
  tmp <- obs %>% select(ID)
  obs <- as.matrix(obs[, 2:19])
  obs[is.na(obs)] <- rowMeans(obs, na.rm = T)[row(obs)][is.na(obs)]
  obs <- as.data.frame(obs)
  tmp <- cbind(tmp, obs)
  
  # join with ID
  xmed <- id %>% 
    select(ID) %>% 
    left_join(tmp)
  remove(tmp, obs)
  
  xmed <- scale(as.matrix(xmed[, 2:19]))
  
  
  ## Environmental data --------------------------------------------------------
  
  env <- env %>%
    mutate(Year = ifelse(Month < 10, Year-1, Year)) %>%
    group_by(Year) %>%
    mutate(Veg = sum(Veg, na.rm = T),
           Dens = mean(Dens, na.rm = T),
           Win = sum(Warn.18, na.rm = T),
           Veg = ifelse(between(Year, 2009, 2023), Veg, NA),
           Dens = ifelse(between(Year, 2008, 2024), Dens, NA),
           Win = ifelse(between(Year, 2008, 2023), Win, NA)) %>%
    distinct(Year, Veg, Dens, Win) %>%
    mutate(VegRoo = Veg/Dens) %>% 
    ungroup()

  env  <- env[3:19,] # [2008:2024,]
  # obs  <- round(as.numeric(scale(obs$n)), 3)
  veg  <- round(as.numeric(scale(env$Veg)), 3)
  dens <- round(as.numeric(scale(env$Dens)), 3)
  vegr <- round(as.numeric(scale(env$VegRoo)), 3)
  win  <- round(as.numeric(scale(env$Win)), 3)

  noVeg  <- c(as.numeric(which(is.na(veg))));  nNoVeg  <- length(noVeg)
  noDens <- c(as.numeric(which(is.na(dens)))); nNoDens <- length(noDens)
  noVegR <- c(as.numeric(which(is.na(vegr)))); nNoVegR <- length(noVegR)
  noWin  <- c(as.numeric(which(is.na(win))));  nNoWin  <- length(noWin)
  
  # write_csv(env, "env.csv")
  
  
  ## Assemble data list --------------------------------------------------------
  
  y <- as.matrix(unname(eh[,-1]))
  
  # y[age == 0] <- 999  # when excluding YAFs
  
  # extract first & last
  # create vector with occasion of first capture
  get.first <- function(x) min(which(x != 999))
  get.last  <- function(x) max(which(x < 4))
  
  first <- apply(y, 1, get.first)
  last  <- apply(y, 1, get.last)
  
  # missing age
  unkAge <- which(rowSums(is.na(age)) == ncol(age))
  
  # died at first occasion
  badFirst <- which(y[cbind(seq_len(nrow(y)), first)] != 1)
  
  # thus, usable individuals are
  # y <- y[-unkAge, ]
  # id <- id[-unkAge, ]
  # age <- age[-unkAge, ]
  # n.inds <- nrow(id)
  
  y <- y[-c(unkAge, badFirst), ]
  id <- id[-c(unkAge, badFirst), ]
  age <- age[-c(unkAge, badFirst), ]
  xmed <- xmed[-c(unkAge, badFirst),]
  eh <- eh[-c(unkAge, badFirst), ]
  
  noX  <- c(as.numeric(which(is.na(xmed))))
  nNoX <- length(noX)
  
  # # check sample size
  # RKtable <- eh %>%
  #   pivot_longer(-ID, names_to = "year", values_to = "eh") %>%
  #   filter(eh == 2) %>%
  # 
  #   left_join(age %>%
  #               as.data.frame() %>%
  #               mutate(ID = eh$ID) %>%
  #               pivot_longer(-ID, names_to = "year", values_to = "age"),
  #             by = c("ID", "year")) %>%
  # 
  #   left_join(id, by = "ID") %>%
  #   count(year, age, sort = TRUE) %>%
  #   mutate(age = age-1) %>%
  #   arrange(year, age)
  # 
  # RKtable
  
  first <- apply(y, 1, get.first)
  last  <- apply(y, 1, get.last)
  
  n.inds <- nrow(id)
  n.ageC <- max(ageC)
  n.occasions <- ncol(y)
  n.obs.states <- 3
  n.true.states <- 4
  
  # assemble list
  dataRK <- list(
    eh = eh,
    id = id,
    y = y,
    age = age,
    ageC = ageC,
    first = first,
    last = last,
    xmed = xmed,
    noX = noX,
    nNoX = nNoX,
    
    veg = veg,
    dens = dens,
    vegr = vegr,
    win = win,
    noVeg = noVeg,
    noDens = noDens,
    noVegR = noVegR,
    noWin = noWin,
    nNoVeg = nNoVeg,
    nNoDens = nNoDens,
    nNoVegR = nNoVegR,
    nNoWin = nNoWin,
    
    n.inds = n.inds,
    n.ageC = n.ageC,
    n.occasions = n.occasions,
    n.obs.states = n.obs.states,
    n.true.states = n.true.states
    # badFirst = badFirst,
    # unkAge = unkAge
  )
  
}

