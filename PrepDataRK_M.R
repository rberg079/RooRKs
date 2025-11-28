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
  # females = TRUE
  
  library(readxl)
  library(tidyverse)
  library(lubridate)
  
  
  ## Load & clean up data ------------------------------------------------------
  
  surv <- read_excel("data/PromSurvivalOct24.xlsx", sheet = "YEARLY SURV")
  if(females){age <- read_csv("data/ageF.csv")}else{age <- read_csv("data/ageM.csv")}
  obs <- read_excel("data/PromObs_2008-2024.xlsx")
  env <- read_csv("data/Env_Mar25.csv")
  
  # modify column names to read as survived to s20XX
  surv <- surv %>%
    select(1, 7, 9:43) %>%
    rename(s2009 = "08to09", s2010 = "09to10", s2011 = "10to11",
           s2012 = "11to12", s2013 = "12to13", s2014 = "13to14",
           s2015 = "14to15", s2016 = "15to16", s2017 = "16to17",
           s2018 = "17to18", s2019 = "18to19", s2020 = "19to20",
           s2021 = "20to21", s2022 = "21to22", s2023 = "22to23",
           s2024 = "23to24", Dead = "Found dead") %>% 
    mutate(Sex = Sex-1)
  
  if(females){
    surv <- surv %>% filter(Sex == 1) # females
  }else{
    surv <- surv %>% filter(Sex == 0) # males
  }
  
  
  ## Encounter history ---------------------------------------------------------
  
  eh <- surv %>%
    select(1,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37) %>%
    mutate(s2008 = NA) %>%
    select(1,18,2:17)
  
  id <- eh %>% select(1)
  
  eh <- eh %>% 
    mutate(s2008 = ifelse(!is.na(s2009), 1, NA),
           s2009 = ifelse(!is.na(s2010) & is.na(s2009), 1, s2009),
           s2010 = ifelse(!is.na(s2011) & is.na(s2010), 1, s2010),
           s2011 = ifelse(!is.na(s2012) & is.na(s2011), 1, s2011),
           s2012 = ifelse(!is.na(s2013) & is.na(s2012), 1, s2012),
           s2013 = ifelse(!is.na(s2014) & is.na(s2013), 1, s2013),
           s2014 = ifelse(!is.na(s2015) & is.na(s2014), 1, s2014),
           s2015 = ifelse(!is.na(s2016) & is.na(s2015), 1, s2015),
           s2016 = ifelse(!is.na(s2017) & is.na(s2016), 1, s2016),
           s2017 = ifelse(!is.na(s2018) & is.na(s2017), 1, s2017),
           s2018 = ifelse(!is.na(s2019) & is.na(s2018), 1, s2018),
           s2019 = ifelse(!is.na(s2020) & is.na(s2019), 1, s2019),
           s2020 = ifelse(!is.na(s2021) & is.na(s2020), 1, s2020),
           s2021 = ifelse(!is.na(s2022) & is.na(s2021), 1, s2021),
           s2022 = ifelse(!is.na(s2023) & is.na(s2022), 1, s2022),
           s2023 = ifelse(!is.na(s2024) & is.na(s2023), 1, s2023)) %>% 
    mutate(ID = 999)
  # figure out how to do this w a for loop someday
  
  # change values before first obs to 999
  for(i in 1:nrow(eh)){
    eh[i, 1:(min(which(eh[i,] == 1)) -1)] <- 999
  }
  
  eh <- eh %>% select(2:18) # fix later
  
  # OBSERVATION STATES
  # 1 - seen on site
  # 2 - seen off site
  # 3 - recovered roadkill
  # 4 - recovered other
  # 5 - undetected
  
  # gather ID & which were found dead
  id <- left_join(id, surv) %>% select(ID, Dead)
  id <- id %>% mutate(Dead = ifelse(!is.na(Dead), 1, 0))
  
  eh <- cbind(id, eh)
  
  # create a Fate column
  # Fate = 3 when found dead of human cause
  # Fate = 4 when found dead of natural cause
  eh <- eh %>%
    rowwise() %>%
    mutate(Fate = ifelse(any(c_across(3:19) == 2), 3, NA)) %>%
    mutate(Fate = ifelse(is.na(Fate) & Dead == 1, 4, Fate))
  
  # SPECIAL CASES: change Dead to 1 for handful of RKs
  # that were not recovered according to survival file
  eh <- eh %>% mutate(Dead = ifelse(!is.na(Fate), 1, Dead))
  id <- id %>% select(ID) %>% left_join(.,eh) %>% select(ID, Dead)
  
  # now that cause of death is stored in the Fate column
  # replace observations of 0, 2, & 3 with observation state 5
  # currently 2s are at first October when inds were not alive to be seen
  eh <- eh %>%
    mutate_at(3:19, ~replace(., . == 0, 5)) %>% 
    mutate_at(3:19, ~replace(., . == 2, 5)) %>%
    mutate_at(3:19, ~replace(., . == 3, 2)) %>%
    mutate_at(3:19, ~replace(., . == 4, 5))
  
  # create column Gone with year where each ID is first unobserved
  # value for Gone is NA for individuals still alive at end of 2023
  suppressWarnings(
    eh <- eh %>% 
    rowwise() %>% 
    mutate(Gone = max(which(c_across(3:19) == 5)),
           Gone = ifelse(Gone == "-Inf" | ID == 832, NA, Gone)) %>%
    ungroup()) # suppose to give -Inf warnings
  
  # split eh into alive vs dead IDs
  live <- eh %>% filter(is.na(Gone))
  dead <- eh %>% filter(!is.na(Gone))
  
  tmp <- dead %>% select(1:2)
  dead <- dead %>% select(3:21)
  
  # when Fate is 3 or 4, replace last obs of 1 or 2 with Fate
  for(i in 1:nrow(dead)) {
    last5 <- max(which(dead[i, 1:17] == 5))
    if(last5 > 0 && !is.na(dead[i, 18])) {
      dead[i, last5-1] <- dead[i, 18]
    }
  }
  
  # group alive & dead IDs back into one eh
  dead <- cbind(tmp, dead)
  eh <- rbind(dead, live) %>% arrange(ID)
  remove(i, last5, dead, live, tmp)
  
  # select relevant columns & replace remaining NAs with obs 5
  eh <- eh %>%
    select(-Dead) %>% 
    mutate(across(2:18, ~replace(., is.na(.), 5)))
  
  # rename columns by actual year
  eh <- eh %>% rename("2008" = s2008, "2009" = s2009, "2010" = s2010, "2011" = s2011,
                      "2012" = s2012, "2013" = s2013, "2014" = s2014, "2015" = s2015,
                      "2016" = s2016, "2017" = s2017, "2018" = s2018, "2019" = s2019,
                      "2020" = s2020, "2021" = s2021, "2022" = s2022, "2023" = s2023,
                      "2024" = s2024)
  
  eh <- eh %>% select(1:18)
  # write_csv(eh, "eh.csv")
  
  # extract first & last
  # create vector with occasion of first capture
  get.first <- function(x) min(which(x != 999))
  get.last  <- function(x) max(which(x < 5))
  
  
  ## Individual data -----------------------------------------------------------
  
  # observer days data
  obs <- obs %>%
    select(Date, Month, Day, Year, Time, ID, X, Y, Observer) %>% 
    mutate(Time = format(as.POSIXct(Time), format = "%H:%M")) %>% 
    distinct(Year, Date, Observer) %>% 
    group_by(Year) %>% 
    summarise(n = n()) %>% 
    ungroup()
  
  # TEMP; until I figure out Emily's obs*days
  obs <- rbind(obs, c(2021, 10)) %>% 
    arrange(Year)
  
  # age data
  # age <- surv %>%
  #   select(1,3,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36) %>%
  #   mutate_at(vars(-ID, -Cohort), ~ifelse(.== "A", NA, .)) %>%
  #   mutate_at(vars(-ID, -Cohort), ~as.numeric(.)) %>%
  #   mutate("2024" = as.numeric(NA)) %>%
  #   rename("2008" = "Age08", "2009" = "Age09", "2010" = "Age10",
  #          "2011" = "Age11", "2012" = "Age12", "2013" = "Age13",
  #          "2014" = "Age14", "2015" = "Age15", "2016" = "Age16",
  #          "2017" = "Age17", "2018" = "Age18", "2019" = "Age19",
  #          "2020" = "Age20", "2021" = "Age21", "2022" = "Age22",
  #          "2023" = "Age23")
  # 
  # tmp <- age %>% select(3:19)
  # 
  # for (j in 1:nrow(tmp)) {
  #   # forward fill
  #   for (i in 2:ncol(tmp)) {
  #     if (is.na(tmp[j, i])) {
  #       tmp[j, i] <- tmp[j, i - 1] + 1
  #     }
  #   }
  # 
  #   # backward fill
  #   for (i in (ncol(tmp) - 1):1) {
  #     if (is.na(tmp[j, i])) {
  #       tmp[j, i] <- tmp[j, i + 1] - 1
  #     }
  #   }
  # }
  # 
  # tmp[tmp < 0] <- NA
  # 
  # age <- tmp
  # remove(tmp)
  # 
  # write_csv(age, "ageM.csv")
  
  age <- as.matrix(age)+1
  ageC <- c(0, rep(1,2), rep(2,4), rep(3,3), rep(4,30))+1
  remove(surv)
  

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
    ungroup()
  
  env  <- env[3:19,] # [2008:2024,]
  obs  <- round(as.numeric(scale(obs$n)), 3)
  veg  <- round(as.numeric(scale(env$Veg)), 3)
  dens <- round(as.numeric(scale(env$Dens)), 3)
  win  <- round(as.numeric(scale(env$Win)), 3)
  
  noVeg  <- c(as.numeric(which(is.na(veg))));  nNoVeg  <- length(noVeg)
  noDens <- c(as.numeric(which(is.na(dens)))); nNoDens <- length(noDens)
  noWin  <- c(as.numeric(which(is.na(win))));  nNoWin  <- length(noWin)
  
  # write_csv(env, "env.csv")
  
  
  ## Assemble data list --------------------------------------------------------
  
  y <- eh[,-1] %>% as.matrix() %>% unname()
  
  n.inds <- nrow(id)
  n.ageC <- max(ageC)
  n.occasions <- ncol(y)
  n.obs.states <- 5
  n.true.states <- 5
  
  # missing age
  unkAge <- which(apply(age, 1, function(x){all(is.na(x))}))
  
  y <- y[-unkAge, ]
  id <- id[-unkAge, ]
  age <- age[-unkAge, ]
  n.inds <- nrow(id)
  
  first <- apply(y, 1, get.first)
  last  <- apply(y, 1, get.last)
  
  # assemble list
  dataRK <- list(
    y = y,
    id = id,
    obs = obs,
    env = env,
    age = age,
    ageC = ageC,
    first = first,
    last = last,
    
    veg = veg,
    dens = dens,
    win = win,
    noVeg = noVeg,
    noDens = noDens,
    noWin = noWin,
    nNoVeg = nNoVeg,
    nNoDens = nNoDens,
    nNoWin = nNoWin,
    
    n.inds = n.inds,
    n.ageC = n.ageC,
    n.occasions = n.occasions,
    n.obs.states = n.obs.states,
    n.true.states = n.true.states,
    unkAge = unkAge
    )
  
}

