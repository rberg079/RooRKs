# 7 August 2023

library(readxl)
library(tidyverse)
library(lubridate)

# setwd("C:/Users/Rachel/OneDrive/Bureau/RK analyses")


## Load & clean up data --------------------------------------------------------

SURV <- read_excel("PromSurvivalOct23.xlsx", sheet = "YEARLY SURV")
surv <- SURV

# modify column names to read as survived to s20XX
surv <- surv %>%
  select(1,7:11,12:39) %>%
  rename(s2009 = "08to09", s2010 = "09to10", s2011 = "10to11",
         s2012 = "11to12", s2013 = "12to13", s2014 = "13to14",
         s2015 = "14to15", s2016 = "15to16", s2017 = "16to17",
         s2018 = "17to18", s2019 = "18to19", s2020 = "19to20",
         s2021 = "20to21", s2022 = "21to22", Dead = "Found dead")

# select only usable cohorts & one sex at a time
surv <- surv %>% filter(Cohort < 2022 | is.na(Cohort))
surv <- surv %>% filter(ID <= 1292)

surv <- surv %>% filter(Sex == 2) # females


## Encounter history -----------------------------------------------------------

eh <- surv %>%
  select(1,8,10,12,14,16,18,20,22,24,26,28,30,32,34) %>%
  mutate(s2008 = NA) %>%
  select(1,16,2:15)

id <- eh %>% select(1)

# modify column of previous year
# if roo survived to 2010 then it was alive in 2009
# for(year in 2008:2021){
#   eh <- eh %>%
#     mutate(across(starts_with(paste0("s", year)),
#            ~ifelse(!is.na(get(paste0("s", year + 1))) & is.na(.), 1, .)))
# }

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
         s2021 = ifelse(!is.na(s2022) & is.na(s2021), 1, s2021))
# figure out how to do this w a for loop someday

# change values before first obs to 999
for(i in 1:nrow(eh)){
  eh[i, 1:(min(which(eh[i,] == 1)) -1)] <- 999
}

eh <- eh %>% select(2:16) # fix later

# OBSERVATION STATES
# 1 - seen
# 2 - recovered roadkill
# 3 - recovered other
# 4 - undetected

# gather ID & which were found dead
id <- left_join(id, surv) %>% select(ID, Dead)
id <- id %>% mutate(Dead = ifelse(!is.na(Dead), 1, 0))

eh <- cbind(id, eh)

# create a Fate column
# Fate = 2 when found dead of human cause
# Fate = 3 when found dead of natural cause
eh <- eh %>%
  rowwise() %>%
  mutate(Fate = ifelse(any(c_across(3:17) == 2), 2, NA)) %>%
  mutate(Fate = ifelse(is.na(Fate) & Dead == 1, 3, Fate))

# SPECIAL CASES: change Fate of IDs 58, 88, 363 to NA
# were off the study area before found dead so treat as not found
eh <- eh %>% mutate(Fate = ifelse(ID == 58 | ID == 88 | ID == 363, NA, Fate))

# SPECIAL CASES: change Dead to 1 for handful of RKs
# that were not recovered according to survival file
eh <- eh %>% mutate(Dead = ifelse(!is.na(Fate), 1, Dead))
id <- id %>% select(ID) %>% left_join(.,eh) %>% select(ID, Dead)

# now that cause of death is stored in the Fate column
# replace observations of 0, 2, & 3 with observation state 4
# currently 2s are at first October when inds were not alive to be seen
eh <- eh %>%
  select(3:18) %>%
  mutate_at(1:15, ~replace(., . == 2, 4)) %>%
  mutate_at(1:15, ~replace(., . == 3, 4)) %>%
  mutate_at(1:15, ~replace(., . == 0, 4))

# when Fate is 2 or 3, replace last obs of 1 with Fate
for(i in 1:nrow(eh)) {
  last1 <- max(which(eh[i, 1:15] == 1))
  if(last1 > 0 && !is.na(eh[i, 16])) {
    eh[i, last1] <- eh[i, 16]
  }
}

# select relevant columns & replace remaining NAs with obs state 4
eh <- eh %>% select(1:15) %>% mutate_all(~replace(., is.na(.), 4))
eh <- eh %>% rename("2008" = s2008, "2009" = s2009, "2010" = s2010, "2011" = s2011,
                    "2012" = s2012, "2013" = s2013, "2014" = s2014, "2015" = s2015,
                    "2016" = s2016, "2017" = s2017, "2018" = s2018, "2019" = s2019,
                    "2020" = s2020, "2021" = s2021, "2022" = s2022)

eh <- cbind(id, eh) %>% select(1,3:17)

# write_csv(eh, "eh.csv")

# extract first & last
# create vector with occasion of first capture
get.first <- function(x) min(which(x != 999))
get.last  <- function(x) max(which(x < 4))
first <- apply(eh[, -1], 1, get.first)
last  <- apply(eh[, -1], 1, get.last)

# REMOVE PROBLEM IDs
# which were only seen once
# therefore have no 1s in the eh
id <- cbind(id, first, last)
bad = (id$first >= id$last)
which(bad)

eh <- eh[!bad,]
id <- id[!bad,]


## Constants -------------------------------------------------------------------

inds <- id$ID
n_inds <- length(inds)

first <- as.numeric(apply(eh[, -1], 1, get.first))
last  <- as.numeric(apply(eh[, -1], 1, get.last))

# create vector for year of first obs of each roo
# first = apply(eh, 1, function(i) min(which(i == 1)))
# last  = apply(eh, 1, function(i) max(which(i == 1)))

