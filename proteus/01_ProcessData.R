############################################################################## #
# Bayesian known-fate model for deer GPS collar 
# Data processing for GPS data
# Author: Rachel Hickcox, Abby Bratt (Proteus)
# Client: DRNSW
# Year: 2023
############################################################################## #

setwd("C:/Users/Rachel/OneDrive/Bureau/Fallow Deer")
dat <- read.csv("CollarData_allsites_simplified.csv")

#### LOAD LIBRARIES ####
install_and_load_packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) { install.packages(new.pkg, dependencies = TRUE)}
  lapply(pkg, library, character.only = TRUE)
}
install_and_load_packages(c("tidyverse", "here", "lubridate"))

#### LOAD DATA #### 
# dat <- readRDS(here("Data", "CollarData_allsites_simplified.RDS"))
# dat <- read.csv(here("CollarData_allsites_simplified.csv"))

#### PROCESS COVARIATE DATA ####
# site (site)
# sex (sex)
# collar manufacturer (mfr)
X <- dat %>% 
  select(deer_number, site, sex, mfr) %>% 
  distinct() %>% 
  arrange(deer_number)

unique(dat$fate)

#### INDIVIDUAL FATE AT LAST OBSERVATION ####
# status` is the status of the animal at the time the data were downloaded
# 'fate` is the ultimate fate of each animal, as far as is known. There is one `dead unknown`. This animal died but couldn't be reached until a few days after the mortality signal was triggered, by which time feral pigs had consumed most of the carcass. `collar fail unknown` indicates a collar that stopped transmitting without any warning.
# note this is at individual level, not the fix level

# OBSERVATION STATES
# 1 - alive and active collar
# 2 - collar failed
# 3 - newly caused natural mortality + active collar
# 4 - newly caused human mortality + active collar
# 5 - undetected

fate <- dat %>% 
  select(deer_number, status, fate) %>% 
  distinct() %>% 
  mutate(fate = case_when(
    status == "active" ~ "active", # status = active, fate = active
    TRUE ~ fate)) %>% # status = dead, fate = dead
  mutate(code = case_when( # TODO check these codes
    fate == "active" ~ 1, 
    fate == "dead accident" ~ 3, # assume this was natural
    fate == "dead unknown" ~ 3, # assume this was natural if not obviously shot
    fate == "dead shot aerial" ~ 4, 
    fate == "dead shot contractor" ~ 4,
    fate == "dead shot hunters" ~ 4, 
    fate == "dead shot landholder" ~ 4,
    fate == "collar fail unknown" ~ 5,
    TRUE ~ 2 # all others are inactive, assumed alive
  )) %>% 
  arrange(deer_number)

table(fate$code)

# TRANSITION STATES
# 1 - alive and active collar
# 2 - alive and inactive collar, observed
# 3 - alive and inactive collar, undetected
# 4 - newly caused natural mortality + active collar
# 5 - newly caused natural mortality + inactive collar
# 6 - newly caused human mortality + active collar
# 7 - newly caused human mortality + inactive collar
# 8 - long dead

#### CLEAN DATA ####
# use deer no as unique id
# `dtg_fix` is the local date and time of the fix.
#  `cday` is the number of days since the deer was collared
# fix number is not right for some reason
dat_clean <- dat %>% 
  select(-c(site, sex, mfr, # remove covariates
            status, fate, # remove fates
            collar_serial_number, # not necessary
            latitude, longitude, hdop, num_sats, # spatial stuff
            dtg_fix_sec # not necessary
  )) %>% 
  arrange(deer_number, dtg_fix) %>% 
  mutate(fix_date = date(dtg_fix), # compute day and week of fix 
         fix_week = week(dtg_fix),
         fix_week = if_else(fix_week == 53, 52, fix_week), # slight correction for leap years
         fix_week = parse_date_time(paste(year(dtg_fix), fix_week, 1, sep = "/"),'Y/W/w'),
         fix_month = month(dtg_fix), 
         fix_month = parse_date_time(paste(year(dtg_fix), fix_month, 1, sep = "/"), "Y/m/d")
  ) %>% 
  group_by(deer_number) %>% 
  mutate(first_fix = min(dtg_fix) %>% ymd_hm(), # compute first and last fix time for detection histories
         collar_time = first_fix - hours(min(fix_num)),
         collar_day = date(collar_time), 
         first_day = collar_day, 
         first_week = min(fix_week), 
         first_month = min(fix_month),
         last_fix = max(dtg_fix) %>% ymd_hm(), # RACHEL EDIT: ymd_hm() was as_datetime() but was not working
         last_day = date(last_fix),
         last_week = max(fix_week), 
         last_month = max(fix_month)
  ) %>%  
  ungroup()

##### Create update fate data frame with timestamps #####
fate_clean <- dat_clean %>% 
  group_by(deer_number) %>% 
  select(first_fix, first_day, first_week, first_month,
         last_fix,  last_day,  last_week,  last_month) %>% 
  distinct() %>% 
  left_join(fate, by = "deer_number")

##### First day/week and last day/week to create detection history matrix that is the right size #####
min_first_daily <- min(dat_clean$collar_day)
max_last_daily <- max(dat_clean$last_day) 
min_first_weekly <- min(dat_clean$fix_week) 
max_last_weekly <- max(dat_clean$fix_week)
min_first_monthly <- min(dat_clean$fix_month)
max_last_monthly <- max(dat_clean$fix_month)

seq_daily <- seq(min_first_daily, max_last_daily, by = "day")
seq_weekly <- seq(min_first_weekly, max_last_weekly, by = "week")
seq_monthly <- seq(min_first_monthly, max_last_monthly, by = "month")

#### CLEAN(ER) DATA - DAILY SCALE ####
dat_clean_daily <- dat_clean %>% 
  group_by(deer_number, fix_date) %>% 
  summarise(n_fixes = n(), # within a day, how many fixes for this individual
            detected = n_fixes > 0) %>% # was the individual detected at least once this day
  complete(fix_date = seq_daily, # fill in missing days for this individual in the time series
           fill = list(n_fixes = 0, detected = FALSE)) %>% 
  mutate(detected = as.numeric(detected)) %>% 
  ungroup() %>% 
  left_join(fate_clean %>% 
              select(deer_number, first_day, last_day, 
                     status, fate, code), 
            by = "deer_number") %>% 
  mutate(obs_code = case_when(
    (detected == 0 & fix_date < first_day) ~ 999, # prior to collaring
    (fix_date == first_day) ~ 1, # active at collaring
    (detected == 1 & fix_date != last_day) ~ 1, # collar detected and active 
    (fix_date == last_day) ~ code, # day of last observation, fill in fate 
    (detected == 0 & fix_date > last_day) ~ 5,
    TRUE ~ 5))

#### DAILY ENCOUNTER HISTORY ####
# Create encounter history, first column is deer ID
eh_daily <- dat_clean_daily %>% 
  select(deer_number, fix_date, obs_code) %>% 
  pivot_wider(id_cols = deer_number, names_from = fix_date, values_from = obs_code)

#### CLEAN(ER) DATA - WEEKLY SCALE ####
dat_clean_weekly <- dat_clean %>% 
  group_by(deer_number, fix_week) %>% 
  summarise(n_fixes = n(), # within a week, how many fixes for this individual
            detected = n_fixes > 0) %>% # was the individual detected at least once this week
  complete(fix_week = seq_weekly, # fill in missing weeks for this individual in the time series
           fill = list(n_fixes = 0, detected = FALSE)) %>% 
  mutate(detected = as.numeric(detected)) %>% 
  ungroup() %>% 
  left_join(fate_clean %>% 
              select(deer_number, first_week, last_week, 
                     status, fate, code), 
            by = "deer_number") %>% 
  mutate(obs_code = case_when(
    (detected == 0 & fix_week < first_week) ~ 999, # prior to collaring
    (fix_week == first_week) ~ 1, # active at collaring
    (detected == 1 & fix_week != last_week) ~ 1, # collar detected and active 
    (fix_week == last_week) ~ code, # day of last observation, fill in fate 
    (detected == 0 & fix_week > last_week) ~ 5,
    TRUE ~ 5))

#### WEEKLY ENCOUNTER HISTORY ####
# Create encounter history, first column is deer ID
eh_weekly <- dat_clean_weekly %>% 
  select(deer_number, fix_week, obs_code) %>% 
  pivot_wider(id_cols = deer_number, names_from = fix_week, values_from = obs_code)

#### CLEAN(ER) DATA - MONTHLY SCALE ####
dat_clean_monthly <- dat_clean %>% 
  group_by(deer_number, fix_month) %>% 
  summarise(n_fixes = n(), # within a week, how many fixes for this individual
            detected = n_fixes > 0) %>% # was the individual detected at least once this week
  complete(fix_month = seq_monthly, # fill in missing weeks for this individual in the time series
           fill = list(n_fixes = 0, detected = FALSE)) %>% 
  mutate(detected = as.numeric(detected)) %>% 
  ungroup() %>% 
  left_join(fate_clean %>% 
              select(deer_number, first_month, last_month, 
                     status, fate, code), 
            by = "deer_number") %>% 
  mutate(obs_code = case_when(
    (detected == 0 & fix_month < first_month) ~ 999, # prior to collaring
    (fix_month == first_month) ~ 1, # active at collaring
    (detected == 1 & fix_month != last_month) ~ 1, # collar detected and active 
    (fix_month == last_month) ~ code, # day of last observation, fill in fate 
    (detected == 0 & fix_month > last_month) ~ 5,
    TRUE ~ 5))

#### MONTHLY ENCOUNTER HISTORY ####
# Create encounter history, first column is deer ID
eh_monthly <- dat_clean_monthly %>% 
  select(deer_number, fix_month, obs_code) %>% 
  pivot_wider(id_cols = deer_number, names_from = fix_month, values_from = obs_code)

#### SAVING DATA ####
# write.csv(eh_weekly, here("weekly_data_processed.csv"))
# write.csv(eh_monthly, here("monthly_data_processed.csv"))
# write.csv(fate_clean, here("fates.csv"))

#### NECESSARY CONSTANTS FOR MODEL ####
inds <- eh_daily$deer_number
n_inds <- length(inds)
# Extract first and last, Create vector with occasion of collaring
get.first <- function(x) min(which(x != 999))
first_daily <- apply(eh_daily[, -1], 1, get.first)
first_weekly <- apply(eh_weekly[, -1], 1, get.first)
first_monthly <- apply(eh_monthly[, -1], 1, get.first)
# Create vector with last ping
get.last <- function(x) max(which(x != 7))
last_daily <- apply(eh_daily[, -1], 1, get.last)
last_weekly <- apply(eh_weekly[, -1], 1, get.last)
last_monthly <- apply(eh_monthly[, -1], 1, get.last)

