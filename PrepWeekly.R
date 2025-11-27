# 27 November 2025
# Attempt to sort survival & migration
# from semi-weekly observation data


## Set up ----------------------------------------------------------------------

library(readxl)
library(tidyverse)
library(lubridate)

years <- 10:25

# fetch each year's data
# join into a single list
surv <- lapply(years, function(y) {
  read_excel("data/PromSurvivalNov25_RB.xlsx", sheet = paste0("Survival", y))
})

# join into a single dataset
# names(surv) <- paste0("surv", years)
surv <- reduce(surv, full_join, by = c("ID", "Sex", "Cohort"))
