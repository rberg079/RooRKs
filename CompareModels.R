#' Compare outputs of different models
#'
#' @param nYear integer. Number of years to consider in the analysis. nYear = 17 by default.
#' @param minYear integer. First year to consider in the analysis. minYear = 2008 by default.
#' @param maxYear integer. Last year to consider in the analysis. maxYear = minYear + nYear - 1 by default.
#' @param nAgeC integer. Number of age classes to consider. nAgeC = 5 by default.
#' @param plotAges 
#' @param plotYears 
#' @param postPaths character vector. Paths to .rds files with posterior samples of models to compare.
#' @param modelNames character vector. User-defined names for models to compare. 
#' @param plotFolder character string. Path to the folder in which to store plots.
#' @param returnSumData logical. If TRUE, returns a data frame containing posterior samples from 
#' all compared models as an object in the R global environment. returnSumData = FALSE by default.
#'
#' @returns a data frame of posterior samples from all compared models, provided returnSumData is TRUE.
#' @export
#'
#' @examples

compareModels_ME <- function(nYear = 17,
                             minYear = 2008,
                             maxYear,
                             nAgeC = 5,
                             # plotAges = 1:nAgeC,
                             # plotYears = 1:nYear,
                             postPaths,
                             modelNames,
                             plotFolder,
                             returnSumData = FALSE){
  
  # for testing purposes
  nYear = 17
  minYear = 2008
  nAgeC = 5
  # plotAges = c(2, 6, 10, 14)
  # plotYears = c(2, 6, 10, 14)
  postPaths = c("results/modelF_tObs_aVeg_atMR_obs.rds",
                "results/modelF_tObs_aVeg_atMR.rds")
  modelNames = c("modF_obs",
                 "modF_og")
  plotFolder = c("figures/tweaks")
  returnSumData = TRUE
  
  
  ## Set up --------------------------------------------------------------------
  
  library(coda)
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(data.table))
  library(NatParksPalettes)
  library(paletteer)
  
  # check that models are specified correctly
  if(missing(postPaths)){
    stop("Models must be specified via file paths (postPaths).")
  }
  
  # make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  nModels <- length(modelNames)
  maxYear <- minYear + nYear - 1
  
  plotCols  <- paletteer_c("grDevices::Temps", nModels)
  plotAges  <- 1:nAgeC
  plotYears <- 1:nYear
  
  ## Reformat posterior samples ------------------------------------------------
  
  post.dat <- data.frame()
  for(i in 1:nModels){
    
    # read in RDS as mcmc.list  
    post <- readRDS(postPaths[i])
    
    # convert to matrix then data.table
    samples <- do.call(rbind, lapply(post, as.matrix))
    rownames(samples) <- 1:nrow(samples)  # add row names
    model.dat <- as.data.table(samples, keep.rownames = "Sample")
    
    # reshape to long format
    model.dat <- melt(model.dat,
                      id.vars = "Sample",
                      variable.name = "Parameter",
                      value.name = "Value")
    
    # add identifier & bind
    model.dat[, Model := modelNames[i]]
    post.dat <- rbindlist(list(post.dat, model.dat))
  }
  
  
  ## Posterior summaries -------------------------------------------------------
  
  # summarise samples as median & 95% CrI
  sum.dat <- post.dat %>%
    group_by(Parameter, Model) %>%
    summarise(Median = median(Value, na.rm = TRUE),
              Lower  = quantile(Value, 0.025, na.rm = TRUE),
              Upper  = quantile(Value, 0.975, na.rm = TRUE),
              .groups = "keep") %>%
    mutate(Parameter = as.character(Parameter)) %>%
    ungroup()
  
  
  ## Extract indices -----------------------------------------------------------
  
  idx.dat <- sum.dat %>%
    distinct(Parameter) %>%
    mutate(Idx = str_extract_all(Parameter, "\\d+"),
           Idx1 = map_dbl(Idx, ~ ifelse(length(.x) >= 1, as.numeric(.x[1]), NA)),
           Idx2 = map_dbl(Idx, ~ ifelse(length(.x) >= 2, as.numeric(.x[2]), NA)),
           
           Age = case_when(grepl('mean.phi|mean.M|mean.R|mu.phi|mu.M|mu.R|B.veg', Parameter) ~ Idx1,
                           # Parameter %in% c("mean.phi", "mean.M", "mean.R", "mu.phi", "mu.M", "mu.R", "B.veg") ~ Idx1,
                           TRUE ~ NA_real_),
           
           Year = case_when(grepl('mean.Pi|mean.Po|mean.rR|mean.rO|veg', Parameter) ~ Idx1,
                            grepl('mean.phi|mean.M|mean.R', Parameter) ~ Idx2,
                            # Parameter %in% c("mean.phi", "mean.M", "mean.R") ~ Idx2,
                            # Parameter %in% c("mean.Pi", "mean.Po", "mean.rR", "mean.rO", "veg") ~ Idx1,
                            TRUE ~ NA_real_),
           
           Year = ifelse(!is.na(Year), Year + minYear - 1, NA_real_),
           
           ParamName = word(Parameter, 1, sep = "\\["),
           ParamName = ifelse(ParamName %in% c("mean.phi", "mean.M", "mean.R") & !is.na(Age),
                              paste0(ParamName, "[", Age, "]"), ParamName))
  
  sum.dat <- left_join(sum.dat, idx.dat, by = "Parameter")
  
  
  ## Parameters to plot --------------------------------------------------------
  
  # posterior density overlaps
  plot.D <- list(
    
    Survival = c(paste0("mu.phi[", plotAges, "]"),
                 paste0("B.veg[", plotAges, "]"),
                 "sigma.phi"),

    Movement = c(paste0("mu.M[", plotAges, "]"),
                 "sigma.M"),

    Roadkill = c(paste0("mu.R[", plotAges, "]"),
                 "sigma.R"),

    Observation = c("mu.rR", "mu.rO",
                    "B.obsR", "B.obsO",
                    "sigma.rR", "sigma.rO"))
  
  
  # time series of vital rates
  plot.TS <- list(
    
    ParamNames = c(expand.grid(a = plotAges) %>%
                     mutate(param = paste0('mean.phi[', a, ']')) %>%
                     pull(param),
                   expand.grid(a = plotAges) %>%
                     mutate(param = paste0('mean.M[', a, ']')) %>%
                     pull(param),
                   expand.grid(a = plotAges) %>%
                     mutate(param = paste0('mean.R[', a, ']')) %>%
                     pull(param),
                   'mean.Pi',
                   'mean.Po',
                   'mean.rR',
                   'mean.rO'
                   # params = paste0('mean.Pi[', plotYears, ']'),
                   # params = paste0('mean.Po[', plotYears, ']'),
                   # params = paste0('mean.rR[', plotYears, ']'),
                   # params = paste0('mean.rO[', plotYears, ']')
                   ),

    ParamLabels = c('Survival probability (YAF)',
                    'Survival probability (1-2)',
                    'Survival probability (3-6)',
                    'Survival probability (7-9)',
                    'Survival probability (10+)',
                    'Movement probability (YAF)',
                    'Movement probability (1-2)',
                    'Movement probability (3-6)',
                    'Movement probability (7-9)',
                    'Movement probability (10+)', 
                    'Roadkill probability (YAF)',
                    'Roadkill probability (1-2)',
                    'Roadkill probability (3-6)',
                    'Roadkill probability (7-9)',
                    'Roadkill probability (10+)',
                    'Detection (on-site)',
                    'Detection (off-site)',
                    'Recovery (roadkill)',
                    'Recovery (other)')
  
    
    # mean.phi = list(name = "Survival probability",
    #                 params = expand.grid(a = plotAges, t = plotYears) %>%
    #                   mutate(p = paste0("mean.phi[", a, ",", t, "]")) %>%
    #                   pull(p)),
    # 
    # mean.M = list(name = "Movement probability",
    #               params = expand.grid(a = plotAges, t = plotYears) %>%
    #                 mutate(p = paste0("mean.M[", a, ",", t, "]")) %>%
    #                 pull(p)),
    # 
    # mean.R = list(name = "Roadkill probability",
    #               params = expand.grid(a = plotAges, t = plotYears) %>%
    #                 mutate(p = paste0("mean.R[", a, ",", t, "]")) %>%
    #                 pull(p)),
    # 
    # mean.Pi = list(name = "Detection (on-site)",
    #                params = paste0("mean.Pi[", plotYears, "]")),
    # 
    # mean.Po = list(name = "Detection (off-site)",
    #                params = paste0("mean.Po[", plotYears, "]")),
    # 
    # mean.rR = list(name = "Recovery (roadkill)",
    #                params = paste0("mean.rR[", plotYears, "]")),
    # 
    # mean.rO = list(name = "Recovery (other)",
    #                params = paste0("mean.rO[", plotYears, "]"))
    
    )
  
  
  ## Plot ----------------------------------------------------------------------
  
  # posterior overlaps
  pdf(file.path(plotFolder, "/PostDensities.pdf"), width = 8, height = 4)
  
  for(x in 1:length(plot.D)){
    print(
      ggplot(subset(post.dat, Parameter %in% plot.D[[x]]),
             aes(x = Value, color = Model, fill = Model)) +
        geom_density(alpha = 1 / nModels) +
        facet_wrap(~Parameter, scales = "free") +
        scale_fill_manual(values = plotCols) +
        scale_color_manual(values = plotCols) +
        theme_bw() +
        theme(panel.grid = element_blank())
    )
  }
  dev.off()
  
  # time series
  pdf(file.path(plotFolder, "PostTimeSeries.pdf"), width = 8, height = 4)
  
  for(x in 1:length(plot.TS$ParamNames)){
    print(
      ggplot(subset(sum.dat,
                    ParamName == plot.TS$ParamNames[x] &
                      Year >= minYear & Year <= maxYear),
             aes(group = Model)) +
        geom_ribbon(aes(x = Year, ymin = Lower, ymax = Upper, fill = Model), alpha = 1 / nModels) +
        geom_line(aes(x = Year, y = Median, color = Model)) +
        # facet_wrap(~Age) +
        scale_fill_manual(values = plotCols) +
        scale_color_manual(values = plotCols) +
        scale_x_continuous(breaks = c(minYear:maxYear)) +
        ggtitle(plot.TS$ParamLabels[x]) +
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major.y = element_blank(),
              axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
  }
  dev.off()
  
  # optional: return summary data
  if(returnSumData){
    return(sum.dat)
  }
}

