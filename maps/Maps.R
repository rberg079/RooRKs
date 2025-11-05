library(tidyverse)
library(lubridate)
library(readxl)
# library(maps)
library(ggmap)
library(ggthemes)
library(RgoogleMaps)

setwd("C:/Users/rberg/OneDrive/Bureau/uSherbrooke/Analyses/Originals")

sex <- read_excel("PromlistAllNov22.xlsx")
sex <- sex %>% dplyr::select(1,2) %>% rename(ID = I.D.) %>%
  mutate(Sex = ifelse(Sex == "M", 1, ifelse(Sex == "F", 2, NA)))

setwd("C:/Users/rberg/OneDrive/Bureau")

# import obs and convert from AGD66 to WSG84

OBS <- read_excel("Combined_2008-2019_220517 Kelly.xlsx")
obs <- OBS

obs <- left_join(obs, sex)

obs <- obs %>% dplyr::select(2:7,24,9,10) %>% mutate(ttime = format(as.POSIXct(Time), format = "%H:%M"))
obs <- obs %>% dplyr::select(1:4,6:10) %>% rename(Time = ttime)

obs$X <- as.numeric(obs$X)
obs$Y <- as.numeric(obs$Y)

obs$X <- obs$X+400121
obs$Y <- obs$Y+5600177

# convert from WSG84 to lat/lon

library(sp)
library(terra)
# library(raster)

obs <- obs %>% mutate(zone = 55)

points <- obs %>% dplyr::select(X,Y) %>% rename(x = X, y = Y)
v <- vect(points, geom = c("x","y"), crs = "+proj=utm +zone=55 +south +datum=WGS84 +units=m")
v

y <- project(v, "+proj=longlat +datum=WGS84")
y

lonlat <- geom(y, wkt = F, hex = F, df = T)
lonlat

lonlat <- lonlat %>% dplyr::select(x,y)
obs <- cbind(obs, lonlat)

remove(points,v,y)

# attempt to plot

register_google(key = "AIzaSyBmrUfzV5bt369r9xb1bhlDx_d9K_FbiXE")

# Sample df with airstrip location
lat <- c(-38.9457, -38.9544)
lon <- c(146.2731, 146.2931)

# using RgoogleMaps
# center = c(mean(lat), mean(lon))
# map <- GetMap(center = center, zoom = 15, type = "osm", destfile = "map.png")
# PlotOnStaticMap(map) # doesn't work well

# using ggmap
prom_map <- get_map(location = c(lon = mean(lon), lat = mean(lat)),
                    maptype = "satellite",
                    zoom = 15) # markers = lonlat

map <- ggmap(prom_map)
map

# map + geom_point(data = obs, aes(x = x, y = y, colour = ID), size = 1, alpha = 0.5)

# select useful ID*Years

obs <- obs %>% filter(!is.na(ID), ID <= 1400,   # limit to known IDs under 1400
                      !is.na(x), !is.na(y),     # limit to observations with coordinates
                      Month >= 7)               # limit to main field season, July - Dec

# limit to IDs seen at least 10 times in a year

obs <- obs %>% group_by(Year, ID) %>% mutate(DaysObs = n_distinct(Date)) %>% ungroup()
obs <- obs %>% arrange(ID, Year)

obs <- obs %>% filter(DaysObs >= 10)

obs <- obs %>% mutate(ID = as.factor(ID))
obs <- obs %>% mutate(Year = as.factor(Year))

library(ggsn)
library(ggspatial)

# select a few useful ID*Years to plot

tmp <- obs %>% filter(Year == 2012)
tmp <- tmp %>% filter(Sex == 1)

# plot 10 random observations from 10 random IDs

# ids <- sample(unique(tmp$ID), 6)                                # select 10 random kangaroos
# ran <- tmp %>% filter(ID %in% ids)                              # limit data set to these 10
# ran <- ran %>% group_by(ID) %>% sample_n(10, replace = FALSE)   # select 10 random observations

# ran <- tmp %>% filter(ID == 3) # 371, 115, 217, 34, 54, 3
# ran <- tmp %>% filter(ID == 536) # 440, 491, 530, 557, 522
ran1 <- tmp %>% filter(ID == 3)
ran2 <- tmp %>% filter(ID == 289)


mapP <- map +
  geom_point(data = ran1, aes(x = x, y = y), colour = "#E57373", size = 1) +
  geom_point(data = ran2, aes(x = x, y = y), colour = "#00ACC1", size = 1) +
  labs(x = "Longitude", y = "Latitude")


# mapP <- map + geom_point(data = ran, aes(x = x, y = y), colour = "#E57373", size = 1)
mapP

mapP <- mapP +
  geom_point(data = ran, aes(x = x, y = y), colour = "#00ACC1", size = 1) +
  scalebar(x.min = 146.2700, x.max = 146.2950, y.min = -38.95880, y.max = -38.94000,
           dist = 100, dist_unit = "m", transform = TRUE, model = 'WGS84',
           height = 0.01, st.dist = 0.02, st.size = 2, st.color = "white",
           box.fill = c("grey", "white"), border.size = 0.2) +
  labs(x = "Longitude", y = "Latitude")

# north2(mapF, x = 0.688, y = 0.18, scale = 0.05, symbol = 1)

mapP +
  annotation_north_arrow(
    which_north = "true",
    location = "br",
    height = unit(1, "cm"), width = unit(1, "cm"),
    pad_x = unit(1.2, "cm"), pad_y = unit(1.8, "cm"),
    # pad_x = unit(1.45, "cm"), pad_y = unit(2.5, "cm"),
    style = north_arrow_fancy_orienteering(
      fill = c("grey", "white"),
      line_col = "black",
      text_col = "white")
  )

ggsave("map2012_3&289.jpeg", scale = 1, width = 18.0, height = 18.0, units = c("cm"), dpi = 600)


