
# Outline of preperation of the bat data for analysis. This part involves organising sample sites into spatial groupings for the measurement of spatial variance.

#Load in the necessary libraries

library(ggplot2)
library(dplyr)
library(sp)
library(dismo)
library(rgeos)
library(tidyverse)

# Read in the data:

dd <- read_csv("Data/Edited_all_files_with_lat_long_VKedits.csv")
str(dd)

dd <- dd %>%
  dplyr::select(Lat, Long, TrapName, Site) %>%
  rename(Spatial.replicate = TrapName) %>%
  filter(!is.na(Lat),
         !is.na(Long))



# Some of the names for Spatial replicates begin with numbers so we will just rename them to avoid problems later on.

dd$Spatial.replicate <- paste0("P_", dd$Spatial.replicate)
dd$Spatial.replicate <- gsub("-", ".", dd$Spatial.replicate)
dd$Spatial.replicate <- as.factor(dd$Spatial.replicate)


#Lets strip out the species information and save it as a new site.info dataframe.

site.info <- dd

# remove duplicates
site.info <- site.info[!duplicated(site.info), ]


#We have between 21 and 87 sampling points per site

#tapply(site.info$Spatial.replicate, site.info$Replicate, length)



# Lets focus in on site DVCA. 

site.DVCA <- site.info[]

# I think this projection is right
safe.proj <- CRS("+proj=longlat +datum=WGS84")

coords <- SpatialPointsDataFrame(data = data.frame(Point = site.DVCA$Spatial.replicate),
                                 SpatialPoints(cbind(site.DVCA$Long, site.DVCA$Lat), 
                                               proj4string = safe.proj))




dists <- matrix(spDists(coords, longlat = TRUE, diagonal = FALSE) * 1000, 
                ncol = nrow(site.DVCA),
                dimnames = list(site.DVCA$Spatial.replicate, site.DVCA$Spatial.replicate))



dists[dists == 0] <- NA



# in a loop get the smallest non-na value for every point in the dists matrix

# Make a vector for our output
nearest_vec <- rep(NA, nrow(dists))
# now get the nearest neighbours
for(i in 1:nrow(dists)){
  #if(i == nrow(dists)){break()}
  nearest_vec[i] <- min(dists[i, -which(is.na(dists[,i]))])
}

# the mean is
mean(nearest_vec)

# the sd is 
sd(nearest_vec)
