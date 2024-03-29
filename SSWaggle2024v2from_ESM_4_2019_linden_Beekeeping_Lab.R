######## This script has been adapted for a new dataset by Douglas Sponsler and Reed Johnson (johnson.5005@osu.edu). 
######## Unused portions of the original script have been commented out but not deleted. Additions are followed by the tag # Sponsler:

# Electronic Supplementary Material 4 - R script to simulate dances from known waggle
# dance durations and headings.
# ------------------------------------------------------------------------------------

# Article title: Incorporating variability in honey bee waggle dance decoding improves
# the mapping of communicated resource locations

# Journal: Journal of Comparative Physiology A

# Authors: Roger Schurch, Margaret J. Couvillon, Dominic D. R. Burns, Kiah
# Tasman, David Waxman and Francis L. W. Ratnieks

# Corresponding author: Roger Schurch, Evolution, Behaviour and
# Environment, School of Life Sciences, University of Sussex, Brighton,
# BN1 9QG, United Kingdom, R.Schuerch@sussex.ac.uk

# script files that go with this script:
# ESM_3.jag

# data files that go with this script:
# ESM_5.csv

# the file will output a comma separated value file and an ASC raster file
# you must create a "data" folder within the folder from which you are running the scripts
# or adapt the paths for the instructions below

# Load Ubuntu packages
# sudo apt install jags r-cran-rgdal r-recommended 
# Load CRAN packages needed in this script
# install.packages(c("circular", "rjags", "sp", "rgdal", "raster", "png", "googlesheets", "magrittr", "oce"))

## Clear everything
rm(list=ls())
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg
library('googlesheets4')  # Read from GoogleDocs
library('magrittr')        # Pipes
library('oce')  # Calculation of azimuth

## Set local working directory (uncomment and change)
setwd("/Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory within the Git repository
## Create data subfolder if it doesn't already exist
dir.create("data", showWarnings = FALSE)

## Import dance data file
gs_title("SureshDance_Thesis") %>%
  gs_download(ws = "Data", to = "SureshDance_Thesis.csv", overwrite = TRUE)


#waggleFile <- "SureshDance_Thesis.csv" # Set the name for the file to use
## Download data from GoogleSheets using 'googlesheets' library.  Store as csv file specified in "wagglefile"
# "SD Tracks" https://docs.google.com/spreadsheets/d/1ZXxYSRK74wCT7xULdF4hyReSCYsl2kUqe04kiROzPUs/edit?usp=sharing
#gs_title("SureshDance_Thesis.csv") %>%
#  gs_download(ws = "Data", to = waggleFile, overwrite = TRUE)
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data

## Remove flagged lines
waggleData <- subset(waggleData, is.na(waggleData$flag)) # Sponsler: a flag field removes empty or incomplete lines

## Fix dates/times and use them to calculate azimuth
waggleData$date <- paste(waggleData$year, sprintf("%02d", waggleData$month), sprintf("%02d", waggleData$day), sep="-")
waggleData$time <- paste(sprintf("%02d", waggleData$hour), sprintf("%02d", waggleData$min), sep=":")
waggleData$dateTime <- as.POSIXct(strptime(paste(waggleData$date, waggleData$time, sep=" "), "%Y-%m-%d %H:%M"), tz="US/Central")
attr(waggleData$dateTime, "tzone") <- "UTC" ## Convert to UTC 
waggleData$azimuth <- sunAngle(waggleData$dateTime, waggleData$lon, waggleData$lat, useRefraction=TRUE)$azimuth

## Calculate heading in radians
waggleData[is.na(waggleData$skew),]$skew <- 0
waggleData$heading.degrees <- (((waggleData$mean_angle-waggleData$skew+waggleData$azimuth)+90) %% 360) ## Something not concordant with GoogleDoc calculation here
waggleData$heading.radians <- (waggleData$heading.degrees*pi)/180

## Designate Bloom vs. Non-Bloom dates
#waggleData$bloom <- FALSE
#waggleData[waggleData$dateTime >= "2016-06-27" & waggleData$dateTime <= "2016-08-04",]$bloom <- TRUE

## Get summary of dances recorded by hive and by date
table(waggleData[,c("hive", "date")])

## Synthetic variable with hive and date
#waggleData$hiveDate <- paste(waggleData$hive, waggleData$date, sep=".")

##################################
## Set calibration
#################################

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

## make the data properly circular
waggleData$heading <- circular(waggleData$heading.radians,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

## how many samples per dance; think carefully about how many samples you really need:
## the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

## preparations to calculate point coords from angle and distance (get from https://epsg.io/)
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 422363.59 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4514220.26 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

## to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

## prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


## prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

## create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

## only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

## prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

#############################################
## Process 
##################################################
## All together

###########################
## Process by date
###########################
#datePick <- "2018-06-25"
#datePick <- "2018-06-28"
datePick <- "2018-07-12"

waggleDataDate <- subset(waggleData, date == datePick)
waggleDataDate 

##########################
## Process by date and color for collected bees
############################
#waggleCollected <- waggleData[waggleData$color != "",]
#waggleCollected$dateColor <- paste(waggleCollected$date, waggleCollected$color,sep=".")
waggleCollected <- waggleData
waggleCollected$dateColor <- paste(waggleCollected$date)

# loop through all dateColors
for (j in unique(waggleCollected$dateColor)) {
  waggleDataDate <- waggleCollected[waggleCollected$dateColor == j,]
  # loop through all the dances

  for(i in 1:length(waggleDataDate$dancer.id)){
    cat(paste(i, "of", length(waggleDataDate$dancer.id), "\n"))
    # choose only the i^th dance
    tempData <- waggleDataDate[i,]

    # prepare the variables for the prediction model
    N2 <- length(tempData$mean_duration.sec)
    x2 <- rep(NA, length(tempData$mean_duration.sec))
    y2 <- tempData$mean_duration.sec
    
    # load the model from file and submit the data
    jags <- jags.model('ESM_3.jag',
                       data = list('x' = x, 'y' = y,
                                   'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                       n.chains = 1,
                       n.adapt = 100)
    
    # update for the burn-in
    update(jags, 100000)
    
    # sample from the posterior
    samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)
    
    # save the samples in a handy variable
    sim.distances <-  samples[,'x2'][[1]]
    
    # the 1000 draws have to be taken according to what is in the posterior samples for distance
    sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                             kappa = 24.9, control.circular = list("radians"))
    
    # calculate the coordinates from the vector with origin of the hives
    rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
    rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)
    
    # save as points for further use
    temp.points <- data.frame(cbind(id = rep(tempData$id, length(rel.dance.easting)),
                                    easting = as.numeric(rel.dance.easting),
                                    northing = as.numeric(rel.dance.northing)))
    
    # save the points in a comma seperated value file
    # csv points can be imported into GIS for further processing
    write.csv(temp.points, paste("data/sim.dance_", tempData$id, ".csv", sep = ""), row.names = FALSE)
    
    if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
      write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
    }else{
      # save the total counts
      tempData2 <- read.csv("data/simAllDances.csv")
      tempData2 <- rbind(tempData2, temp.points)
      write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
    }
    
    # georeference the points on the UK grid
    coordinates(temp.points) = c("easting", "northing")
    #proj4string(temp.points) = CRS("+init=epsg:27700")
    proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)
    
    # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
    # e.g. probability that a dance has been on a certain raster square
    #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
    temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:
    
    # convert the raster to a SpatialGridDataFrame
    #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
    g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:
    
    # save the file to disk (can be imported into GIS)
    currentFileName <- paste("data/raster_", tempData$id, ".asc", sep = "")
    write.asciigrid(g, currentFileName)
    
    if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
      #total.temp.rast <- temp.rast.UKGRID
      total.temp.rast <- temp.rast.UTM17N # Sponsler:
    }else{
      # calculate the probability that a field has been visited
      #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
      total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
    }
  }
  
  # save the combined dances as one raster file ready to be imported in ArcGIS
#  total.temp.rast <- total.temp.rast$id
#  g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
#  write.asciigrid(g.total, "data/totalRaster.asc")
  
  
  # for plotting, we can crop the extent to our needs
  
  # prepare an georeferenced extent (in GIS terms) for the cropping
  coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                                  hiveEasting + 10000),
                                      northing = c(hiveNorthing, hiveNorthing - 10000,
                                                   hiveNorthing + 10000)))
  coordinates(coordPanel.crop) <- c("easting", "northing")
  #proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
  proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:
  
  
  # prepare a raster with the extent to crop the aerial photography
  numCoordPanel.crop <- as.data.frame(coordPanel.crop)
  crop.rast <- raster(ncols = noCells, nrows = noCells)
  extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
  #proj4string(crop.rast) = CRS("+init=epsg:27700")
  proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:
  
  # we crop the data raster to size
  new.data.rast <- crop(total.temp.rast, crop.rast)
  writeRaster(new.data.rast, filename = paste("data/Wooster_Linden_",j,".tif", sep=""), format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer
#  writeRaster(new.data.rast, filename = paste("data/Wooster_Linden_",datePick,".tif", sep=""), format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer  
} # End of by day

#########################################
#########################################

### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)







############################
##### BROKEN UP BY DAY #####
############################

##### MAY 4 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 4) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may4.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 5 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("/Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 5) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may5.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 6 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 6) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may6.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 7 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 7) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may7.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 8 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 8) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may8.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 10 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 10) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may10.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 14 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 14) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may14.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())

##### MAY 24 ###############################################################################################################

# loading packages needed in this script
library('circular')	# for circular stats
library('rjags')	# interface with JAGS
library('sp')		# spatial stats, coordinates etc
library('rgdal')	# convert between different coodrinate systems, spTransform
library('raster')       # for plotting spatial data
library('png')          # to save figures to jpeg

setwd("Users/dolezallab-macpro/Desktop/WaggleDance") # Sponsler: path to the working directory on my laptop

#read the calibration data from ESM_5.csv
calibDataAgg <- read.csv("ESM_5.csv", row.names = 1)
calibDataAgg$heading <- circular(calibDataAgg$heading,
                                 type = "angle",
								 unit = "radian",
								 rotation = "clock",
								 zero = pi/2)

# prepare your own data, i.e. create or read from file a
# data frame with "duration" in sec and "heading" in radians
# DO ONE DANCE AT A TIME, else the prior will overwhelm the data

# as an example, we will use our own data for dances that have gone to
# feeder at 1 km, using only the first dance of every bee.

# dance ids to that feeder
#dance.ids <- c(237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
 #              247, 248, 249, 250, 251, 252, 253, 254, 255, 256,
  #             257, 258, 259, 260, 261, 262, 263, 264, 265, 266,
   #            267, 285, 286, 268, 269, 270)

# store the subset of dances going to that feeder in waggleData
#waggleData <- calibDataAgg[calibDataAgg$dance.id %in% dance.ids,]
waggleData <- read.csv("SureshDance_Thesis.csv") # Sponsler: path to our dance data
waggleData <- subset(waggleData, flag == 1 & day == 24) # Sponsler: a flag field removes empty or incomplete lines

# we only want the first dance of every individual bee, and we prepare a function to achieve that
#getFirstElement <- function(pVector){
	#reslt <- NA
	#if(length(pVector) > 0){
  	  #reslt <- pVector[1]
  	#}
  #reslt
  #}

# now only select the first dance of every bee
#waggleData <- aggregate(cbind(dance.id = waggleData$dance.id,
 #                             duration = waggleData$duration,
                             #distance = waggleData$distance,
  #                            heading = waggleData$heading), #by = list(bee.id = waggleData$bee.id),
                       # getFirstElement)

# make the data properly circular
waggleData$heading <- circular(waggleData$heading,
                               type = "angle",
							   unit = "radian",
							   rotation = "clock",
							   zero = pi/2)

# how many samples per dance; think carefully about how many samples you really need:
# the more samples, the longer it will take to simulate your dances
finalSampleSize <- 1000
thinning <- 100
noJagsSamples <- thinning*finalSampleSize

# preparations to calculate point coords from angle and distance
#hiveEasting <- 534939				# the UK grid easting of the hives in meters
hiveEasting <- 292263.656365 # Sponsler: the UTM 17N (EPSG:26917) easting of the hives in meters
#hiveNorthing <- 108900				# the UK grid northing of the hives in meters
hiveNorthing <- 4426271.356893 # Sponsler: the UTM 17N (EPSG:26917) northing of the hives in meters

# to calculate the rasters
distanceToHives <- 10000			# how far should the rasters extend from the hives in meters
gridCellSize <- 25				# grid size in meters
noCells <- 2*distanceToHives/gridCellSize	# the number of cols and rows needed to get grid of meters
defaultMatrix <- matrix(data = 0,		# create a default matrix as a basis for our dance counts
                        ncol = noCells, nrow = noCells)

# prepare an georeferenced extent (in GIS terms); adapt to fit your coordinate system
coordPanel <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - distanceToHives,
                                 hiveEasting + distanceToHives),
                               northing = c(hiveNorthing, hiveNorthing - distanceToHives,
                                 hiveNorthing + distanceToHives)))
coordinates(coordPanel) <- c("easting", "northing")
#proj4string(coordPanel) = CRS("+init=epsg:27700")
proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)


# prepare a raster with the extent; adapt to fit your coordinate system
numCoordPanel <- as.data.frame(coordPanel)
temp.rast <- raster(ncols = noCells, nrows = noCells)
extent(temp.rast) <- extent(c(numCoordPanel[2:3,1], numCoordPanel[2:3,2]))
#proj4string(temp.rast) = CRS("+init=epsg:27700")
proj4string(temp.rast) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

# create a raster that will store the actual probability-visited data
total.temp.rast <- temp.rast

# only select tagged bees for calibration dances
calibDataAggBees <- calibDataAgg[!is.na(calibDataAgg$bee.id),]

# prepare the variables for the calibration model
N1 <- length(calibDataAggBees$duration)
x <- calibDataAggBees$distance
y <- calibDataAggBees$duration

K <- length(unique(calibDataAggBees$bee.id))
bee <- factor(calibDataAggBees$bee.id)

# loop through all the dances
for(i in 1:length(waggleData$dance.id)){
  cat(paste(i, "of", length(waggleData$dance.id), "\n"))
  # choose only the i^th dance
  tempData <- waggleData[i,]

  # prepare the variables for the prediction model
  N2 <- length(tempData$duration)
  x2 <- rep(NA, length(tempData$duration))
  y2 <- tempData$duration

  # load the model from file and submit the data
  jags <- jags.model('ESM_3.jag',
                     data = list('x' = x, 'y' = y,
                       'N1' = N1, 'K' = K, 'bee' = bee, 'N2' = N2, 'x2' = x2, 'y2' = y2),
                     n.chains = 1,
                     n.adapt = 100)

  # update for the burn-in
  update(jags, 100000)

  # sample from the posterior
  samples <- coda.samples(jags, c('x2'), noJagsSamples, thin = thinning)

  # save the samples in a handy variable
  sim.distances <-  samples[,'x2'][[1]]

  # the 1000 draws have to be taken according to what is in the posterior samples for distance
  sim.heading <- rvonmises(finalSampleSize, mu = tempData$heading,
                           kappa = 24.9, control.circular = list("radians"))

  # calculate the coordinates from the vector with origin of the hives
  rel.dance.easting <- as.numeric(hiveEasting + cos(-(sim.heading- pi/2))*sim.distances)
  rel.dance.northing <- as.numeric(hiveNorthing + sin(-(sim.heading - pi/2))*sim.distances)

  # save as points for further use
  temp.points <- data.frame(cbind(dance.id = rep(tempData$dance.id, length(rel.dance.easting)),
                                  easting = as.numeric(rel.dance.easting),
                                  northing = as.numeric(rel.dance.northing)))

  # save the points in a comma seperated value file
  # csv points can be imported into GIS for further processing
  write.csv(temp.points, paste("data/sim.dance_", tempData$dance.id, ".csv", sep = ""), row.names = FALSE)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    write.csv(temp.points, "data/simAllDances.csv", row.names = FALSE)
  }else{
    # save the total counts
    tempData2 <- read.csv("data/simAllDances.csv")
    tempData2 <- rbind(tempData2, temp.points)
    write.csv(tempData2, "data/simAllDances.csv", row.names = FALSE)
  }

  # georeference the points on the UK grid
  coordinates(temp.points) = c("easting", "northing")
  #proj4string(temp.points) = CRS("+init=epsg:27700")
  proj4string(coordPanel) = CRS("+init=epsg:26917") # Sponsler: our CRS is UTM 17N (EPSG:26917)

  # create a new raster with the temp.rast extent and sample the points on the raster / final sample size
  # e.g. probability that a dance has been on a certain raster square
  #temp.rast.UKGRID <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize
  temp.rast.UTM17N <- rasterize(temp.points, temp.rast, fun = "count", background = 0)/finalSampleSize # Sponsler:

  # convert the raster to a SpatialGridDataFrame
  #g <- as(temp.rast.UKGRID, 'SpatialGridDataFrame')
  g <- as(temp.rast.UTM17N, 'SpatialGridDataFrame') # Sponsler:

  # save the file to disk (can be imported into GIS)
  currentFileName <- paste("data/raster_", tempData$dance.id, ".asc", sep = "")
  write.asciigrid(g, currentFileName)

  if(i <= 1){ # on the first pass create a new file, else add the data to the existing file
    #total.temp.rast <- temp.rast.UKGRID
    total.temp.rast <- temp.rast.UTM17N # Sponsler:
  }else{
    # calculate the probability that a field has been visited
    #total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UKGRID)
    total.temp.rast <- 1 - (1 - total.temp.rast)*(1 - temp.rast.UTM17N)
  }
}

# save the combined dances as one raster file ready to be imported in ArcGIS
total.temp.rast <- total.temp.rast$dance.id
g.total <- as(total.temp.rast, 'SpatialGridDataFrame')
write.asciigrid(g.total, "data/totalRaster.asc")


# for plotting, we can crop the extent to our needs

# prepare an georeferenced extent (in GIS terms) for the cropping
coordPanel.crop <- data.frame(cbind(easting = c(hiveEasting, hiveEasting - 10000,
                                 hiveEasting + 10000),
                               northing = c(hiveNorthing, hiveNorthing - 10000,
                                 hiveNorthing + 10000)))
coordinates(coordPanel.crop) <- c("easting", "northing")
#proj4string(coordPanel.crop) = CRS("+init=epsg:27700")
proj4string(coordPanel.crop) = CRS("+init=epsg:26917") # Sponsler:


# prepare a raster with the extent to crop the aerial photography
numCoordPanel.crop <- as.data.frame(coordPanel.crop)
crop.rast <- raster(ncols = noCells, nrows = noCells)
extent(crop.rast) <- extent(c(numCoordPanel.crop[2:3,1], numCoordPanel.crop[2:3,2]))
#proj4string(crop.rast) = CRS("+init=epsg:27700")
proj4string(crop.rast) = CRS("+init=epsg:26917") # Sponsler:

# we crop the data raster to size
new.data.rast <- crop(total.temp.rast, crop.rast)
writeRaster(new.data.rast, filename = "fsr_may24.tif", format = "GTiff", overwrite = T) # Sponsler: this geotiff can be loaded in QGIS to overlay on landscape layer


### PLOTTING WITHOUT AERIAL PHOTOGRAPHY
# choose colours of your liking, and symbol size that matched the figure
myAlpha <- c(0, seq(0.005, 0.5, 0.005) + 0.3)
myCols <- rev(rainbow(100, alpha = rev(myAlpha)))
charEx <- 2

# Sponsler: plot vector landscape layer -- this doesn't quite work yet; alignment of the two plots is off
#mooreman <- readOGR("/Users/dougsponsler/Documents/Research/CDRC_dance_analysis", "2015_sites_1500m_buffer_squares_site__MO_clip")
#plot(mooreman, col = c("yellow", "green", "gray"))
#par(new = T)

# plot the data raster
plot(new.data.rast, col = myCols, legend.width = charEx/2,
     smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)


# Sponsler: still need to work out the aerial photography part, though we may end up plotting with the vector layers instead

### PLOTTING WITH AERIAL PHOTOGRAPHY
# if you have access to aerial photography (i.e. a GeoTIFF), you can load
# that now and crop it to size
# new.raster <- brick("someAerialPhotography.tif")
# new.raster <- crop(new.raster, crop.rast)

# plot the aerial photography
# plotRGB(new.raster, r=1, g=2, b=3)

# plot a background for the colour scale
# polygon(c(hiveEasting + 50, hiveEasting + 50, hiveEasting + 250, hiveEasting + 250, hiveEasting + 50),
    #    c(hiveNorthing + 50, hiveNorthing + 1075, hiveNorthing + 1075, hiveNorthing + 50, hiveNorthing + 50), col = "white", border = "black", lwd = 2)

# plot the data raster
# plot(new.data.rast, add = TRUE, col = myCols, legend.width = charEx/2,
   #  smallplot=c(.875,.9,0.33,0.9), axis.args = list(cex.axis = charEx/2))

# plot the hive location
#points(hiveEasting, hiveNorthing, pch = 17, col = "white", cex = 2*charEx)
# points(hiveEasting, hiveNorthing, pch = 20, col = "black", cex = charEx/2) # Sponsler:
#points(hiveEasting, hiveNorthing, pch = 17, col = "red", cex = charEx)

# plot the feeder location
#points(533999, 109254, pch = 19, col = "white", cex = 2*charEx)
#points(533999, 109254, pch = 19, col = "black", cex = charEx)



# get rid of the simulation data in memory
rm(sim.distances)
rm(sim.heading)
rm(rel.dance.easting)
rm(rel.dance.northing)
rm(samples)
rm(temp.points)
rm(g.total)
rm(total.temp.rast)
rm(temp.rast.UKGRID)

rm(list = ls())
