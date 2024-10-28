
# Read in the query file from 2007 - 2017 
NDVIfiles <- read.csv("Spatial_Layers/LAADS_query.2018-01-09T16_59.csv")

# Make a vector of file names
FilesToGet <- paste0("https://ladsweb.modaps.eosdis.nasa.gov",NDVIfiles[,2])

# Save file names
FileNames <- substr(FilesToGet,75,119)

# Destination files
destFileNames <- paste0("Spatial_Layers/NDVI/hdf/",FileNames)

# Download the NDVI files
mapply(function(x,y) download.file(x,y,mode = 'wb', extra = getOption("download.file.extra")),FilesToGet,destFileNames)

# Where to save tiff
GeoFiles <- substr(FileNames,1,16)

library(gdalUtils)

# Get the hdf files 
hdffiles <- list.files("Spatial_Layers/NDVI/hdf",pattern = "*.hdf",full.names = TRUE)
geotiff <- paste0("Spatial_Layers/NDVI/geotiff/",GeoFiles,".tiff")

# Convert to geoTiff
mapply(function(x,y){ 
       b <- get_subdatasets(x)
       gdal_translate(b[1], dst_dataset = y)}, x = hdffiles, y = geotiff)

#
NDVI<- sapply(list.files("Spatial_Layers/NDVI/geotiff", pattern = ".tiff", full.names = TRUE),raster::raster)
stack(NDVI)

JAM <- raster::getData("GADM",country = "Jamaica", level = 0)
JAM <- sp::spTransform(JAM, sp::CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))


library(raster)
library(sp)

# Read in shapefiles of the plots #

plotfiles <- list.files("Spatial_Layers/Jamaica Plots",pattern = "Point_ge.shp",recursive = TRUE, full.names = TRUE)

expplots <- list.files("Spatial_Layers/Jamaica Plots/Jamaica Plots", pattern = "Exp_Plots.shp", full.names = TRUE)

# use the raster package to read in the shapefiles #
JAMplots <- lapply(c(plotfiles,expplots), shapefile)

# Name the shapefiles #
names(JAMplots) <- c("Logwood_II","Logwood_I","MainMangrove","Mangrove_I","Mangrove_II","Experimental")

# Put the name used in the Banding data into the shapefiles #
JAMplots[[1]]$PLOT <- "LII"
JAMplots[[2]]$PLOT <- "LI"
JAMplots[[3]]$PLOT <- "MM"
JAMplots[[4]]$PLOT <- "MI"
JAMplots[[5]]$PLOT <- "MII"

# For the first 5 plots - create variables and reorder so we can stack them all into a single file #
JAMplots[c(1:5)] <- lapply(JAMplots[c(1:5)], FUN = function(x){
                          x$LONG <- x@coords[,1]
                          x$LAT <- x@coords[,2]
                          x@coords <- x@coords[,1:2]
                          y <- x[,c("PLOT","Comment","LONG","LAT")]
                          names(y) <- c("PLOT","LOCATION","LONG","LAT")
                          return(y)})

# Change the projection from WGS84 to Equidistant #
JAMplots <- lapply(JAMplots, spTransform, CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"))

# Merge all the plots into a single shapfile #
JAMplots <- do.call("rbind", JAMplots) 

JAMplots$LONG_UTM <- JAMplots@coords[,1]
JAMplots$LAT_UTM <- JAMplots@coords[,2]

JAMplots <- sp::spTransform(JAMplots, sp::CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"))

# Extract values to each point #

plotvals <- lapply(NDVI, FUN = function(x){y <- extract(x,JAMplots)
                                           yearday <- gsub(x = names(x),pattern = "MOD13Q1.A", replacement = "")
						       year <- substr(yearday,1,4)
                                           day <- as.numeric(substr(yearday,5,7))
                                           z <- data.frame(year = year,
                                                           day = day,
                                                           date = as.Date(day,paste0(year,"-01-01")),
                                                           plot = JAMplots$PLOT,
                                                           ndvi = y)
                                           return(z)})

NDVIplots <- do.call("rbind",plotvals)

NDVIplots$month <- format(NDVIplots$date,"%m")
rownames(NDVIplots)<-c()


library(dplyr)

NDVImean <- data.frame(group_by(NDVIplots,year,month) %>%
               summarize(mean.ndvi = mean(ndvi/100000000,na.rm = TRUE)))

#write.csv(NDVImean,"Data/NDVImonthly2007_2017.csv",row.names = FALSE)



####################### PRECIP - TRMM - TROPICAL RAINFALL #############################
library(raster)
# Get the hdf files 
hdffiles <- list.files("E:/TRMM",pattern = "*.HDF",full.names = TRUE)
geotiff <- gsub(x=hdffiles,pattern = "HDF",replacement = "tiff")

# Convert to geoTiff
mapply(function(x,y){ 
       b <- gdalUtils::get_subdatasets(x)
       gdalUtils::gdal_translate(b[1], dst_dataset = y)}, x = hdffiles, y = geotiff)

#
geotiff <- list.files("E:/TRMM",pattern = "*.tiff", full.names = TRUE)
TRMM <- lapply(geotiff,raster::raster)
TRMM <- raster::stack(TRMM)

newExtent <- extent(-50,50,-180,180)

TRMM1 <- raster::setExtent(TRMM,ext=newExtent)
TRMM2 <- t(TRMM1)

plot(TRMM2[[1]])
crs(TRMM2)<-sp::CRS("+init=epsg:4326")

JAMplots <- shapefile("Spatial_Layers/JAMplots_all.shp")
JAMplots <- sp::spTransform(JAMplots,sp::CRS("+init=epsg:4326"))

JAM <- getData("GADM",country = "Jamaica",level = 1)

str(rain03)
rain01 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0101.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain02 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0201.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain03 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0301.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain04 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0401.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain05 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0501.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain06 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0601.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain07 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0701.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain08 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0801.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain09 <- extract(TRMM2[[grep(x=hdffiles,pattern = "0901.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain10 <- extract(TRMM2[[grep(x=hdffiles,pattern = "1001.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain11 <- extract(TRMM2[[grep(x=hdffiles,pattern = "1101.7")]],JAM[JAM$NAME_1=="Westmoreland",])
rain12 <- extract(TRMM2[[grep(x=hdffiles,pattern = "1201.7")]],JAM[JAM$NAME_1=="Westmoreland",])

Year_rain <- rain01[[1]][,1:20]+rain02[[1]]+rain03[[1]]+rain04[[1]]+rain05[[1]]+rain06[[1]]+
             rain07[[1]]+rain09[[1]]+rain10[[1]]+rain11[[1]]+rain12[[1]]

polygon(x = c(10,10,20,20),y = c(0,0.5,0.5,0))

plot(apply(Year_rain,2,mean),type = "o", xaxt = "n",las = 2, ylab = expression(mm^-h), pch = 19, cex = 1.5)
axis(1, at = 1:20,labels = 1998:2017)

plot(TRMM2[[1]],ext = extent(JAMplots)+5)
url <- "ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/TRMM_3B43M/TRMM_3B43M_"
YrsToGet <- 1998:2018
MonthsToGet <- formatC(1:12,width = 2, flag = "0")
years <- rep(YrsToGet,each = length(MonthsToGet))
months <- rep(MonthsToGet,length(YrsToGet))

FilesToGet <- paste0(url,years,"-",months,".FLOAT.TIFF")
DestFiles <- paste0("E:/TRMM/TRMM_",years,months,"_neo.tiff")

# Download the NDVI files
mapply(function(x,y) download.file(x,y,mode = 'wb', extra = getOption("download.file.extra")),FilesToGet,DestFiles)

TRMM <- lapply(list.files("Spatial_Layers/TRMM/",pattern = "neo.tiff", recursive = TRUE, full.names = TRUE),raster)

plotvals <- lapply(TRMM, FUN = function(x){y <- extract(x,JAMplots)
                                           yearday <- gsub(x = names(x),pattern = "TRMM_", replacement = "")
						       year <- substr(yearday,1,4)
                                           month <- as.numeric(substr(yearday,5,6))
                                           z <- data.frame(year = year,
                                                           month = month,
                                                           plot = JAMplots$PLOT,
                                                           monthlyPrecip = y)
                                           return(z)})

precip <- do.call("rbind",plotvals)

library(dplyr)
precipMonth <- data.frame(group_by(precip,year,month) %>%
               summarize(MonthlyPrecip = mean(monthlyPrecip,na.rm = TRUE)))

#write.csv(precipMonth,"Data/PrecipMonthly2007_2017.csv",row.names = FALSE)



# Download surface temperature data from NEO - Get the ones that are currenlty available

mon <- c("01","02","03","04","05","06","07","08","09","10","11","12")
yr <- 2000:2018

mons <- rep(mon,length(yr))
years <- rep(yr, each = length(mon))

temp2get <- paste0("ftp://neoftp.sci.gsfc.nasa.gov/geotiff.float/MOD11C1_M_LSTDA/MOD11C1_M_LSTDA","_",years,"-",mons,".FLOAT.TIFF")
temp2get <- temp2get[2:206]

outtemp <- paste0("E:/MOD11C1_M_LSTDA/MOD11C1_M_LSTDA_",years[2:206],"_",mons[2:206],".FLOAT.TIFF")


mapply(x = temp2get[141:205], y = outtemp[141:205], FUN = function(x,y){
          download.file(url = x, destfile = y)})

               summarize(MonthlyPrecip = mean(monthlyPrecip,na.rm = TRUE)))
