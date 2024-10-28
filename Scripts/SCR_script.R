#############################################################################
# Organize data to create a capture history for each individual for each
# mist-net location within the sampling grid. These data are required
# to estimate population size and density. Once these data are collated 
# we can extend the analyses to include 'open-population' models that account
# for recruitment between years (when the population is open). 
#############################################################################

# First, we'll try this approach with American Redstart's then apply it 
# across other taxa once the format / structure is correct. 

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

# shapefile(JAMplots,"Spatial_Layers/JAMplots_all.shp")

plot(JAMplots, pch = 20, col = as.numeric(as.factor(JAMplots$PLOT)))

plo <- data.frame(point = as.character(paste(JAMplots$PLOT,gsub(x = JAMplots$LOCATION,pattern = " ", replacement = "-"))),
           x  = JAMplots@coords[,1],
           y = JAMplots@coords[,2])

# Read in the plots #
JAMplots <- shapefile("Spatial_Layers/JAMplots_all.shp")

# Read in the spatial locations of the net lanes #
LTREBnets <- raster::shapefile("Spatial_Layers/LTREBnetlanes.shp")

# Create the general location of each net set by creating a dummy variable #
LTREBnets$NetLoc <- paste0(LTREBnets$Plot,"-",substr(LTREBnets$Location,1,3))

# Convert it to numeric # 
LTREBnets$NetLocNum <- as.numeric(as.factor(LTREBnets$NetLoc))

# Create a Plot Location Net ID that matches the banding data #
LTREBnets$PlotLoc <- paste0(LTREBnets$Plot," ",gsub(" ","-",LTREBnets$Location))

# Find the number of net locations #
netlocs <- max(LTREBnets$NetLocNum)

# Generate a polygon for each net location that represents the area surveyed #
EffortPolys <- vector("list",netlocs)

for(n in seq(netlocs)){
poly <- as(raster::extent(LTREBnets[LTREBnets$NetLocNum == n,]), "SpatialPolygons")
polybuff <- rgeos::gBuffer(poly,width = 5,byid = TRUE)
EffortPolys[[n]] <- polybuff
}

LTREBnetLocs <- paste0(LTREBnets$Plot," ",gsub(x=LTREBnets$Location,pattern = " ", replacement = "-"))

# crs(LTREBnets) <- "+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"

# Randomly sample points along a line #
sp::spsample(LTREBnets[1,],n = 10, type = "random")


# Read in data from Github#
LTREB <- read.csv("https://raw.githubusercontent.com/Dossman/LTREB-Survival-Analysis/master/data/complete%20banding.csv")

# Subset to only include those captures by the LTREB team #
LTREB <- subset(LTREB, BNDCREW == "LTREB")

# Take a look at the plots that were used through time #
table(LTREB$YEAR,LTREB$CAPSITE)

# Covert dates and get month of capture #
LTREB$CAPDATE <- as.POSIXct(LTREB$CAPDATE,format = "%d-%b-%y")
LTREB$CAPMON <- as.numeric(format(LTREB$CAPDATE, "%m"))

# Combine capture site and capture location # 
LTREB$CAPLOC <- paste(LTREB$CAPSITE,LTREB$LOCATION)

# How much of the data are ready to use without fixing capture location? #
sum(LTREB$CAPLOC %in% LTREBnetLocs)/length(LTREB$CAPLOC)

# Find capture locations that aren't in the standard format ("LII K-10-N")
CapLocs_to_change <- unique(LTREB$CAPLOC[!(LTREB$CAPLOC %in% LTREBnetLocs)])

CapLocs_to_change[grep(x = CapLocs_to_change, pattern = "LI ")]

# Change unrecognized locations into standard format 
################################ LOGWOOD I ########################################
par(mar = c(0,0,0,0))
plot(JAMplots[JAMplots$PLOT == "LI",],pch = ".")
text(JAMplots@coords, JAMplots$LOCATION, cex = 0.6)
plot(LTREBnets,add = TRUE)
#text(maptools::SpatialLinesMidPoints(LTREBnets),LTREBnets$PlotLoc,cex = 0.5)

# LI A-5-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A-4[.][0-9]$")]<-"LI A-5-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI Z[.][4-9]-5$")]<-"LI A-5-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A-5[.][0-9]$")]<-"LI A-5-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A[.][0-9]-[5]$")]<-"LI A-5-W"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A[.][0-9]-[4-6][.][0-9]$")]<-"LI A-5-W"

# LI B-3-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B-2[.][0-9]$")]<-"LI B-3-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A[.][4-9]-3$")]<-"LI B-3-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A[.][4-9]-[2-3][.][0-9]$")]<-"LI B-3-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B-3[.][0-9]$")]<-"LI B-3-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][0-9]-[3]$")]<-"LI B-3-W"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][0-9]-[2-4][.][0-9]$")]<-"LI B-3-W"

# LI C-5-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C-4[.][0-9]$")]<-"LI C-5-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][4-9]-5$")]<-"LI C-5-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][4-9]-[4-5][.][0-9]")]<-"LI C-5-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C-5[.][0-9]$")]<-"LI C-5-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C[.][0-9]-[5]$")]<-"LI C-5-W"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C[.][0-9]-[4-6][.][0-9]$")]<-"LI C-5-W"

# LI C-7-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C-6[.][0-9]$")]<-"LI C-7-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][4-9]-7$")]<-"LI C-7-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][4-9]-[6-8][.][0-9]$")]<-"LI C-7-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C-7[.][0-9]$")]<-"LI C-7-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C[.][0-9]-[7]$")]<-"LI C-7-W"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI C[.][0-9]-[6-8][.][0-9]$")]<-"LI C-7-W"

# LI B-8-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B-7[.][0-9]$")]<-"LI B-8-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A[.][4-9]-8$")]<-"LI B-8-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI A[.][4-9]-[7-9][.][0-9]$")]<-"LI B-8-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B-8[.][0-9]$")]<-"LI B-8-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][0-9]-[8]$")]<-"LI B-8-W"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LI B[.][0-9]-[7-9][.][0-9]$")]<-"LI B-8-W"


# Find capture locations that aren't in the standard format ("LII K-10-N")
CapLocs_to_change <- unique(LTREB$CAPLOC[!(LTREB$CAPLOC %in% LTREBnetLocs)])

CapLocs_to_change[grep(x = CapLocs_to_change, pattern = "LI ")]

# Fixing the stragglers ##
LTREB$CAPLOC <- gsub(x=LTREB$CAPLOC,pattern = "LI B-9", replacement = "LI B-8")

################################ LOGWOOD II ########################################
par(mar = c(0,0,0,0))
plot(JAMplots[JAMplots$PLOT == "LII",],pch = ".")
text(JAMplots@coords, JAMplots$LOCATION, cex = 0.6)
plot(LTREBnets,add = TRUE)
#text(maptools::SpatialLinesMidPoints(LTREBnets),LTREBnets$PlotLoc,cex = 0.5)

# LII K-6-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII K[.][0-9]-[6]$")]<-"LII K-6-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII K-6[.][0-9]$")]<-"LII K-6-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII J[.][4-9]-6$")]<-"LII K-6-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII K-5[.][0-9]$")]<-"LII K-6-W"

# LII K-10-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII K[.][0-9]-[10]$")]<-"LII K-10-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII K-10[.][0-9]$")]<-"LII K-10-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII J[.][4-9]-10$")]<-"LII K-10-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII K-9[.][0-9]$")]<-"LII K-10-W"

# LII J-8-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII J[.][0-9]-[8]$")]<-"LII J-8-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII J-8[.][0-9]$")]<-"LII J-8-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I[.][4-9]-8$")]<-"LII J-8-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII J-7[.][0-9]$")]<-"LII J-8-W"

# LII I-10-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I[.][0-9]-[10]$")]<-"LII I-10-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I-10[.][0-9]$")]<-"LII I-10-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII H[.][4-9]-10$")]<-"LII I-10-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I-9[.][0-9]$")]<-"LII I-10-W"

# LII I-6-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I[.][0-9]-[6]$")]<-"LII I-6-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I-6[.][0-9]$")]<-"LII I-6-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII H[.][4-9]-6$")]<-"LII I-6-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^LII I-5[.][0-9]$")]<-"LII I-6-W"


# Find capture locations that aren't in the standard format ("LII K-10-N")
CapLocs_to_change <- unique(LTREB$CAPLOC[!(LTREB$CAPLOC %in% LTREBnetLocs)])

CapLocs_to_change[grep(x = CapLocs_to_change, pattern = "LII ")]

################################ MANGROVE I ########################################
par(mar = c(0,0,0,0))
plot(JAMplots[JAMplots$PLOT == "MI",],pch = ".")
text(JAMplots@coords, JAMplots$LOCATION, cex = 0.6)
plot(LTREBnets,add = TRUE)
#text(maptools::SpatialLinesMidPoints(LTREBnets),LTREBnets$PlotLoc,cex = 0.5)

# MI B-2-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI B[.][0-9]-[2]$")]<-"MI B-2-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI B-2[.][0-9]$")]<-"MI B-2-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI A[.][4-9]-2$")]<-"MI B-2-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI B-1[.][0-9]$")]<-"MI B-2-W"

# MI D-3-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI D[.][0-5]-[3]$")]<-"MI D-3-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI D-3[.][0-9]$")]<-"MI D-3-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI C[.][4-9]-3$")]<-"MI D-3-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI D-2[.][0-9]$")]<-"MI D-3-W"

# MI E-3-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI E[.][0-9]-[3]$")]<-"MI E-3-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI E[.][0-9]-[2-3][.][0-9]$")]<-"MI E-3-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI E-3[.][0-9]$")]<-"MI E-3-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI D[.][6-9]-3$")]<-"MI E-3-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI E-2[.][0-9]$")]<-"MI E-3-W"

# MI G-3-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI G[.][0-9]-[3]$")]<-"MI G-3-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI G[.][0-9]-[3][.][0-9]$")]<-"MI G-3-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI G-3[.][0-9]$")]<-"MI G-3-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI F[.][0-9]-3$")]<-"MI G-3-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI F[.][0-9]-[2-3][.][0-9]$")]<-"MI G-3-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI G-2[.][0-9]$")]<-"MI G-3-W"

# MI H-2.5-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI H[.][0-9]-[2.5]$")]<-"MI H-2.5-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI H[.][0-9]-[2-3][.][0-9]$")]<-"MI H-2.5-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI H-2[.][5-9]$")]<-"MI H-2.5-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI G[.][0-9]-2.5$")]<-"MI H-2.5-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI H-2[.][0-5]$")]<-"MI H-2.5-W"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MI H[.][0-9]-[1-2][.][0-9]$")]<-"MI H-2.5-W"

# Find capture locations that aren't in the standard format ("LII K-10-N")
CapLocs_to_change <- unique(LTREB$CAPLOC[!(LTREB$CAPLOC %in% LTREBnetLocs)])

CapLocs_to_change[grep(x = CapLocs_to_change, pattern = "MI ")]


################################ MANGROVE II ########################################
par(mar = c(0,0,0,0))
plot(JAMplots[JAMplots$PLOT == "MII",],pch = ".")
text(JAMplots@coords, JAMplots$LOCATION, cex = 0.6)
plot(LTREBnets,add = TRUE)
#text(maptools::SpatialLinesMidPoints(LTREBnets),LTREBnets$PlotLoc,cex = 0.5)

# MII F-5-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII F[.][0-9]-[5]$")]<-"MII F-5-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII F-5[.][0-9]$")]<-"MII F-5-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII E[.][4-9]-5$")]<-"MII F-5-S"
#LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII A-6[.][0-9]$")]<-"MII F-5-W"

# MII F-9-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII F[.][0-9]-[9]$")]<-"MII F-9-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII F-9[.][0-9]$")]<-"MII F-9-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII E[.][4-9]-9$")]<-"MII F-9-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII F-8[.][0-9]$")]<-"MII F-9-W"

# MII E-8-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII E[.][0-9]-[8]$")]<-"MII E-8-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII E-8[.][0-9]$")]<-"MII E-8-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII D[.][4-9]-8$")]<-"MII E-8-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII E-7[.][0-9]$")]<-"MII E-8-W"

# MII D-9-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII D[.][0-9]-[9]$")]<-"MII D-9-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII D-9[.][0-9]$")]<-"MII D-9-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII C[.][4-9]-9$")]<-"MII D-9-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MII D-8[.][0-9]$")]<-"MII D-9-W"

# Find capture locations that aren't in the standard format ("LII K-10-N")
CapLocs_to_change <- unique(LTREB$CAPLOC[!(LTREB$CAPLOC %in% LTREBnetLocs)])

CapLocs_to_change[grep(x = CapLocs_to_change, pattern = "MII ")]

################################ MAIN MANGROVE  ########################################
par(mar = c(0,0,0,0))
plot(JAMplots[JAMplots$PLOT == "MM",],pch = ".")
text(JAMplots@coords, JAMplots$LOCATION, cex = 0.6)
plot(LTREBnets,add = TRUE)
#text(maptools::SpatialLinesMidPoints(LTREBnets),LTREBnets$PlotLoc,cex = 0.5)

# MM B-4-*
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MM B-3[.][0-9]$")]<-"MM B-4-N"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MM A[.][0-9]-4$")]<-"MM B-4-E"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MM B-4[.][0-9]$")]<-"MM B-4-S"
LTREB$CAPLOC[grep(x = LTREB$CAPLOC, pattern = "^MM B[.][0-9]-[4]$")]<-"MM B-4-W"

# Change all MII B-4 to MM
LTREB$CAPLOC <- gsub(x = LTREB$CAPLOC, pattern = "MII B-4", replacement = "MM B-4")

# Find capture locations that aren't in the standard format ("LII K-10-N")
CapLocs_to_change <- unique(LTREB$CAPLOC[!(LTREB$CAPLOC %in% LTREBnetLocs)])

CapLocs_to_change[grep(x = CapLocs_to_change, pattern = "MM ")]

######################################################################################
#
# Generate a random point along the netlane a bird was caught in. 
# 
######################################################################################

# Randomly sample points along the net an individual was captured in #
# At most this adds 12m of uncertainty to movements                  #

LTREB$LONG <- NA
LTREB$LAT <- NA
for(i in 1:nrow(LTREB)){
if(LTREB[i,grep(names(LTREB),pattern = "CAPLOC")] %in% LTREBnets$PlotLoc){
net <- which(LTREBnets$PlotLoc == LTREB[i,grep(names(LTREB),pattern = "CAPLOC")])
pt<-sp::spsample(LTREBnets[net,],n = 1, type = "random")
LTREB$LONG[i]<-pt@coords[1,1]
LTREB$LAT[i]<-pt@coords[1,2]
}else{next}
}


# Take a look at how many captures for each species in the dataset #
table(LTREB$SPECIES)

# spp <- unique(LTREB$SPECIES)

# Keep only the species with more than 10 captures #
spp <- names(which(table(LTREB$SPECIES) > 10))
years <- unique(LTREB$YEAR)
years <- years[order(years)]
month <- as.numeric(unique(LTREB$CAPMON))
month <- month[order(month)]

# How many nets are there in LTREB - 80
n.nets <- length(LTREBnets)

netID <- LTREBnets$PlotLoc
# Get location for the center of the net #
netLocs <- maptools::SpatialLinesMidPoints(LTREBnets)@coords

###################################################################################
#
# Define the state-space (patchy habitat - Mike Meredith approach
#
###################################################################################

# identify plots as shapefiles #
MM <- JAMplots[JAMplots$PLOT == "MM",]
LI <- JAMplots[JAMplots$PLOT == "LI",]
MI <- JAMplots[JAMplots$PLOT == "MI",]
LII <- JAMplots[JAMplots$PLOT == "LII",]
MII <- JAMplots[JAMplots$PLOT == "MII",]

# create a single list of the plots #
plots <- list(MM,LI,MI,LII,MII)

# generate convex hull around each plot 
plotpolys <- lapply(plots,rgeos::gConvexHull)

plotspecificarea <- lapply(plotpolys,FUN = function(x){rgeos::gArea(x)/10000})

# 5 x 5 m raster of plot extent
JAMraster <- raster(extent(JAMplots),res=c(5,5))

# 5 x 5 m raster of plot extent
#JAMraster <- raster(extent(JAMplots),res=c(10,10))

# set values to 0 - i.e., not sampled 
JAMraster[]<-0

# make rasters for each plot and give them values 1 - 5
plotrasters <- vector("list",5)
for(i in 1:5){
plotrasters[[i]] <- rasterize(plotpolys[[i]],JAMraster)
plotrasters[[i]][plotrasters[[i]]==1]<-i
}

# Calculate area of the 5 plots in hectares
(area(plotpolys[[1]])+area(plotpolys[[2]])+area(plotpolys[[3]])+area(plotpolys[[4]])+area(plotpolys[[5]]))/10000

# habmatrix 
study.area <- sum(stack(c(JAMraster,plotrasters)),na.rm = TRUE)

# create matrix
habmat <- as.matrix(study.area)

L <- rasterToPoints(study.area)

centerNetLocs <- rowColFromCell(JAMraster,extract(JAMraster,SpatialPoints(cbind(netLocs[,1],netLocs[,2])),cellnumbers=TRUE)[,1])

rowcol <- rowColFromCell(study.area,extract(study.area,SpatialPoints(cbind(L[,1],L[,2])),cellnumbers=TRUE)[,1])

# Confirm that the plot raster and centerNetLocs don't match up 
plot(raster(habmat,xmn = 0,xmx = ncol(habmat),ymn=0,ymx=nrow(habmat)))
points(centerNetLocs ,pch = 19)

# Flip the plot raster to match the net locations 
study.flip <- flip(raster(habmat,xmn=0, xmx=ncol(habmat), ymn=0, ymx=nrow(habmat)),"y")
habmat <- as.matrix(study.flip)

# Confirm that flipping the plots in y direction fixed the problem 
plot(study.flip)
points(centerNetLocs[,1]~centerNetLocs[,2])

ylimit <- nrow(habmat)
xlimit <- ncol(habmat)

#habmat[habmat>1]<-1

# Yearly temp, precip and SOI 

MoBay <- read.csv("Data/MontegoBay_WeatherStation.csv")

MoBay$date <- as.POSIXct(MoBay$DATE, format = "%Y-%m-%d")

MoBay$Year <- format(MoBay$date,"%Y")
MoBay$Month <- format(MoBay$date, "%m")


MonthTemp <- tapply(MoBay$TAVG, list(MoBay$Year,MoBay$Month), FUN = mean, na.rm = TRUE)
MonthPrecip <- tapply(MoBay$PRCP, list(MoBay$Year,MoBay$Month), FUN = mean, na.rm = TRUE)

# SOI <- as.matrix(read.table("https://crudata.uea.ac.uk/cru/data/soi/soi.dat",header = FALSE))
SOI <- read.csv("Data/SOIindex.csv")

colnames(SOI) <- c("Year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Annual")

SOI[SOI==-99.99]<-NA

SOI <- SOI[SOI[,1] > 2007,]

# Years 2008-2017 mean Jan-Mar
T<-apply(MonthTemp[11:20,1:3],1,mean)
P<-apply(MonthPrecip[11:20,1:3],1,sum,na.rm = TRUE)
S <- SOI[1:10,14] 

StandTemp <- (T-mean(T))/sd(T)
StandPrecip <- (P-mean(P))/sd(P)
SOIindex <- S
SOIindex_poly <- poly(S,2)
######################################################################################################################
#
# Single Species - save win.data & inits for each species separately
#
######################################################################################################################

# How many nets are there in LTREB - 80
n.nets <- length(LTREBnets)

netID <- LTREBnets$PlotLoc

# Get location for the center of the net #
netLocs <- maptools::SpatialLinesMidPoints(LTREBnets)@coords

# Get species that have more than 10 captures
spp <- names(which(table(LTREB$SPECIES)>10))

# How many species is that?
nspp <- length(spp)

ncaptured <- tapply(LTREB$SPECIES[LTREB$SPECIES %in% spp],list(LTREB$SPECIES[LTREB$SPECIES %in% spp]),FUN=length)
ncaptured <- ncaptured[!is.na(ncaptured)]

M <- ifelse(ncaptured>400,800,500)


# Bird x month x year - but instead of 1 and 0 it the matrix should contain the 
# net number it was captured in. If not captured (CH = 0) then it should get net
# 81 - (n.nets + 1). 
# Create an empty capture history filling with 81 #
for(s in 1:nspp){

MS_scr_CH <- array(81, dim = c(M[s],length(month),length(years)))

ltreb <- LTREB[LTREB$SPECIES == spp[s],]

ninds <- length(unique(ltreb$BAND))
indid <- unique(ltreb$BAND)

# Find where in CH to place 1 for capture
rowID <- pmatch(ltreb$BAND,indid,duplicates.ok = TRUE)
monID <- pmatch(as.numeric(ltreb$CAPMON),as.numeric(month),duplicates.ok = TRUE)
yearID <- pmatch(ltreb$YEAR,years,duplicates.ok = TRUE)
locID <- pmatch(ltreb$CAPLOC,netID, duplicates.ok = TRUE)

SEX <- ltreb$SEX
sex <- ifelse(SEX=="M",1,2)

for(i in 1:length(rowID)){
MS_scr_CH[rowID[i],monID[i],yearID[i]]<-locID[i]
}


# For now some birds have NA for capture location - 
# those need to be fixed but I gave them net = 81 (i.e, not captured)
MS_scr_CH[is.na(MS_scr_CH)]<-81

win.data <- list(nnets = n.nets,
                 Ycat = MS_scr_CH,
                 nmonths = length(month),
                 nyears = length(years),
                 xlimit = xlimit,
                 ylimit = ylimit,
                 NetLoc = centerNetLocs,
                 habmat = habmat,
                 M = M[s],
                 sex = sex,
                 precip = StandPrecip,
                 temp = StandTemp,
                 soi = SOIindex,
                 soi2 = SOIindex_poly[,2])

# initial values for survival #
z_start <- array(NA,c(win.data$M,win.data$nyears))

z_start[,1:10] <- c(rep(1,ninds),rep(0,win.data$M-ninds)) #BTBW


# initial values for activity center 
Sst <- array(NA,c(win.data$M,2,length(years)))

# Sst[,1:2,] <- cbind(runif(M, 0, win.data$xlimit), 
#                     runif(M, 0, win.data$ylimit))

in_hab <- which(habmat != 0, arr.ind = TRUE)

Sst[,1:2,] <- in_hab[sample(x = 1:nrow(in_hab),size = win.data$M,  replace = FALSE),]


saveRDS(win.data,paste0("Data/",spp[s],"_windata_10m.rds"))
saveRDS(z_start,paste0("Data/",spp[s],"_z_start_10m.rds"))
saveRDS(Sst,paste0("Data/",spp[s],"_Sst_10m.rds"))
}

inits <- function(){list(z = z_start,
                  S = Sst)}

params <- c("psi","sp.phi","backphi","intercept", "beta", "sigma", "N", "D","R","z","S")



AMRE_ltreb <- LTREB[LTREB$SPECIES == "AMRE",]

# Find where in CH to place 1 for capture
rowID <- pmatch(AMRE_ltreb$BAND,AMREid,duplicates.ok = TRUE)
monID <- pmatch(as.numeric(AMRE_ltreb$CAPMON),as.numeric(month),duplicates.ok = TRUE)
yearID <- pmatch(AMRE_ltreb$YEAR,years,duplicates.ok = TRUE)
locID <- pmatch(AMRE_ltreb$CAPLOC,netID, duplicates.ok = TRUE)


# Bird x month x year - but instead of 1 and 0 it the matrix should contain the 
# net number it was captured in. If not captured (CH = 0) then it should get net
# 81 - (n.nets + 1). 
# Create an empty capture history filling with 81 #

M <- 1000

AMRE_scr_CH <- array(81, dim = c(M,length(month),length(years)))
# give names to the capture history 
dimnames(AMRE_scr_CH)[[1]] <- c(as.character(AMREid),as.character(paste0("M",seq(1,M-nAMRE,1))))
dimnames(AMRE_scr_CH)[[2]] <- month
dimnames(AMRE_scr_CH)[[3]] <- years

for(i in 1:length(rowID)){
AMRE_scr_CH[rowID[i],monID[i],yearID[i]]<-locID[i]
}

# For now some birds have NA for capture location - 
# those need to be fixed but I gave them net = 81 (i.e, not captured)
AMRE_scr_CH[is.na(AMRE_scr_CH)]<-81


win.data <- list(nnets = n.nets,
                 Ycat = AMRE_scr_CH,
                 nmonths = length(month),
                 nyears = length(years),
                 xlimit = xlimit,
                 ylimit = ylimit,
                 NetLoc = centerNetLocs,
                 habmat = habmat,
                 M = dim(AMRE_scr_CH)[1])

# initial values for survival #
z_start <- array(NA,c(win.data$M,win.data$nyears)))
z_start[,1:10] <- c(rep(1,nAMRE),rep(0,M-nAMRE))

# initial values for activity center 
Sst <- array(NA,c(M,2,length(years)))

# Sst[,1:2,] <- cbind(runif(M, 0, win.data$xlimit), 
#                     runif(M, 0, win.data$ylimit))

in_hab <- which(habmat != 0, arr.ind = TRUE)

Sst[,1:2,] <- in_hab[sample(x = 1:10021,size = win.data$M,  replace = FALSE),]

inits <- function(){list(z = z_start,
                  S = Sst)}

params <- c("psi","phi","intercept", "beta", "sigma", "N", "D","R","pOK")

Sys.time()
a<-Sys.time()
M1 <- jagsUI::jags.basic(model = "C:/Users/hallworthm/Dropbox (Smithsonian)/LTREB_Jamaica/Models/multiseason_scr.txt",
             data = win.data,
             inits = inits,
             parameters.to.save = params,
             n.chains = 3,
             n.iter = 10,
             n.burn = 1,
             n.thin = 1,
             n.cores = 3,
             parallel = TRUE,
             seed = 12345)
Sys.time()-a
Sys.time()

source("C:/Users/hallworthm/Dropbox (Smithsonian)/Cost_of_Reproduction/Functions/sims.list.R")

M2 <- process.output(M1)
str(M2$sims.list)
plot(M2$mean$D)

plot(study.flip)
points(centerNetLocs[,1]~centerNetLocs[,2],pch = 19, col = "red", cex =1.25)
points(M2$mean$S[,1,10]~M2$mean$S[,2,10])

plot(M2$mean$N)



######################################################################################################################
#
# Multiple Species 
#
######################################################################################################################

# How many nets are there in LTREB - 80
n.nets <- length(LTREBnets)

netID <- LTREBnets$PlotLoc

# Get location for the center of the net #
netLocs <- maptools::SpatialLinesMidPoints(LTREBnets)@coords

# Get species that have more than 10 captures
spp <- names(which(table(LTREB$SPECIES)>10))

keepers <- c("AMRE","BANA","BAWW","BTBW","BWVI",
             "GABU","JAEU","JAOR","JWEV","LOKI",
             "NOPA","NOWA","OVEN","PRAW","SAFL",
             "SWWA","WCTH","WEWA","YEWA","YSGR")

spp <- spp[spp %in% keepers]

# How many species is that?
nspp <- length(spp)

# define the migrants
migs <- c("AMRE","BAWW","NOWA","OVEN","PRAW","PROW","SWWA")

# define migrant vs resident
mig.res <- ifelse(spp %in% migs, "Migrant","Resident")

M <- 1000

# Bird x month x year - but instead of 1 and 0 it the matrix should contain the 
# net number it was captured in. If not captured (CH = 0) then it should get net
# 81 - (n.nets + 1). 
# Create an empty capture history filling with 81 #

MS_scr_CH <- array(81, dim = c(M,length(month),length(years),nspp))

ninds <- rep(NA,nspp)

for(s in 1:nspp){
ltreb <- LTREB[LTREB$SPECIES == spp[s],]

ninds[s] <- length(unique(ltreb$BAND))
indid <- unique(ltreb$BAND)

# Find where in CH to place 1 for capture
rowID <- pmatch(ltreb$BAND,indid,duplicates.ok = TRUE)
monID <- pmatch(as.numeric(ltreb$CAPMON),as.numeric(month),duplicates.ok = TRUE)
yearID <- pmatch(ltreb$YEAR,years,duplicates.ok = TRUE)
locID <- pmatch(ltreb$CAPLOC,netID, duplicates.ok = TRUE)

for(i in 1:length(rowID)){
MS_scr_CH[rowID[i],monID[i],yearID[i],s]<-locID[i]
}
}

# For now some birds have NA for capture location - 
# those need to be fixed but I gave them net = 81 (i.e, not captured)
MS_scr_CH[is.na(MS_scr_CH)]<-81

win.data <- list(nnets = n.nets,
                 Ycat = MS_scr_CH,
                 nmonths = length(month),
                 nyears = length(years),
                 xlimit = xlimit,
                 ylimit = ylimit,
                 NetLoc = centerNetLocs,
                 habmat = habmat,
                 M = M,
                 nspp = nspp,
                 migstatus = ifelse(mig.res == "Migrant",1,2))

# initial values for survival #
z_start <- array(NA,c(win.data$M,win.data$nyears),win.data$nspp))

for(s in 1:win.data$nspp){
z_start[,1:10,s] <- c(rep(1,ninds[s]),rep(0,win.data$M-ninds[s])) #BTBW
}

# initial values for activity center 
Sst <- array(NA,c(win.data$M,2,length(years),win.data$nspp))

# Sst[,1:2,] <- cbind(runif(M, 0, win.data$xlimit), 
#                     runif(M, 0, win.data$ylimit))

in_hab <- which(habmat != 0, arr.ind = TRUE)

for(s in 1:win.data$nspp){
Sst[,1:2,,s] <- in_hab[sample(x = 1:10021,size = win.data$M,  replace = FALSE),]
}

inits <- function(){list(z = z_start,
                  S = Sst)}

params <- c("psi","sp.phi","backphi","intercept", "beta", "sigma", "N", "D","R","z","S")


#saveRDS(win.data,"Data/win_data.rds")
#saveRDS(z_start,"Data/z_start.rds")
#saveRDS(Sst,"Data/Sst.rds")

Sys.time()
a<-Sys.time()
M1 <- jagsUI::jags.basic(model = "C:/Users/hallworthm/Dropbox (Smithsonian)/LTREB_Jamaica/Models/MS_multiseason_scr.txt",
             data = win.data,
             inits = inits,
             parameters.to.save = params,
             n.chains = 3,
             n.iter = 10,
             n.burn = 1,
             n.thin = 1,
             n.cores = 3,
             parallel = TRUE,
             seed = 12345)
Sys.time()-a
Sys.time()

source("C:/Users/hallworthm/Dropbox (Smithsonian)/Cost_of_Reproduction/Functions/sims.list.R")

M2 <- process.output(M1)

plot(M2$mean$D[1,])
points(1:10,M2$mean$D[2,],pch = 19)














########################################################################
#
# Post Hydra-processing
#
########################################################################
###################################################################################
#
# Define the state-space (patchy habitat - Mike Meredith approach
#
###################################################################################
# Read in shapefiles of the plots #
library(raster)

JAMplots <- raster::shapefile("Spatial_Layers/JAMplots_all.shp")

# identify plots as shapefiles #
MM <- JAMplots[JAMplots$PLOT == "MM",]
LI <- JAMplots[JAMplots$PLOT == "LI",]
MI <- JAMplots[JAMplots$PLOT == "MI",]
LII <- JAMplots[JAMplots$PLOT == "LII",]
MII <- JAMplots[JAMplots$PLOT == "MII",]

# create a single list of the plots #
plots <- list(MM,LI,MI,LII,MII)

# generate convex hull around each plot 
plotpolys <- lapply(plots,rgeos::gConvexHull)

plotspecificarea <- lapply(plotpolys,FUN = function(x){rgeos::gArea(x)/10000})

# 5 x 5 m raster of plot extent
JAMraster <- raster(extent(JAMplots),res=c(5,5))

# set values to 0 - i.e., not sampled 
JAMraster[]<-0

# make rasters for each plot and give them values 1 - 5
plotrasters <- vector("list",5)
for(i in 1:5){
plotrasters[[i]] <- rasterize(plotpolys[[i]],JAMraster)
plotrasters[[i]][plotrasters[[i]]==1]<-i
}

# Calculate area of the 5 plots in hectares
(area(plotpolys[[1]])+area(plotpolys[[2]])+area(plotpolys[[3]])+area(plotpolys[[4]])+area(plotpolys[[5]]))/10000

# habmatrix 
study.area <- sum(stack(c(JAMraster,plotrasters)),na.rm = TRUE)

# create matrix
habmat <- as.matrix(study.area)

# Flip the plot raster to match the net locations 
study.flip <- flip(raster(habmat,xmn=0, xmx=ncol(habmat), ymn=0, ymx=nrow(habmat)),"y")
habmat <- as.matrix(study.flip)


M2 <- readRDS("Results/AMRE_survival.rds")

str(M2$sims.list$S)
sims <- 750
inds <- 800
years <- 10

plot_yr <- array(NA,c(sims,inds,years))
plot_density <- array(NA,c(sims,years,5))

area <- c(4.085429,5.474032,4.090979,6.484236,4.914097)

for(i in 1:sims){
for(y in 1:years){
for(n in 1:800){
plot_yr[i,n,y] <- habmat[floor(M2$sims.list$S[i,n,2,y]+1),floor(M2$sims.list$S[i,n,1,y]+1)]*M2$sims.list$z[i,n,y]
}
for(p in 1:5){
plot_density[i,y,p] <- length(which(plot_yr[i,,y]==p))/area[p]
}
}
}

densities <- apply(plot_density,c(2,3),quantile,prob = c(0.025,0.5,0.975))

plot(NA,ylim = c(0,5),xlim = c(1,10),type = "o",pch = 19)
for(p in 1:5){
points(densities[2,,p],type = "o", pch = 19, col = p)
segments(x0 = 1:10,x1 = 1:10, 
         y0 = densities[1,,p],
         y1 = densities[3,,p])
}