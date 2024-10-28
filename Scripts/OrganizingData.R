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


SpeciesList <- vector("list",length(spp))

names(SpeciesList) <- spp

temp <- subset(LTREB, BNDCREW == "LTREB")


Sys.time()
a <- Sys.time()
for(s in seq(spp)){

temp.spp <- subset(LTREB,SPECIES == spp[s] & BNDCREW == "LTREB")

inds <- unique(temp.spp$BAND)

# Create empty vectors to store data #
FirstCap <- LastCap <- Sex <- InitAge <- rep(NA,length(inds))

# give vectors names (Bandnumbers) #
names(FirstCap) <- names(LastCap) <- names(Sex) <- names(InitAge) <- inds

# Create empty arrays to store data #
ch.spp.mon <- array(0,c(length(inds),length(month),length(years)))
ch.spp.mon.loc <- array(NA,c(length(inds),length(month),length(years),2))
ch.yearly <- array(NA,c(length(inds),length(years)))
manipulated <- array(0,c(length(inds),length(month),length(years)))
habitat <- array(NA,c(length(inds),length(month),length(years)))
hab.yearly <- array(NA,c(length(inds),length(years)))

# Give array names #
dimnames(ch.spp.mon)[[1]] <- dimnames(ch.spp.mon.loc)[[1]] <- dimnames(ch.yearly)[[1]] <- dimnames(manipulated)[[1]] <- inds
dimnames(ch.spp.mon)[[2]] <- dimnames(ch.spp.mon.loc)[[2]] <- dimnames(manipulated)[[2]] <- month
dimnames(ch.spp.mon)[[3]] <- dimnames(ch.spp.mon.loc)[[3]] <- dimnames(manipulated)[[3]] <- dimnames(ch.yearly)[[2]] <- years

# Loop through individuals to create monthly capture histories #
for(i in seq(inds)){

# subset the dataset for each individual 
temp.ind <- subset(temp.spp,BAND == inds[i])

# determine which year the individual was first captured #
FirstCap[i] <- which(years == temp.ind$YEAR[1]) 

# determine which year the individual was last captured #
LastCap[i] <- which(years == max(temp.ind$YEAR,na.rm = TRUE))

# Sex and Age of individual at inital capture #
Sex[i] <- temp.ind$SEX[1]
InitAge[i] <- temp.ind$AGE[1]

# For each capture in the dataset #
for(n in seq(nrow(temp.ind))){

# enter 1 if captured # 
ch.spp.mon[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n])] <- 1

# fill in capture locations #
ch.spp.mon.loc[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n]),1] <- temp.ind$LONG[n]
ch.spp.mon.loc[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n]),2] <- temp.ind$LAT[n]

#if(temp.ind$CAPLOC[n] %in% LTREBnets$PlotLoc){
#xycoords <- sp::spsample(LTREBnets[LTREBnets$PlotLoc == temp.ind$CAPLOC[n],], n = 1, type = "random")
#ch.spp.mon.loc[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n]),1] <- xycoords@coords[1,1]
#ch.spp.mon.loc[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n]),2] <- xycoords@coords[1,2]
#}else{
#ch.spp.mon.loc[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n]),1] <- NA
#ch.spp.mon.loc[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n]),2] <- NA
#}

manipulated[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n])] <- ifelse(temp.ind$MANIPULATED[n]== "Y",1,0)

habitat[i,which(month == temp.ind$CAPMON[n]),which(years == temp.ind$YEAR[n])] <- ifelse(temp.ind$HABITAT[n] == "L",2,1)

} # end captures

# Change Sex 
# Sex M = 5, F = 1, U = 7 

Sex[Sex == 5] <- 2
Sex[Sex %in% c(6,7)] <- 3
Sex[is.na(Sex)] <- 3

# Change Initial Age 
# AGE
AHY = 1 AHY? = 2 ASR = 3 ASY = 4 ASY = 5  ASY? = 6 HY = 7 HY? = 8 SY = 9 SY  = 10 TY = 11 U = 12

InitAge[InitAge %in% c(7,8,9,10)]<-1
InitAge[InitAge %in% c(1,2,3,4,5,6)]<-2
InitAge[InitAge %in% c(11)] <- 3
InitAge[InitAge %in% c(12)] <- 4

# Yearly capture history #
# Some captures were made in Sept-December 
# those captures correspond to the non-breeding season 
# of Jan- May. Therefore, I summed the captures from
# Sept-December in year t-1 and Jan-May in year t 

#for(y in 2:length(years)){
#ch.yearly[i,y] <- sum(ch.spp.mon[i,6:9,y-1],ch.spp.mon[i,1:5,y])
#hab.yearly[i,y] <- max(c(habitat[i,6:9,y-1],habitat[i,1:5,y]),na.rm = TRUE)
#} # end years 

for(y in 1:length(years)){
ch.yearly[i,y] <- sum(ch.spp.mon[i,1:2,y])
hab.yearly[i,y] <- max(habitat[i,1:2,y],na.rm = TRUE)
} # end years 

for(n in FirstCap[i]:LastCap[i]){
hab.yearly[i,n] <- ifelse(hab.yearly[i,n]== "-Inf",hab.yearly[i,n-1],hab.yearly[i,n])
}
} # end individuals

FoundLocs <- array(NA,dim(ch.spp.mon.loc[,,,1]))
FoundLocs[,1,] <- apply(ch.spp.mon.loc[,,,1],c(1,3),min,na.rm = TRUE)
FoundLocs[,2,] <- apply(ch.spp.mon.loc[,,,2],c(1,3),min,na.rm = TRUE)

#which(ch.yearly == 1 & is.na(FoundLocs[,1,]), arr.ind = TRUE)

hab.yearly[hab.yearly == "-Inf"] <- NA

SpeciesList[[s]] <- list(FirstCap = FirstCap,
                         LastCap  = LastCap,
                         ch.monthly = ch.spp.mon,
                         monthly.locs = ch.spp.mon.loc,
                         FoundLocs = FoundLocs,
                         ch.yearly = ch.yearly,
                         manipulated = manipulated,
                         InitAge = InitAge, 
                         Sex = Sex,
                         habitat.m = habitat,
                         hab.yearly = hab.yearly)
} # end species

Sys.time()
Sys.time()-a


migs <- c("AMRE","BAWW","NOWA","OVEN","PRAW","PROW","SWWA")

mig.res <- ifelse(spp %in% migs, "Migrant","Resident")

# Get data ready for Multiple Species Survival analysis #

nbirds <- unlist(lapply(SpeciesList, FUN = function(x){nrow(x$ch.yearly)}))

# Save only AMRE that were first banded in 2008 or later #
#nbirds[1] <- length(which(SpeciesList[[1]]$FirstCap >= 21))

# make empty capture history for each species #
ch.spp <- array(0,c(max(nbirds),10,length(nbirds)))
hab.spp <- array(NA,c(max(nbirds),10,length(nbirds)))
FoundLocs <- array(NA,c(max(nbirds),10,2,length(nbirds)))
# fill the species capture history #
for(i in seq(nbirds)){
ch.spp[1:nbirds[i],,i] <- SpeciesList[[i]]$ch.yearly[,1:10]
hab.spp[1:nbirds[i],,i] <- SpeciesList[[i]]$hab.yearly[,1:10]
FoundLocs[1:nbirds[i],1:10,1,i] <- SpeciesList[[i]]$FoundLocs[,1,1:10]
FoundLocs[1:nbirds[i],1:10,2,i] <- SpeciesList[[i]]$FoundLocs[,2,1:10]
}


which(ch.spp == 1 & is.na(FoundLocs[,,1,]),arr.ind = TRUE)
# Change all values > 1 to 1 #
ch.spp[ch.spp > 1] <- 1

# Re-center the FirstCap information 
FirstCap <- LastCap <- array(NA,c(max(nbirds),length(nbirds)))
for(i in seq(nbirds)){
FirstCap[1:nbirds[i],i] <- apply(ch.spp[1:nbirds[i],,i],1,FUN = function(x){min(which(x!=0))})
LastCap[1:nbirds[i],i] <- apply(ch.spp[1:nbirds[i],,i],1,FUN = function(x){max(which(x!=0))})
}

for(s in seq(nbirds)){
for(i in 1:nbirds[s]){
for(n in FirstCap[i,s]:LastCap[i,s]){
hab.spp[i,n,s] <- ifelse(is.na(hab.spp[i,n,s]),hab.spp[i,n-1,s],hab.spp[i,n,s])
}
for(n in LastCap[i,s]:length(hab.spp[i,,s])){
hab.spp[i,n,s] <- ifelse(is.na(hab.spp[i,n,s]),hab.spp[i,n-1,s],hab.spp[i,n,s])
}
}
}

sex1 <- sapply(SpeciesList,FUN = function(x){x$Sex})
sex <- array(NA,c(max(nbirds),length(nbirds)))
for(i in seq(nbirds)){
sex[1:nbirds[i],i] <- sex1[[i]]
}
sex[sex==4]<-3

#####################################################################################################
#
# Find individuals that were captured but don't have location data associated with the capture
#
#####################################################################################################
which(win.data$CH[,,7]==1 & is.na(win.data$FoundLocs[,,1,7]),arr.ind = TRUE)




# READ IN NDVI and PRECIP data

NDVI <- read.csv("Data/NDVImonthly2007_2017.csv")
years <- 2007:2017
months <- 1:12
NDVIarray <- array(NA,c(length(years),length(months)))
for(i in seq(years)){
for(m in seq(months)){
tmp <- NDVI$mean.ndvi[NDVI$year == years[i] & NDVI$month == months[m]]
NDVIarray[i,m] <- ifelse(length(tmp)>0,tmp,NA)
}
}

AprilNDVI <- NDVIarray[2:nrow(NDVIarray),4]
DeltaNDVI <- rep(NA,10)
for(y in 2:11){
DeltaNDVI[y-1] <- NDVIarray[y-1,10]-NDVIarray[y,4]
}


# Precip #
Precip <- read.csv("Data/PrecipMonthly2007_2017.csv")
years <- 2007:2017
months <- 1:12
Preciparray <- array(NA,c(length(years),length(months)))
for(i in seq(years)){
for(m in seq(months)){
tmp <- Precip$MonthlyPrecip[Precip$year == years[i] & Precip$month == months[m]]
Preciparray[i,m] <- ifelse(length(tmp)>0,tmp,NA)
}
}



AprilPrecip <- Preciparray[2:nrow(Preciparray),4]
SumPrecip <- rep(NA,10)
for(y in 2:11){
SumPrecip[y-1] <- sum(c(Preciparray[y-1,10:12],Preciparray[y,1:4]),na.rm = TRUE)
}
SumPrecip[10] <- mean(SumPrecip[1:9])

EffortExtent <- lapply(EffortPolys,extent)

grid.study.boundaries <- do.call(rbind,lapply(EffortExtent, function(x){
                                               d <- cbind(x[1],x[2],x[3],x[4])
                                               return(d)}))


win.data <- list(CH = ch.spp,
                 First = FirstCap,
                 n = nbirds,
                 sex = sex,
                 nspp = length(nbirds),
                 n.occasions = dim(ch.spp)[2],
                 migstatus = ifelse(mig.res == "Migrant",1,2),
                 migrants = which(mig.res == "Migrant"),
                 residents  = which(mig.res == "Resident"),
                 ngrids = nrow(grid.study.boundaries),
                 grid.study.boundaries = grid.study.boundaries,
                 xlim = extent(JAMplots)[1:2],
                 ylim = extent(JAMplots)[3:4],
                 FoundLocs = FoundLocs,
                 NDVIstart = scale(NDVIarray[,10]),
                 NDVImar = scale(AprilNDVI),
                 NDVIdelta = scale(DeltaNDVI),
                 PPTmar = scale(AprilPrecip),
                 PPTsum = scale(SumPrecip),
                 PPTstart = scale(Preciparray[,10]),
                 ncovs = 6)

win.data$PPTstart[is.na(win.data$PPTstart)]<-0
win.data$FoundLocs[is.infinite(win.data$FoundLocs)]<-NA
#### ------------------ Function for surival state init ------------------- ####

known.state.cjs <- function(x,FirstCap){
   state <- x
   for (i in 1:dim(x)[1]){
      n1 <- min(which(x[i,]==1))
      n2 <- max(which(x[i,]==1))
      state[i,n1:n2] <- 1
      state[i,1:FirstCap[i]]<- NA
      state[i,n1] <- NA
      }
   state[state==0] <- NA
   return(state)
   }


#### ------------------ Function for inits  ------------------- ####
z <- array(NA,c(max(nbirds),10,length(nbirds)))

for(i in seq(nbirds)){
z[1:nbirds[i],,i] <- known.state.cjs(x = ch.spp[1:nbirds[i],,i], FirstCap = FirstCap[1:nbirds[i],i])
}

known.state.cjs(ch.spp[1:nbirds[1],,1],FirstCap[1:nbirds[1],1])
# Set Initial values
inits <- function(){
              list(z = z)
  
}

library(jagsUI)


n.chains <- 3
n.thin <- 1
n.adapt <- 10
n.iter <- 10
n.burnin <- 5

Sys.time()
a <- Sys.time()
M <- jagsUI::jags(model.file = "Models/MS_survival_spatial.txt",
             data = win.data, 
             inits = inits, 
             n.iter = n.iter, 
             n.thin = n.thin, 
             n.chains = n.chains,
             n.burnin = n.burnin,
             parameters.to.save = c("BetaHyperPrior","mean.phi.t","p.mean","res.phi","mig.phi","beta","Life","sp.phi","backphi","hyperphi"),
             parallel = TRUE)
Sys.time()
Sys.time() - a


plot(apply(M$mean$Life,2,mean,na.rm = TRUE))
plot(density(M$mean$Life[,1]))


M$mean$beta
plot(NA, ylim = c(0,1), xlim = c(1,2))
points(1,M$mean$backphi[1])
points(2,M$mean$backphi[2])
#par(bty = "l",mfcol = c(2,2),mar = c(4,4.1,1,1))
# Migrant vs Resident # Males
#plot(M$mean$mig.phi[2,1,],type = "b", ylim = c(0.2,0.8),xlim = c(1,11),pch = 19, cex = 1.2,
#     ylab = expression(phi[Male]),xlab = "",yaxt = "n",xaxt = "n",main = "Mangrove",cex.lab = 2)
#axis(2,las = 2)
#axis(1, at = 1:10, las = 2, labels = c("'07-'08","'08-'09","'09-'10","'10-'11","'11-'12","'12-'13","'13-'14","'14-'15","'15-'16","'16-'17"), cex.lab = 0.8)
#points(1.2:10.2,M$mean$res.phi[2,1,],type = "b", pch = 21, cex =1.2)
#segments(x0 = 1:10, x1 = 1:10, y0 = M$q2.5$mig.phi[2,1,], y1 = M$q97.5$mig.phi[2,1,])
#segments(x0 = 1.2:10.2, x1 = 1.2:10.2, y0 = M$q2.5$res.phi[2,1,], y1 = M$q97.5$res.phi[2,1,])
#legend(2.8,0.8,legend = c("Resident","Migrant"),pch = c(21,19),bty = "n",ncol = 2)

#plot(M$mean$res.phi[1,1,],type = "b", ylim = c(0.2,0.8),xlim = c(1,11),pch = 2, cex = 1.2,
#     ylab = expression(phi[Female]),xlab = "",yaxt = "n",xaxt = "n",cex.lab = 2)
#axis(2,las = 2)
#axis(1, at = 1:10, las = 2, labels = c("'07-'08","'08-'09","'09-'10","'10-'11","'11-'12","'12-'13","'13-'14","'14-'15","'15-'16","'16-'17"), cex.lab = 0.8)
#segments(x0 = 1.2:10.2, x1 = 1.2:10.2, y0 = M$q2.5$mig.phi[1,1,], y1 = M$q97.5$mig.phi[1,1,])
#segments(x0 = 1:10, x1 = 1:10, y0 = M$q2.5$res.phi[1,1,], y1 = M$q97.5$res.phi[1,1,])
#points(1.2:10.2,M$mean$mig.phi[1,1,],type = "b", pch = 17, cex =1.2, col = "black")


# Migrant vs Resident # Males - LOGWOOD
#plot(M$mean$mig.phi[2,2,],type = "b", ylim = c(0.2,0.8),xlim = c(1,11),pch = 19, cex = 1.2,
#     ylab = "",xlab = "",yaxt = "n",xaxt = "n",main = "Logwood")
#axis(2,las = 2)
#axis(1, at = 1:10, las = 2, labels = c("'07-'08","'08-'09","'09-'10","'10-'11","'11-'12","'12-'13","'13-'14","'14-'15","'15-'16","'16-'17"), cex.lab = 0.8)
#points(1.2:10.2,M$mean$res.phi[2,2,],type = "b", pch = 21, cex =1.2)
#segments(x0 = 1:10, x1 = 1:10, y0 = M$q2.5$mig.phi[2,2,], y1 = M$q97.5$mig.phi[2,2,])
#segments(x0 = 1.2:10.2, x1 = 1.2:10.2, y0 = M$q2.5$res.phi[2,1,], y1 = M$q97.5$res.phi[2,2,])
#legend(2.8,0.8,legend = c("Resident","Migrant"),pch = c(21,19),bty = "n",ncol = 2)

# Migrant vs Resident # Females - LOGWOOD
#plot(M$mean$res.phi[1,2,],type = "b", ylim = c(0.2,0.8),xlim = c(1,11),pch = 2, cex = 1.2,
#     ylab ="",xlab = "",yaxt = "n",xaxt = "n")
#axis(2,las = 2)
#axis(1, at = 1:10, las = 2, labels = c("'07-'08","'08-'09","'09-'10","'10-'11","'11-'12","'12-'13","'13-'14","'14-'15","'15-'16","'16-'17"), cex.lab = 0.8)
#segments(x0 = 1.2:10.2, x1 = 1.2:10.2, y0 = M$q2.5$mig.phi[1,2,], y1 = M$q97.5$mig.phi[1,2,])
#segments(x0 = 1:10, x1 = 1:10, y0 = M$q2.5$res.phi[1,2,], y1 = M$q97.5$res.phi[1,2,])
#points(1.2:10.2,M$mean$mig.phi[1,2,],type = "b", pch = 17, cex =1.2, col = "black")










# Function to return m.array format for each species. Function is from 
# Kery and Schaub Bayesian Population Analysis Using WinBUGS #
#
# Using m.array is much faster computationally but you lose the ability
# to model individual movements 

Convert2marray <- function(CaptureHistory){
 ninds <- dim(CaptureHistory)[1]
 n.occasions <- dim(CaptureHistory)[2]
 m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)

for(t in 1:n.occasions){
 m.array[t,1] <- sum(CaptureHistory[,t])
 } # t

for(i in 1:ninds){
 pos <- which(CaptureHistory[i,]!=0)
 g <- length(pos)
  for(z in 1:(g-1)){
    m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]]+1
   } # z
 } # i
for(t in 1:n.occasions){
   m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:(n.occasions+1)])
   } # t

out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
return(out)
}

# Covert individual capture histories to m.array format # 
Convert2marray(SpeciesList[[1]]$ch.yearly)

