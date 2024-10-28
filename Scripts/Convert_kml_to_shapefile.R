kmlfiles <- list.files("Jamaica Plots/Jamaica Plots", pattern = ".kml", full.names = TRUE)

getCoords <- function(kmlfile){
# read in the kml file #
kml <- readLines(kmlfile)

# get the names of the points #
getname <- "\t\t\t\t\t<SimpleData name=\"Comment\"> *([^<]+?)</SimpleData>"

# get the needed info #
temp <- kml[grep(getname,kml)]
temp1 <- gsub(pattern = "\t\t\t\t\t<SimpleData name=\"Comment\">","",temp)
temp2 <- gsub(pattern = "</SimpleData>","",temp1)

# extract the coordinates using maptools #
tkml <- maptools::getKMLcoordinates(kmlfile=kmlfile, ignoreAltitude=T)

coords <- t(sapply(tkml,unlist))

out <- data.frame( LOCATION = temp2,
                  LONG = coords[,1], 
                  LAT = coords[,2])
return(out)
}

# Run the function #
plots <- lapply(kmlfiles,getCoords)

# get names of the plots to make shapefile with #
plotnames <- gsub(list.files("Jamaica Plots/Jamaica Plots", pattern = ".kml"),pattern = ".shp.kml",replacement = "")

# Add plot names and reorder columns #
for(i in seq(plotnames)){
plots[[i]]$PLOT <- rep(plotnames[i],nrow(plots[[i]]))
plots[[i]] <- plots[[i]][,c("PLOT","LOCATION","LONG","LAT")]
}

cplots <- do.call("rbind",plots)

plots <- sp::SpatialPointsDataFrame(sp::SpatialPoints(cbind(cplots$LONG,cplots$LAT)),cplots)

raster::crs(plots) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

#raster::shapefile(plots,"Jamaica Plots/Jamaica Plots/Exp_Plots.shp")