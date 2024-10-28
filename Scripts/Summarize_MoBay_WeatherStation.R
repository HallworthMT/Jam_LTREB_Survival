
MoBay <- read.csv("Data/MontegoBay_WeatherStation.csv")

MoBay$date <- as.POSIXct(MoBay$DATE, format = "%Y-%m-%d")

MoBay$Year <- format(MoBay$date,"%Y")
MoBay$Month <- format(MoBay$date, "%m")


MonthTemp <- tapply(MoBay$TAVG, list(MoBay$Year,MoBay$Month), FUN = mean, na.rm = TRUE)
MonthPrecip <- tapply(MoBay$PRCP, list(MoBay$Year,MoBay$Month), FUN = mean, na.rm = TRUE)



T<-apply(MonthTemp[10:nrow(MonthTemp),1:3],1,mean)
P<-apply(MonthPrecip[10:nrow(MonthTemp),1:3],1,sum,na.rm = TRUE)
S<-SOI[,14]

plot(P~S)
# Get SOI data 

SOI <- as.matrix(read.table("https://crudata.uea.ac.uk/cru/data/soi/soi.dat",header = FALSE))

colnames(SOI) <- c("Year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec","Annual")

SOI[SOI==-99.99]<-NA

SOI <- SOI[SOI[,1] > 2006,]


# write.csv(SOI,"Data/SOIindex.csv",row.names = FALSE)

