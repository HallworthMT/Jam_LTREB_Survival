#### -------------------- set path to library dir --------------------------- #### 

.libPaths("Packages")

(n.iter <- as.numeric(commandArgs()[6]))
(n.burn <- as.numeric(commandArgs()[7]))
(n.thin <- as.numeric(commandArgs()[8]))
(n.adapt <- as.numeric(commandArgs()[9]))
(species <- as.numeric(commandArgs()[10]))

#species <- 18
library(jagsUI)

# Read in data 
winfiles <- list.files("Jamaica_LTREB/Data", pattern = "*_windata.rds", full.names = TRUE)
z_start_files <- list.files("Jamaica_LTREB/Data", pattern = "*_z_start.rds", full.names = TRUE)
Sst_files <- list.files("Jamaica_LTREB/Data", pattern = "*Sst.rds", full.names = TRUE)

#winfiles <- list.files("Data/", pattern = "*_windata_10m.rds", full.names = TRUE)
#z_start_files <- list.files("Data/", pattern = "*_z_start.rds", full.names = TRUE)
#Sst_files <- list.files("Data/", pattern = "*Sst.rds", full.names = TRUE)

species_codes <- gsub(".*Data/\\s*|_windata.rds.*","",winfiles)

species <- 1
(species_code <- species_codes[species])

win.data <- readRDS(winfiles[species])
Sst <- readRDS(Sst_files[species])
z_start <- readRDS(z_start_files[species])

win.data$M <- 600



str(win.data,1)
str(Sst,1)
str(z_start,1)

#params <- c("psi","phi","intercept","beta.phi",
#            "beta","trend","gamma", "sigma",
#            "sigma.beta","tau.beta","sigma.trend",
#            "tau.trend","N","R","D","S","z","Nplots","Dplots",
#			"Longevity")

params <- c("psi","phi","intercept","beta.phi",
            "beta","trend","gamma", "sigma",
            "sigma.beta","tau.beta","sigma.trend",
            "tau.trend","N","R","D","Nplots","Dplots",
			"Longevity","sigmaXY")


n0 <- max(which(apply(win.data$Ycat, 1,FUN=function(x){all(x==81)})==FALSE))

Ycat <- win.data$Ycat[1:n0,,]



win.data$Ycat <- Ycat
win.data$n0 <- n0
win.data$M <- win.data$n0+10
zeros <- array(NA,c(win.data$M,win.data$nyears))
zeros[(n0+1):win.data$M,]<-0
win.data$zeros <- zeros 

n.chains <- 3
n.cores <- 3

n.iter = 20
n.burn = 2
n.thin = 1
n.adapt = 200

inits <- function(){list(z = z_start[1:win.data$M,],
                         S = Sst[1:win.data$M,,],
                         psi = 0.5,
                         phi = runif(0.1,0.9),
                         gamma = runif(10,0.4,1),
                         sigma = runif(1,100,300))}




Sys.time()
a<-Sys.time()
# M1 <- jagsUI::jags(model = "Models/multiseason_scr_PhiCov_streamline.txt",
# M1 <- jagsUI::jags(model = "Models/multiseason_scr_PhiCov_streamline_simplify_n0.txt",
#M1 <- jagsUI::jags(model = "Models/multiseason_scr_reverse_order.txt",
M1 <- jagsUI::jags(model = "Models/multiseason_scr_covs.r",
                   module = c("glm"),
                   data = win.data,
                   inits = inits,
                   parameters.to.save = params,
		       codaOnly = c("N","R","S","D","Nplots","Dplots"),
                   n.chains = n.chains,
                   n.iter = n.iter,
                   n.burn = n.burn,
                   n.thin = n.thin,
                   n.adapt = n.adapt,
                   n.cores = n.cores,
                   parallel = TRUE,
                   seed = 12345,
                   DIC = FALSE,
	             store.data = TRUE)
Sys.time()-a
Sys.time()
hist(M1$sims.list$sigmaXY)
saveRDS(M1,paste0("Jamaica_LTREB/Results/",species_code,"_survival_scr.rds"))

