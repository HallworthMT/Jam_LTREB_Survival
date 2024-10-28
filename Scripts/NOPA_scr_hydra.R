#### -------------------- set path to library dir --------------------------- #### 

library(parallel,lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(coda, lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(grid)
library(lattice,lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(rjags,lib.loc="/home/hallworthm/CostOfRepro/Packages")
library(jagsUI,lib.loc="/home/hallworthm/CostOfRepro/Packages")
#library(abind,lib.loc="/home/hallworthm/CostOfRepro/Packages")
#library(R2jags,lib.loc="/home/hallworthm/CostOfRepro/Packages")
#library(runjags, lib.loc="/home/hallworthm/CostOfRepro/Packages")

#### -------------------- set path to library dir --------------------------- ####

my.path <- "/home/hallworthm/CostOfRepro/Packages"

.libPaths(c(my.path,.libPaths()))

win.data <- readRDS("Data/NOPA_windata.rds")

Sst <- readRDS("Data/NOPA_Sst.rds")
z_start <- readRDS("Data/NOPA_z_start.rds")


inits <- function(){list(z = z_start,
                  S = Sst)}


params <- c("psi","phi","intercept","beta.phi","beta","gamma", "sigma", "N","R","D","S","z","Nplots","Dplots")


n0 <- max(which(apply(win.data$Ycat, 1,FUN=function(x){all(x==81)})==FALSE))

Ycat <- win.data$Ycat[1:n0,,]

zeros <- array(NA,c(win.data$M,win.data$nyears))

zeros[(n0+1):win.data$M,]<-0

win.data$Ycat <- Ycat
win.data$n0 <- n0
win.data$zeros <- zeros 

n.chains <- 3
n.iter <- 10000
n.burn <- 5000
n.thin <- 5
n.cores <- 3


Sys.time()
a<-Sys.time()
M1 <- jagsUI::jags.basic(model = "Models/multiseason_scr_PhiCov_streamline.txt",
             data = win.data,
             inits = inits,
             parameters.to.save = params,
             n.chains = n.chains,
             n.iter = n.iter,
             n.burn = n.burn,
             n.thin = n.thin,
             n.cores = n.cores,
             parallel = TRUE,
             seed = 12345,
             DIC = FALSE)
Sys.time()-a
Sys.time()


source("Scripts/sims.list.R")

M2 <- process.output(M1)

saveRDS(M2,"NOPA_survival_streamline.rds")
saveRDS(M1,"NOPA_survival_streamline_coda.rds")

