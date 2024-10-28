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

win.data <- readRDS("Data/BTBW_windata.rds")

Sst <- readRDS("Data/BTBW_Sst.rds")
z_start <- readRDS("Data/BTBW_z_start.rds")


inits <- function(){list(z = z_start,
                  S = Sst)}


params <- c("psi","phi","intercept","beta.phi","beta","gamma", "SIGMA", "N","R","D","S","z","Nplots","Dplots")


n.chains <- 3
n.iter <- 1000
n.burn <- 500
n.thin <- 2
n.cores <- 3


Sys.time()
a<-Sys.time()
M1 <- jagsUI::jags.basic(model = "Models/multiseason_scr_PhiCov_bvn.txt",
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

saveRDS(M2,"BVN_BTBW_survival.rds")
saveRDS(M1,"BVN_BTBW_survival_coda.rds")

