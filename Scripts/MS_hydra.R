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

win.data <- readRDS("Data/win_data.rds")
Sst <- readRDS("Data/Sst.rds")
z_start <- readRDS("Data/z_start.rds")

inits <- function(){list(z = z_start,
                  S = Sst)}

params <- c("psi","sp.phi","backphi","intercept", "beta", "sigma", "N","R","D")


n.chains <- 3
n.iter <- 10
n.burn <- 5
n.thin <- 1
n.cores <- 3


Sys.time()
a<-Sys.time()
M1 <- jagsUI::jags.basic(model = "Models/MS_multiseason_scr.txt",
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

saveRDS(M2,"multispecies_survival.rds")
saveRDS(m1,"multispecies_survival_coda.rds")

