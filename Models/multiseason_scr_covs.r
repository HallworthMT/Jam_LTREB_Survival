model{

#########################################################################################
#
# Priors
#
#########################################################################################
# psi - detection 
psi ~ dbeta(1,1) 

# intercept - intercept 
intercept ~ dnorm(0,0.01)

# sigma - standard deviation 
sigma ~ dunif(0,500)

# beta - distance related decay for detection
beta <- 1/(2*sigma*sigma)

# betas for survival models
for(nc in 1:5){
beta.phi[nc] ~ dnorm(0,0.01)
}

# precision measure for movements in X & Y directions
tauXY <- 1/pow(sigmaXY,2)

# standard deviation for movements in X & Y directions
sigmaXY ~ dunif(0,500)

for(y in 1:nyears){
# Prior for year specific recruitment probability
gamma[y] ~ dbeta(1,1) 

# Prior for year specific survival
#phi[y] ~ dbeta(1,1) 

# intercept of phi model 
#lphi[y] <- log(phi[y]/(1-phi[y]))

# Linear model for survival covariates
#phi.s[y] <- (1 / (1+exp(-(lphi[y]+beta.phi[1]*precip[y]+
#                                  beta.phi[2]*temp[y]+
#                                  beta.phi[3]*(precip[y]*temp[y])+
#                                  beta.phi[4]*soi[y]+
#								   beta.phi[5]*soi2[y]))))
# Linear model for survival covariates
phi.s[y] <- (1 / (1+exp(-(lphi+beta.phi[1]*precip[y]+
                               beta.phi[2]*temp[y]+
                               beta.phi[3]*(precip[y]*temp[y])+
                               beta.phi[4]*soi[y]+
							   beta.phi[5]*soi2[y]))))
} # END years

# intercept of phi model 
lphi <- log(phi/(1-phi))

# Prior for year specific survival
phi ~ dbeta(1,1) 

##################################################################################
#
# Likelihood 
#
##################################################################################
# FIRST YEAR

for(i in 1:M){
# Survival at time 1 is bernouli of detection
z[i,1]~dbern(psi)

# B indicates if individual is available to be recruited at time 1
B[i,1]<-(1-z[i,1])

# Activity centers for first year
S[i,1,1] ~ dunif(0,xlimit)
S[i,2,1] ~ dunif(0,ylimit)

# SECOND + YEARS locations
for(y in 2:nyears){
# Activity centers 

# this works but is somewhat unrealistic 
S[i,1,y] ~ dnorm(S[i,1,y-1],tauXY)T(0,xlimit) # truncate to 0 and xlimit
S[i,2,y] ~ dnorm(S[i,2,y-1],tauXY)T(0,ylimit) # truncate to 0 and ylimit

# Alive state at y is dependent on phi and gamma
 mu[i,y] <- (phi.s[y]*z[i,y-1])+(gamma[y]*B[i,y-1])
# mu[i,y] <- (phi.s[y]*z[i,y-1])+(gamma*B[i,y-1])

# alive or dead bernouli trial 
z[i,y] ~ dbern(mu[i,y])

# sum over alive states from 1 to y
still_alive[i,y] <- sum(z[i,1:y]) 

# B indicates if individual is available to be recruited at time y
B[i,y] <- 1-step(still_alive[i,y]-1)
} # END nyears

for(y in 1:nyears){
# Is the location within the plots 
# habmat = 1 if 'on' plot and habmat = 0 if 'off' plot
plotLoc[i,y] <- habmat[trunc(S[i,2,y]+1), trunc(S[i,1,y]+1)]

# Calculate the distance from activity center to each net
d2[i,1:nnets,y] <- pow(S[i,1,y]-NetLoc[1:nnets,1],2) + pow(S[i,2,y]-NetLoc[1:nnets,2],2)

for(j in 1:nmonths){
 for(n in 1:nnets){
 logitp[i,j,n,y] <- exp(intercept - beta*d2[i,n,y])*z[i,y]
 capturep[i,j,n,y] <- logitp[i,j,n,y]/(1+sum(logitp[i,j,,y]))
 } # end nets
 
# the last cell probability indicates not captured
capturep[i,j,nnets+1,y] <- 1-sum(capturep[i,j,1:nnets,y])
} # end months
} # end years

# Estimate number of recruits per year. 
# This ensures that an individual can only recruit into the population once
r[i,1] <- z[i,1]
r[i,2]<-(1-z[i,1])*z[i,2]
r[i,3]<-(1-z[i,1])*(1-z[i,2])*z[i,3]
r[i,4]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*z[i,4]
r[i,5]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*z[i,5]
r[i,6]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*(1-z[i,5])*z[i,6]
r[i,7]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*(1-z[i,5])*(1-z[i,6])*z[i,7]
r[i,8]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*(1-z[i,5])*(1-z[i,6])*(1-z[i,7])*z[i,8]
r[i,9]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*(1-z[i,5])*(1-z[i,6])*(1-z[i,7])*(1-z[i,8])*z[i,9]
r[i,10]<-(1-z[i,1])*(1-z[i,2])*(1-z[i,3])*(1-z[i,4])*(1-z[i,5])*(1-z[i,6])*(1-z[i,7])*(1-z[i,8])*(1-z[i,9])*z[i,10]


} # end M

for(i in 1:n0){
for(y in 1:nyears){
for(j in 1:nmonths){
# Ycat is a capture history with net number where bird was captured 
# instead of a 1 and (nnets+1)(81 in this case) if a bird was not
# captured. 
Ycat[i,j,y] ~ dcat(capturep[i,j,,y])
CaptProb[i,j,y]<-sum(capturep[i,j,1:nnets,y])
} # END months
} # END years
} # END n0

for(i in (n0+1):M) {
for(y in 1:nyears){
  zeros[i,y] ~ dbern(sum(capturep[i,1:nmonths,1:nnets,y])*z[i,y]) # capturep = probability of being captured at least once
  } # END nyears
}# END augmented n0+1:M
#############################################################################################
#
# DERIVED PARAMETERS
#
#############################################################################################

# Total area of all study plots
Area <- 25.04877 # hectares
area <- c(4.085429,5.474032,4.090979,6.484236,4.914097)

for(y in 1:nyears){

for(pp in 1:5){
# Determine how many birds that are alive are on plots
Nplots[pp,y] <- sum(equals(plotLoc[1:M,y],pp)*z[1:M,y])
# Determine the density of birds on each plot
Dplots[pp,y] <- Nplots[pp,y]/area[pp]
}

# calculate abundance for each year
N[y] <- sum(z[1:M,y])

# estimate density in each year
D[y] <- sum(Nplots[,y])/Area 

# estimate recruitment in each year
R[y] <- sum(r[1:M,y])

} # END nyears

for(i in 1:n0){
# get the longevity of each individual
Longevity[i] <- sum(z[i,])
}
# END model
}