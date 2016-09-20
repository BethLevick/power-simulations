### Code to run analysis from paper on GLMM simulations and power calculations ###
### Bethany Levick, University of Liverpool 2016 ###
#####################################################################################
### Set up ###
## setwd
setwd("C:/Users/Bethany/Dropbox/PhD/methods-paper/power-simulation")
## load functions
source( "methods-functions.r" )
## load set up
source( "methods-setup.r" )
## This will produce warnings due to some data coercion
#####################################################################################
## Run initial model using the data as collected
## colonisation
response <- col.r
## km dists
colkm2 <- glmmadmb( response ~ (hvals) + (sec.occ) + (1|sector),
	family="binomial", data=coltab )
#####################################################################################
## Alter coefficient value without changing intercept
dataset <- matrix(ncol=5, c(NA,NA,NA,NA,NA))
out <- matrix( nrow=1, ncol=5, NA )
simresp <- matrix(ncol=2, nrow=1, NA)
simssameint <- multiSim(coefvec=seq(0.03,0.3,length.out=10), object=colkm2, resp.data=col.r, size=10, use.data=coltab, nsim=10, 
	newTrials=FALSE, trialsvec=0, 
	newint=FALSE, newcoef=TRUE, coefname="hvals")
	
save(out, file="sameintout.RData" )
save( dataset, file="sameintdata.RData" )
save( simresp, file="sameintresp.RData" )

## plot of power
plot( simssameint[,2], simssameint[,4], xlab="Coefficient", ylab="P Value" )

#####################################################################################
## Fit intercept values
dataset <- matrix(ncol=5, c(NA,NA,NA,NA,NA))
out <- matrix( nrow=1, ncol=5, NA )
simresp <- matrix(ncol=2, nrow=1, NA)
simsfitint <- multiSim(coefvec=seq(0.03,0.3,length.out=10), object=colkm2, resp.data=col.r, size=10, use.data=coltab, nsim=10, 
	newTrials=FALSE, trialsvec=0, 
	newint=TRUE, newcoef=TRUE, coefname="hvals")
	
save(out, file="fitintout.RData" )
save( dataset, file="fitintdata.RData" )
save( simresp, file="fitintresp.RData" )

#####################################################################################
## Fit intercept values and increase sample size
colnames(newtab)[3] <- "sector" 
newtab$sector <- as.factor(newtab$sector)
## TODO FIX THE ABOVE BEFORE RUNNING
dataset <- matrix(ncol=5, c(NA,NA,NA,NA,NA))
out <- matrix( nrow=1, ncol=5, NA )
simresp <- matrix(ncol=2, nrow=1, NA)	
simsfitintnewrows <- multiSim(coefvec=seq(0.03,0.3,length.out=10), object=colkm2, resp.data=col.r, size=10, use.data=newtab, nsim=10, 
	newTrials=TRUE, trialsvec=as.vector(newtrialsr), newint=TRUE, 
	intval=0, newcoef=TRUE,
	coefname="hvals")
	
save(out, file="morerowsout.RData" )
save( dataset, file="morerowsdata.RData" )
save( simresp, file="morerowsresp.RData" )







