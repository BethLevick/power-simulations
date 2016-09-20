##################################################################################################
## Functions for GLMM power analysis methods ##
##################################################################################################
## Function for environment set up
### is.installed ###
## check for installed packages
## if required package is not installed, install
## otherwise, require to load package
is.installed <- function(pkg="packagename"){
if( is.element( pkg, installed.packages()[,1] )==FALSE ){
install.package( pkg )
}else{
require(package=pkg, character.only=TRUE)
}
}
##################################################################################################
## Functions for data set up, calculating connectivity values
#### getH ####
## To get H value estimates
## takes b value and the matrix of distances between burrows
getH <- function(bval, dist.mat){
	## empty vector
	hvals <- c(rep(NA,nrow(dist.mat)))
	## for each row, the corresponding vector value in the empty vector
	## is the h value for that burrow
	for(j in 1:nrow(dist.mat)){
		hvals[j] <- sum( exp( (-abs(dist.mat[j,]))/bval ), na.rm=T )
	}
	return(hvals)
} 

### runModel ####
## function to optimise
## b value is the value of b at this iteration
## dist.mat the distance matrix
## response a 2 column matrix with a col of successes and a col of events
## and rand the random effect vector
runModel <- function(bvalue, dist.mat, response, rand){
	#print(bvalue)
	## estimate hvalues for this value of b
	hvalues <- getH( bvalue, dist.mat )

	## run model
	## put predictors together in a data frame (random effect and h values)
	## with correct data types
	tab <- data.frame( h=as.numeric(hvalues), rand=as.factor(rand) )
	## then run GLMM through ADMB function
	mod <- glmmadmb( response ~ (h) + (1|rand), family="binomial", data=tab )

	## return neg log lik
	## just logLik(mod) returns a logLik object
	## using as.numeric() converts this to just a number variable
	return( -as.numeric(logLik(mod)) ) 
}
##########################################################################
## Functions to generate simulations and perform analysis
### coefInt ###
## wrapper to iterate through coefficient values, find best fit intercept and find resulting p value
## coefvec is a vector of potential coefficient values
## int is a vector of potential intercept values
## resp.data is the original response data, a matrix of successes and trials
## size is number of intercept values to use
## use.data is the original predictor data
coefInt <- function(coefvec, object, resp.data, size, use.data){

outfr <- data.frame( coef=coefvec, int=c(rep(NA,length(coefvec))), pval=c(rep(NA,length(coefvec))),
	stde=c(rep(NA,length(coefvec))), glmp=c(rep(NA,length(coefvec))) )
#sendout <<- outfr
for(i in 1:length(coefvec) ){
#print( paste("coef at position", i, sep=": "))
## find the best fitting intercept at each coefficient value
newint <- fitInt( object, resp.data, size, coefval=coefvec[i], use.data )[[3]]
## store the intercept and the p value of the hvals term at this coefficient
outfr$int[i] <-  newint
## generate new response data at current coefficient with best fit intercept
response <<- simulate.con(object, use.data, newcoef=TRUE,coefval=coefvec[i],newint=TRUE, intval=newint, type="response", newTrials=TRUE, trialsvec=as.vector(newtrialsr))
#print(head(response))
## then re fit the model
## TODO find way to pass the function as an argument
fit <- glmmadmb( response ~ (hvals) + (sec.occ) + (1|sector),
	family="binomial", data=use.data )
## again find a way to reference the correct p value rather than giving a direct location
outfr$pval[i] <- summary(fit)$coefficients[2,4]
## save standard error of hvals term
outfr$stde[i] <- summary(fit)$coefficients[2,2]
## get pval from ANOVA
## doesnt appear to work with glmmadmb
##remove the model
rm(fit)
## rebuild as glm
fit <- glm( response ~ hvals * sec.occ, family="binomial" )
## save p val
outfr$glmp[i] <- summary(fit)$coefficients[2,4]
}

return(outfr)

}

#################################################################
### fitInt ###
## function to find best fit intercept value at current coefficient value
## object is model object
## newTrials is a boolean stating whether the sample size is to be extended
## if newTrials is TRUE then a vector of new numbers of trials should be passed as trialsvec
fitInt <- function(object, resp.data, size, coefval, use.data, newTrials, trialsvec){

## find the proportion of successes over all the events in the original data
ps.orig <- sum( resp.data[,1] )/sum( resp.data[,2] )

## for a range of intercept values
## generate new data at current coefficient value and inspect properties of the data
## range of intercept values centred at value in original model
## extract intercept value
int <- as.numeric(object$b[names(object$b)=="(Intercept)"])
## vector of size given above
## mean/sd set based on a sample of 100 to give a mix of pos and neg values (ie so would end up as mean&sd=10
int.vals <- rnorm(size, mean=coefval)

#print( int.vals )
#print( paste("Coefficient", coefval, sep="; ") )

## collect data to estimate suitability of data
datafit <- data.frame( int=int.vals, ps.new=c(rep(NA,length(int))), w.new=c(rep(NA,length(int))) )

for(i in 1:length(int.vals)){
## generate new data at given coef and intercept
newdat <- simulate.con(object, use.data, newcoef=TRUE,coefname="hvals", coefval,newint=TRUE, intval=int.vals[i], type="response",newTrials, trialsvec)
#print(head(newdat))
datafit$ps.new[i] <- sum( newdat[,1] )/sum( newdat[,2] )
## weight for observation at this row is the total number of trials across the observations
datafit$w.new[i] <- sum( newdat[,2] )
}


#datafit1 <<- datafit

## fit logistic regression model
## this uses the MASS glm function
## MASS is the basis of glmmADMB so this should not be too deviated
## this will return a warning that can be ignored
## weights - total number of trials included
fit <- glm( ps.new ~ int.vals, data=datafit, family=binomial(logit), weights=w.new )

modfit <<- fit

## use this to identify the intercept value for the original proportion successes
opt.int <- (ps.orig/(1-ps.orig)) - summary(fit)$coefficients[1,1]/summary(fit)$coefficients[2,1]

return( list( datafit, fit, opt.int ) )

}

##################################################################################
### simulateCon ##
## function to generate new response data at a given coefficient value
## newcoef is a boolean, stating whether the coefficient value should be altered
## if TRUE then coefval gives the value of the coefficient to use
## coefname is the name of the coefficient of interest, must match the name given in the data
## newint is  a boolean, stating whether to estimate a better fitting intercept value
## type gives the intended output of data - if "response" data is returned in original format
simulate.con <- function(object, use.data, newcoef=FALSE,coefval=0,coefname="",newint=FALSE, intval=0, type="response", newTrials=FALSE, trialsvec=0){
#print("1")
## names of fixed terms
fixed <- attr( attr(object$frame, "terms"), "term.labels" )
random <- attr( object$q, "names" )
## variance can be extracted, convert to sd
random.sigma <- sqrt( as.numeric(object$S[[1]]) )
## random effect values using ranef function
#random.vals <- ranef(object)[[1]]
## random effect values using random normal dist of sd random.sigma
random.vals <- matrix( ncol=1, rnorm( nrow(object$U[[1]]), sd=random.sigma ) )
rownames(random.vals) <- rownames(object$U[[1]])
coefs <- summary(object)$coefficients[,1]
## extract intercept value
int <- as.numeric(object$b[names(object$b)=="(Intercept)"])

## can't see how to pull vector of random effect groups for each row from the data
## have to take data input for now
## if new data is passed use this as the data frame
## else use the data stored in the model object
#if( newdata ){
#	dat.f <- newdata
#}else{
#	dat.f <- object$frame
#}

## New data frame of variables altered by coefficient vals
moddf <- matrix( nrow=nrow(use.data), ncol=length(c(fixed, random)) )
colnames(moddf) <- c(fixed, random)
## needs to be a df to support non numeric data classes
moddf <- as.data.frame(moddf)

##data.fr <<- moddf

#print( head( moddf ) )
#print(ncol(moddf))
## add in each column's data
for(i in 1:ncol(moddf)){
#print(i)
## if describing an interaction term
if( length(grep(":", colnames(moddf)[i]))>0 ){
	## get the involved variables by splitting the term name
	vars <- strsplit( colnames(moddf)[i], ":" )[[1]]
	#print(vars)
	## assume only 2 factors interacting for now
	## would have to find a way to generalise
	## return data as the two sets of data multipled
	moddf[,i] <- use.data[,colnames(use.data)==vars[1]] * use.data[,colnames(use.data)==vars[2]]
}else{
	## else return the matching data from the data frame inputted
	moddf[,i] <- use.data[,colnames(use.data)==colnames(moddf)[i]]
}
}

#print( head( moddf ) )

## add in values for random effect terms
moddf <- cbind(moddf, r.v=c(rep(NA,nrow(moddf))) )
#print(head(moddf))

for(i in 1:length(random.vals)){
moddf$r.v[moddf[,colnames(moddf)==random]==rownames(random.vals)[i]] <- random.vals[i]
}
#print("2")
moddf.adj <- moddf

## if there is to be a new value for a coefficient
if( newcoef==TRUE ){
	for(i in 1:ncol(moddf.adj)){
		## if this is the predictor values for the coefficient of interest, use the new value
		if(colnames(moddf.adj)[i] %in% fixed & colnames(moddf.adj)[i]==coefname){
			moddf.adj[,i] <- moddf.adj[,i] * coefval
		}else if( colnames(moddf.adj)[i] %in% fixed & colnames(moddf.adj)[i]!=coefname ){
		## else adjust the predictor values by the coefficient from the model object
			moddf.adj[,i] <- moddf.adj[,i] * coefs[names(coefs)==colnames(moddf.adj)[i]]
		}
	}
}else{
	## else generate values from the coefficients as is
	for(i in 1:ncol(moddf.adj)){
		if( colnames(moddf.adj)[i] %in% fixed ){
		## adjust the predictor values by the coefficient from the model object
			moddf.adj[,i] <- moddf.adj[,i] * coefs[names(coefs)==colnames(moddf.adj)[i]]
		}
	}
}

## if there is to be a new value for the intercept
if( newint==TRUE ){
int <- intval
}else{
int <- int
}

#datafr <<- moddf.adj
#intvals <<- int

## logit output is then the rowSums of the fixed terms
## plus the value of the random effect
out <- int + rowSums( moddf.adj[,colnames(moddf.adj) %in% fixed ] ) + moddf.adj$r.v


## if linear predictor requested, leave linear predictor as is
if( type=="linear predictor" ){
resp <- out

## if response requested (default), 
}else if( type=="response" ){
## convert linear predictor by antilogit
resp <- (exp(out))/(1+(exp(out)))
respprobs <- resp
#print(head(resp))
## if we have multiple trials (data is n success n trials) convert
	if( ncol(object$frame$response==2) ){
	resp1 <- matrix( nrow=nrow(object$frame$response), ncol=2, NA )
	#print(nrow(resp1))
	if( newTrials==TRUE ){
	#print(length(trialsvec))
	resp1 <- matrix( nrow=length(trialsvec), ncol=2, NA )
	resp1[,2] <- trialsvec
#	print(head(resp1))
	}else{
	resp1[,2] <- object$frame$response[,2]
	}
	#print(nrow(resp1))
	#print(head(resp1[,2]))
	#print(head(resp))
	#print(head(rbinom( n=nrow(resp1), size=resp1[,2], prob=resp )))
	resp1[,1] <- rbinom( n=nrow(resp1), size=resp1[,2], prob=resp )
	resp <- resp1
	}

}

return(resp)

}

#################################################################
### multiSim ###
## to allow simulateCon to be run over a range of coefficient values
## find p vals at each coeff for multiple generations of data
## nsim is number of simulations to get at each coefficient size
multiSim <- function(coefvec, object, resp.data, size, use.data, nsim, newTrials, trialsvec, newint=TRUE, intval=0, newcoef=TRUE,
	coefname="hvals"){
	
for(i in 1:length(coefvec) ){
#print(coefvec[i])
#print( paste("coef at position", i, sep=": "))
## find the best fitting intercept at each coefficient value
if( newint==TRUE ){
newint <- fitInt( object, resp.data, size, coefval=coefvec[i], use.data, newTrials, trialsvec )[[3]]
}

	for(j in 1:nsim){
		## generate new response data at current coefficient with best fit intercept
		response <<- simulate.con(object, use.data, newcoef,coefname, coefval=coefvec[i],newint, intval=newint, type="response",
			newTrials, trialsvec)
	#	print(nrow(response))
		#print(head(response))
		## then re fit the model
		## TODO find way to pass the function as an argument
		fit <- glmmadmb( response ~ (hvals) + (sec.occ) + (1|sector),
			family="binomial", data=use.data )
			
		#print(summary(fit))
		## run drop 1 to get p value for removing h vals term
		comp <- drop1(fit, test="Chisq")
		#print(comp)
		#if( useGLM==TRUE ){
			## fit for glm
			## new glm response
		#	response.glm <- simulate.nm(object.glm, use.data, newcoef=TRUE,coefval=coefvec[i], changecoef="hvals", newint=FALSE, intval=0, type="response")
			## re fit glm
		#	fit.glm <- glm( response.glm ~ hvals + sec.occ, family="binomial", data=use.data )
			## run drop 1 to get p value for removing h vals term
		#	comp.glm <- drop1(fit.glm, test="Chisq")
			## again find a way to reference the correct p value rather than giving a direct location
			## "coef", "run", "int", "pval", "stde", "glmp"
			#out <<- rbind( out, c(coefvec[i], j, newint, comp$Pr[rownames(comp)=="hvals"], 
			#	summary(fit)$coefficients[2,2], comp.glm$Pr[rownames(comp.glm)=="hvals"]) )
		#}else{
		## force it to bind outside of the function environment for now
		out <<- rbind( out, c(coefvec[i], j, newint, comp$Pr[rownames(comp)=="hvals"], 
				summary(fit)$coefficients[2,2]) )
		#}		
		#if( comp$Pr[rownames(comp)=="hvals"] == 1 ){
			new <- cbind(response, cbind(c(rep(coefvec[i],nrow(response))),c(rep(j,nrow(response))) ),use.data[,colnames(use.data)=="hvals"] )
			dataset <<- rbind( dataset, new )
			simresp <<- rbind( simresp, response )
		#}
	}

}

return(out[2:nrow(out),])


}