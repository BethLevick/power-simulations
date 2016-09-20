##################################################################################################
## Set up Kazakhstan data sets for GLMM power analysis methods ##
##################################################################################################
is.installed(pkg="glmmADMB")

## load in data set: df called "tab"
load("events-table.RData")
## load in distance matrix for calculating h values: matrix called "dmat"
load("distMat-km.RData")
## transitions data set for expanding sample sizew
load("transitionsData.RData")
## factors data for larger sample
load("factors.RData")

## save and re load instead of re doing
load( "trials-vec.RData" )
load( "input-new.RData" )

## split this into colonisation response, extinction response, random effect
## random effect is 3rd col
secs <- as.factor(tab[,3])
## occupancy is first 2 columns (extinction data is in this data frame as well, but we ignore that here)
col.r <- data.matrix( tab[,1:2] )
##########################################
## adjust by average sector occupancy
## add sector info to transitions data frame

## build h values
coltab <- data.frame( sector=secs, hvals=as.numeric(getH(0.11, dmat)), prop.suc=col.r[,1]/col.r[,2] )

trans.s <- cbind( as.data.frame(trans), sector=as.character(factor_df$sector) )

sec.occ <- data.frame( sector=unique(trans.s$sector), av.occ=c(rep(NA,length(unique(trans.s$sector)))) )

for(i in 1:nrow(sec.occ)){
tmp <- trans.s[trans.s$sector==sec.occ$sector[i],]

#print( head(tmp) )

sec.occ$av.occ[i] <- mean( (rowSums(tmp[,1:7], na.rm=T)/nrow(tmp)), na.rm=T)
}

## then send these calculated averages back to factor_df
## for each sector
factor_df$sec.occ <- c(rep(NA,nrow(factor_df)))
for(i in 1:nrow(sec.occ)){
factor_df$sec.occ[factor_df$sector==sec.occ$sector[i]] <- sec.occ$av.occ[i]
}


## add onto model dfs
coltab$sec.occ <- factor_df$sec.occ
###################################################################
## re name columns in the new data to allow match across
colnames(newtab)[3] <- "sector" 
newtab$sector <- as.factor(newtab$sector)
###################################################################