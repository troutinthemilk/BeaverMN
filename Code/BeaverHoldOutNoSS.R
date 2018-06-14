##Jake Ferguson 
#runs a gompertz population dynamic model on a number of beaver populations in JAGS

library(R2jags)
library(dplyr)
library(RColorBrewer)
library(readxl)
library(coda)
#library(superdiag)
#library(mcmcplots)

curr.hyp		<- 'notall' #hypotheses as denoted by sean
SS 				<- F
PDSI 			<- T
holdout.prop	<- 0.0#0.25

#1. Dispersal @2y/o; winter affects reproduction; rainfall/May temp affects reproduction.
#Tmp.wint (lag 2); tmax.5, ppt.grow, ppt.fall (all lag 3); harvest (lag 1); wolf presence (no lag); habitat, ppt.grow (lag 0)
#2. Dispersal @2y/o; winter affects juvenile survival; rainfall/May temp affects reproduction.
#Tmp.wint (no lag); tmax.5, ppt.grow, ppt.fall (all lag 3); harvest (lag 1); wolf presence (no lag); habitat, ppt.grow (lag 0)
#3. Dispersal @2y/o; winter affects reproduction; rainfall/May temp affects kit body condition.
#Tmp.wint, tmax.5, ppt.grow, ppt.fall (all lag 2); harvest (lag 1); wolf presence (no lag); habitat, ppt.grow (lag 0)
#4. Dispersal @2y/o; winter affects juvenile survival; rainfall/May temp affects kit body condition.
#Tmp.wint (no lag); tmax.5, ppt.grow, ppt.fall (all lag 2); harvest (lag 1); wolf presence (no lag); habitat, ppt.grow (lag 0)

##read in data

colony.dat 		<- as.data.frame(read_excel("../Data/Annual_census_5.31.18.xlsx", sheet=1, na="NA"))
#colony.dat 		<- as.data.frame(read_excel("../Data/Annual_census_4.9.18.xlsx", sheet=1, na="NA"))
envvars.dat 	<- as.data.frame(read_excel("../Data/Annual_census_5.31.18.xlsx", sheet=2, na="NA"))[,-c(1,4)]
habvars.dat 	<- as.data.frame(read_excel("../Data/Annual_census_5.31.18.xlsx", sheet=3, na="NA"))

routename.vec   <- c("Blackduck", "Cass",  "Cass_crow", "C_st_louis", "Ely_finger", "Hay_kelliher", "Itasca", "Kabetogama", "Kanabec", "Kawishiwi", "Kooch_north", "Northome", "Red_lake", "Southern_pine", "West_vermillion") #removed S_st_louis

colony.dat 	<- colony.dat %>% filter(rte.name %in% routename.vec) %>% droplevels
envvars.dat <- envvars.dat %>% filter(rte.name %in% routename.vec) %>% droplevels


DroughtSensPCA 	<- prcomp(habvars.dat[,24:29], center=T, scale.=T)
habvars.dat		<- cbind(Route_name=habvars.dat$Route_name, data.frame(Quality=habvars.dat$Decid_pct + habvars.dat$Mixed_pct, veh.dens=habvars.dat$Veh_dens, PC1=predict(DroughtSensPCA)[,1]))

habvars.dat	 	<- habvars.dat[order(habvars.dat$Route_name),]

pc1.vec			<- rep(habvars.dat$PC1, each=31)

if(curr.hyp=="all") 
{	if(PDSI) {
		curr.ind 		<- c(3, 29,29,29, 28,28,28, 44,44,44,44, 45,45,45,45)
	} else {
		curr.ind 		<- c(3, 29,29,29, 28,28,28, 30,30,30,30, 31,31,31,31)
	}	
} else {
	if(PDSI) {
		curr.ind 		<- c(3, 29,29, 28,28, 44,44, 45,45, 46)
	} else {
		curr.ind 		<- c(3, 29,29, 28,28, 30,30, 31,31, 47)
	}

}
lag.vals		<- switch(curr.hyp, "1"=c(1,2,3,3,0,3,0), "2"=c(1,0,3,3,0,3,0), "3"=c(1,2,2,2,0,2,0), "4"=c(1,0,2,2,0,2,0), "all"=c(1, 0,1,2, 1,2,3, 0,1,2,3, 0,1,2,3, 0, 0), "notall"=c(1, 0,2, 2,3, 2,3, 2,3, 0, 0)) 

envvars.dat		<- envvars.dat[, c(1,2,curr.ind)]
envvars.dat		<- cbind(envvars.dat, Quality=rep(habvars.dat$Quality, each=31))
envvars.dat 	<- cbind(envvars.dat[,1:2], scale(envvars.dat[,-(1:2)]))

#need to keep track of how to incude pc1 and vehicle density.
#envvars.dat[,3] <- envvars.dat[,3]*rep(habvars.dat$veh.dens, each=31)
#envvars.dat[,8] <- envvars.dat[,8]*pc1.vec
#envvars.dat[,9] <- envvars.dat[,9]*pc1.vec
#envvars.dat[,10] <- envvars.dat[,10]*pc1.vec
#envvars.dat[,11] <- envvars.dat[,11]*pc1.vec

#envvars.dat		<- cbind(envvars.dat, PC1=pc1.vec, Quality=rep(habvars.dat$Quality, each=31))

year.min	<- min(colony.dat$year)
year.max 	<- max(colony.dat$year)
Dens.mat 	<- matrix(NA, year.max - year.min + 1, length(routename.vec))

dimnames(Dens.mat)[[1]] <- year.min:year.max
dimnames(Dens.mat)[[2]] <- routename.vec

Count.mat 		<- Observer.mat <- ObsType.mat <- km.mat <- Wolf.mat <- Dens.mat
init.index 		<- final.index <- vector('numeric', length(routename.vec))

#construct design matrices for observation error model
options(na.action='na.pass')
ObsX 			<- model.matrix(year ~ obs.type + observer1, data=colony.dat)
observCov.mat 	<- ObsX[,-1]
ObsCov.array	<- array(NA, c(year.max - year.min + 1, dim(ObsX)[2]-1, length(routename.vec))) 
EnvCov.array 	<- array(NA, c(year.max - year.min + 1, dim(envvars.dat)[2]-1, length(routename.vec))) 

ObsX 			<- model.matrix(year ~ obs.type + observer1, data=colony.dat)

dimnames(ObsCov.array)[[1]] <- year.min:year.max
dimnames(ObsCov.array)[[2]] <- dimnames(observCov.mat)[[2]]
dimnames(ObsCov.array)[[3]] <- routename.vec

dimnames(EnvCov.array)[[1]] <- year.min:year.max
dimnames(EnvCov.array)[[2]] <- c(dimnames(envvars.dat)[[2]][-c(1:2)], "Wolf")
dimnames(EnvCov.array)[[3]] <- routename.vec

obsna.mat	 <- densnotna.mat <- matrix(NA, 30, length(routename.vec))#store indices that are NA in observation covariates
obsna.length <- densnotna.length <- holdout.obs <- vector('numeric', length(routename.vec))

##format input datasets for the jags model

for(i in 1:length(routename.vec)) {

	currCol 		<- colony.dat[which(colony.dat$rte.name == routename.vec[i]),]
	currObsCov 		<- observCov.mat[which(colony.dat$rte.name == routename.vec[i]),] 
		
	if(is.na(currCol$col.km[2])) {
		currCol 		<- currCol[-(1:3),]
		currObsCov		<- currObsCov[-(1:3),]
		curr.yearmin 	<- min(currCol$year)
		curr.yearmax 	<- max(currCol$year) 

	} else {
		curr.yearmin 	<- min(currCol$year)
		curr.yearmax 	<- max(currCol$year)
	}	

	init.index[i] 	<- curr.yearmin-year.min + 1
	final.index[i] 	<- curr.yearmax - year.min + 1
	curr.indices 	<- init.index[i]:final.index[i]

	holdout.obs[i] 					<- round((final.index[i]-init.index[i]+1)*holdout.prop)

	Dens.mat[curr.indices, i]		<- currCol$col.km 
	Wolf.mat[curr.indices, i]		<- currCol$wolf_ubi #as.numeric(currCol$wolf == "present")
	Count.mat[curr.indices, i]		<- currCol$num.col
	ObsCov.array[curr.indices,,i]	<- currObsCov
	km.mat[curr.indices, i]			<- currCol$rte.km 

	#format environmental covariates
	envTemp 	<- filter(envvars.dat, rte.name==routename.vec[i]) #filter for current route
	First.index	<- which(envTemp$year == curr.yearmin) - lag.vals
	Last.index 	<- which(envTemp$year == curr.yearmax) - lag.vals
	envTemp		<- envTemp[,-c(1,2)] #remove year and routename from covariate matrix

	temp <- matrix(NA, length(First.index[1]:Last.index[1]), dim(envTemp)[2])
	for(j in 1:dim(temp)[2]) {
		envVec	<- envTemp[First.index[j]:Last.index[j],j]
		if(sd(envVec) == 0){
			temp[,j] <- envVec
			next()
		} #one case is all zeros so can't use scale function
		temp[,j] <- envVec#(envVec - mean(envVec))/sd(envVec)
	}

	EnvCov.array[curr.indices,1:dim(temp)[2],i]	<-  temp #cast to matrix, scale covariates

	#missing covariate observations
	na.vec <- union(which(is.na(observCov.mat[which(colony.dat$rte.name == routename.vec[i]),1])), which(is.na(observCov.mat[which(colony.dat$rte.name == routename.vec[i]),2])))

	obsna.mat[1:length(na.vec), i] <- na.vec
	obsna.length[i] <- length(na.vec)

	#missing density observations
	notna.vec <- which(!is.na(Dens.mat[curr.indices, i]))
	densnotna.mat[1:length(notna.vec), i] <- notna.vec + init.index[i] - 1
	densnotna.length[i] <- length(notna.vec)

}

#now add in wolf covariate
EnvCov.array[,dim(EnvCov.array)[2],] <- Wolf.mat

size.mat <- observCov.mat
size.mat[which(observCov.mat==0)] <- 1
obsfreq <- apply(observCov.mat, 2, sum, na.rm=T)/apply(size.mat, 2, sum, na.rm=T)

for(i in 1:dim(EnvCov.array)[2]) {
	dimnames(EnvCov.array)[[2]][i] <- paste(dimnames(EnvCov.array)[[2]][i], '-lag', lag.vals[i], sep='')
}


#HarCov.array <- EnvCov.array[,1,]
#EnvCov.array <- EnvCov.array[,-1,]
#state space model on the log scale. 
##Includes habitat covariates as influence on intrinsic growth rate (betaHab)
##Includes observer covariates as influence on the observed popualtion state (betaObs)

NumPops		<- length(routename.vec)
numEnvVars	<- dim(EnvCov.array)[2]
datTable 	<- list(lD=log(Dens.mat), D=Dens.mat, NumPops=NumPops, init.index=init.index, final.index=final.index, ObsCov.array=ObsCov.array, NumObsVars=dim(ObsCov.array)[2], obsna.mat=obsna.mat, obsna.length=obsna.length, EnvCov.array=EnvCov.array, NumEnvVars=numEnvVars, HarCov.array=NULL, holdout.obs=holdout.obs)
overall.dens 	<- apply(Count.mat, 1, sum, na.rm=T)/apply(km.mat, 1, sum, na.rm=T)
l.dens 			<- log(overall.dens)

#overall.fit <- lm(l.dens[-1] ~ l.dens[-length(l.dens)])
initTable 	<- function() {
	list(a=rnorm(1, -0.1, 0.1), b=rnorm(1, -0.1, 0.01))
}

#stop()

if(SS) {
	par.est	<- c('a', 'sigmaInt', 'b', 'b2', 'b3', 'sigmaSlope', 'sigObs', 'lambdaEnv', 'lambdaHar', 'lambdaObs', 'betaEnv', 'betaHar', 'Xtrue', 'Xpred', 'aRE')
	jags.fit  		<- jags(model.file='BeaverSSmodels.jags', n.iter=1e4, n.chains=4, data=datTable, inits=initTable, parameters.to.save=par.est)
} else {
	par.est	<- c('a', 'sigmaInt', 'b', 'b2', 'sigmaSlope', 'sigInt', 'sigObs', 'lambdaEnv', 'lambdaHar', 'lambdaObs', 'betaEnv', 'betaHar', 'Xpred', 'aRE')
	jags.fit  		<- jags(model.file='Beavermodels.jags', n.iter=1e4, n.chains=4, data=datTable, inits=initTable, parameters.to.save=par.est)
}

jags.fit 		<- update(jags.fit, n.iter=1e4, n.thin=1e2)
#dic.samples(jagsfit.upd$model, n.iter=1e4, type="pD")
#jags.fit 		<- recompile(jags.fit)
#jags.fit 		<- autojags(jags.fit, n.iter=1e4, n.thin=100, n.update=100)

cat('current hypotheses', curr.hyp, '\n')
cat("DIC", jags.fit$BUGSoutput$DIC, "\n")

jags.fit.mcmc 	<- as.mcmc(jags.fit)

#acfdiag <- autocorr.diag(jags.fit.mcmc, lags=c(1))

#densityplot(jags.fit.mcmc)
##model output figure code##
jags.fit.summ 	<- summary(jags.fit.mcmc)[[1]]
names.vec 		<- dimnames(jags.fit.summ)[[1]]
index 			<- NULL
Xpred	 		<- Xtrue <- Dens.mat

##extract and order Xtrue and Xpred values##
##Xtrue is the true underlying (log)population size
##Xpred is the model prediction from the Gompertz model
for(i in 1:length(names.vec)) {
	if(length(grep("Xtrue", names.vec[i]))>0) {
		index <- c(index, i)
		temp <- as.numeric(unlist(regmatches(names.vec[i], gregexpr("[[:digit:]]+", names.vec[i]))))
		Xtrue[temp[1], temp[2]] <- jags.fit.summ[i,1]
	}
	if(length(grep("Xpred", names.vec[i]))>0) {
		index <- c(index, i)
		temp <- as.numeric(unlist(regmatches(names.vec[i], gregexpr("[[:digit:]]+", names.vec[i]))))
		Xpred[temp[1], temp[2]] <- jags.fit.summ[i,1]
	}

}

#plot stuff
#graphics.off()
X11()
par(mfrow=c(4,4)) 
Year.vec <- as.numeric(dimnames(Dens.mat)[[1]])
for(i in 1:NumPops) {
	plot(Year.vec, Dens.mat[,i], type='b', pch=19, cex=0.75, xlab="Year", ylab="Density", cex.lab=1.4, main=routename.vec[i], col="darkgray")
	lines(Year.vec, exp(Xtrue[,i]), lwd=2)
	lines(Year.vec, exp(Xpred[,i]), lwd=2, col="cornflowerblue")
}

X11()
par(mfrow=c(4,4)) 
resid.mat <- Xpred - Xtrue
for(i in 1:NumPops) {
	plot(Year.vec, resid.mat[,i], type='l', cex=0.75, xlab="Year", ylab="Density", cex.lab=1.4, main=routename.vec[i], col="black")
}


X11()
#par(mfrow=c(1,2))
library(MCMCvis)
MCMCplot(jags.fit.mcmc, params="betaEnv", main="Environmental covariates", ref_ovl=TRUE)

#MCMCplot(jags.fit.mcmc, params="betaHar", labels="Harvest:veh_dens", main="Interaction covariates")
#MCMCplot(jags.fit.mcmc, params="betaObs", labels=dimnames(observCov.mat)[[2]], main="Observation variance covariates")

X11()
MCMCtrace(jags.fit.mcmc, params = c('lambda'), ISB=FALSE, ind=TRUE, pdf=FALSE)
X11()
MCMCtrace(jags.fit.mcmc, params = c('a', 'b', 'b2'), ISB=TRUE, ind=TRUE, pdf=FALSE)

#stop()
if(SS) {
	if(PDSI) {
		#next()
		save.image(file=paste("BeaverGompSS_QuarterHoldout_PDSI_Holdout", holdout.prop, ".Rdata", sep=''))
	} else {
		#save.image(file=paste("BeaverGompSS_QuarterHoldout_PPT_", curr.hyp, ".Rdata", sep=''))
	}
} else {
	#save.image(file=paste("BeaverGomp_Hyp", curr.hyp, ".Rdata", sep=''))
}

#plot(Xpred[holdout.indices,i], Xpred[holdout.indices,i] - Xpred.nocov[holdout.indices,i], pch=19); points(Xpred.best[holdout.indices,i], Xpred.best[holdout.indices,i] - Xpred[holdout.indices,i], pch=19, col='red'); abline(a=0, b=0)

#mse.nocov <- mse.cov <- vector('numeric', 15)
#for(i in 1:15) {
#	mse.nocov[i] <- sum(abs(Xpred.best[holdout.indices,i] - Xpred.nocov[holdout.indices,i]), na.rm=T)
#	mse.cov[i] <- sum(abs(Xpred.best[holdout.indices,i] - Xpred[holdout.indices,i]), na.rm=T)
#}

a <- jags.fit.summ[which(names.vec=="a"),1]
b <- jags.fit.summ[which(names.vec=="b"),1]
aRE.vec  <- vector('numeric', 15)
beta.vec <- vector('numeric', 44)
for(i in 1:(dim(envvars.dat)[[2]]-3)) {
	beta.vec[i] <- jags.fit.summ[which(names.vec==paste("betaEnv[",i,"]",sep='')), 1]
}

for(i in 1:15) {
	aRE.vec[i] <- jags.fit.summ[which(names.vec==paste("aRE[",i,"]",sep='')), 1]
}


X11()
par(mfrow=c(3,5))
mse.vec <- vector('numeric', 15)

for(i in 1:15) {

	holdout.indices <- (final.index[i]-holdout.obs[i]):final.index[i]
	pred.start <- (final.index[i]-holdout.obs[i])-1

	X <- log(Dens.mat[,i])

	plot(X[holdout.indices], type='b', ylim=c(-2,0.5), main=i, pch=19)

	pred.val <- vector('numeric', length(holdout.indices))
	for(j in 1:length(holdout.indices)) {
		pred.val[j] <- X[holdout.indices[j]-1] + a + aRE.vec[i] + b*X[holdout.indices[j]-1] + sum(beta.vec*EnvCov.array[holdout.indices[j],,i])
		
		if(is.na(pred.val[j])) {
			if(j == 1) {next()}
			pred.val[j] <- pred.val[j-1] + a + aRE.vec[i] + b*pred.val[j-1] + sum(beta.vec*EnvCov.array[holdout.indices[j],,i])	
		}
	}

	lines(pred.val, lwd=2, col='cornflowerblue')

	mse.vec[i] <- mean(abs(X[holdout.indices][-1] - pred.val[-1]), na.rm=T)
}


#no state-space model
#no cov mean(mse.vec) = 0.2424127
#lag 0 mean(mse.vec) = 0.2251269
#lag 1 mean(mse.vec) =  0.2252731
#lag 2 mean(mse.vec) = 0.2121801
#lag 1,2 mean(mse.vec) = 0.2153151
#lag 0,2 mean(mse.vec) = 
#lag 2 mean(mse.vec), lasso = 0.2296801
#lag 2 mean(mse.vec), ridge = 
#PDSI mean(mse.vec) = 0.2076699, -13.83146
#PPT mean(mse.vec) = 0.231497, -2.482535