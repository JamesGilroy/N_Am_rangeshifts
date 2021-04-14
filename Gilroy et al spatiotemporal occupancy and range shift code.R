library(nlme)#,lib.loc = 'F:/R-3.4.0/library')
library(mgcv)#,lib.loc = 'F:/R-3.4.0/library')
require(rjags)#,lib.loc = 'E:/R-3.1.1/library')
require(R2jags)#,lib.loc = 'E:/R-3.1.1/library')
library(igraph)#,lib.loc = 'E:/R-3.1.1/library')

# Code to fit hierarchical spatiotemporal occupancy model, generate mean occupancy rasters for the time period, and generate posterior samples of estimated range limit shifts 
# Example data is provided for CBC surveys for Setophaga pinus
# Species observation and site data objects were a priori filtered to only include locations with coverage in eeach 3-year window across the study period (see Methods)

survey <- "CBC"

years <- 1986:2015
wind <- 3 #length of  occupancy closure window
windows <- seq(min(years),max(years),wind)
nwindows <- length(windows)
findInterval(years,windows)


#####################################################
## Rread in data objects
# get all survey points
  covconst <- read.csv(paste(survey,"_covermatrix.csv",sep=""),row.names="X") #annual coverage of all sites
  locs <- read.csv(paste(survey,"_locations.csv",sep=""),row.names="X") #all site locations
  windconst<-read.csv(paste(survey,"_wind.csv",sep=""),row.names="X") #survey wind conditions at retained sites
  cloudconst <- read.csv(paste(survey,"_cloud.csv",sep=""),row.names="X")#survey cloud conditions
if(survey=="CBC"){
	obsconst <-read.csv(paste(survey,"_obsmatrix.csv",sep=""),row.names="X")
}
 #restrict data to the years included in study 
colpick <- which(names(covconst)==paste("y",min(years),sep="")): which(names(covconst)==paste("y",max(years),sep=""))
covconst <- covconst[,colpick]
windconst <- windconst[,colpick]
cloudconst <- cloudconst[,colpick]
if(survey=="CBC"){  obsconst <- obsconst[,colpick]}

n <- nrow(cloudconst) #retained sites
n2 <- nrow(covconst) #all sites
# Build annual dataset:
covpts <- locs[1:n,] #all locations retained for modelling
covpts$id1 <- paste(covpts$lon,"_",covpts$lat,sep="")
datx1 <- data.frame(lon=locs$lon[1:n],lat=locs$lat[1:n],y=1) # Put into data frame
datx1 <- rbind(datx1,data.frame(lon=locs$lon[(n+1):n2],lat=locs$lat[(n+1):n2],y=1)) # Put into data frame
datx<-data.frame()
for(i in 1:nwindows){datx <- rbind(datx,data.frame(datx1,window=i,point=1:n2))} #add replicate spatial data for each window

#standardize predictors
datx$window2 <- (datx$window-mean(datx$window))/sd(datx$window)
datx$lon <- (datx$lon-mean(datx$lon))/sd(datx$lon)
datx$lat <- (datx$lat-mean(datx$lat))/sd(datx$lat)

#Now read in species data
dat <- read.csv(paste(survey,"_","Setophaga_pinus.csv",sep="")) #get species data
dat$id1 <- paste(dat$long,dat$lat,sep="_") #location ID to match to sampling data
k1 <- c(35,18) #long/lat knots
k2 <- c(floor(k1/2), nwindows-1) #temporal interaction knots
#Prepare spline basis for spatiotemporal site occupancy:
jags.ready <- jagam(y~ te(lon,lat,bs="tp",k=k1) +ti(lon,lat,window2,bs="tp",k=k2),	data=datx, family = "binomial",file="jagam_simpleF.bug")
S1<-jags.ready$jags.data$S1
S2<-jags.ready$jags.data$S2
X <- jags.ready$jags.data$X
zero <- jags.ready$jags.data$zero

#remove background /non-retained sites from object for JAGS 
X <- X[datx$point<=n,]
datx <- datx[datx$point<=n,]

# Now match  sampling events to occupancy matrix 
y <- key <- wind <- cloud <-obs <- c()
zobs <- rep(0,nrow(datx)) 
for(i in 1:n){ #loop through sites
  evs <- which(covconst[i,]==T) #how many sampling events 
  counts <-vapply(evs,FUN=function(x) sum(dat$counts[dat$id1==covpts$id1[i] & dat$Year==years[x]],na.rm=T),FUN.VALUE = 1)  # number of individuals observed per event
  windz<-findInterval(years[evs],windows) #time windows for each event
  evs <- evs[windz<=nwindows] 
  counts <- counts[windz<=nwindows]
 	windz <- windz[windz<=nwindows]
 	y<- c(y,counts)
  key <- c(key,vapply(windz,FUN=function(x) which(datx$point==i & datx$window==x),FUN.VALUE = 1)) #identifies occupancy site/window for each event
  zobs[vapply(windz[counts>0],FUN=function(x) which(datx$point==i & datx$window==x),FUN.VALUE = 1)]<-1 #identifies observed presences
  cloud <- c(cloud,as.numeric(cloudconst[i,evs]))  #conditions for events
  wind <- c(wind,as.numeric(windconst[i,evs]))
  if(survey=="CBC"){obs <- c(obs,as.numeric(obsconst[i,evs]))}
}
y[y>1]<-1 #convert to binary data

#standardise predictors
wind <- (wind-mean(wind))/sd(wind)
cloud <- (cloud-mean(cloud))/sd(cloud)
if(survey=="CBC"){obs<-  (obs-mean(obs))/sd(obs)}


#Now prepare the data objects that will be sent to the JAGS model:

kj1<- k1[1]*k1[2] # square of k1 - used in JAGS code for priors matrix dimensions
kj2<- length(zero)-kj1 # #also used for priors matrix dimensions


n2a <- n2
n2 <- length(zobs)
n3 <- length(y)
n2s <- 1:n2
n2s <- n2s[n2s%in%key]
zobs <- zobs[1:max(n2s)] #used to initiate model
zobs[which((1:max(n2s))%in%n2s==F)]<-NA 

# Specify JAGS models for BBS and CBC:
setwd(data)
sink("BBS_simple_SPATEMP.txt")			 		
cat("
    model {
    #hyper-priors for mean presence spline
		eta ~dnorm(0,0.001) #random intercepts
	  for (i in 1:7) {
    lambda[i] ~ dgamma(.05,.005)
	  }
    #latlong spline prior covariance matrices:
    K1 <- S1[1:(kj1-1),1:(kj1-1)] * lambda[1]  + S1[1:(kj1-1),kj1:(2*(kj1-1))] * lambda[2] + S1[1:(kj1-1),(1+2*(kj1-1)):(3*(kj1-1))] * lambda[3] 
    #temporal interaction prior covariance matrices:
  	K2 <- S2[1:(kj2),1:(kj2)] * lambda[4]  + S2[1:(kj2),(kj2+1):(2*(kj2))] * lambda[5] + S2[1:(kj2),(2*kj2+1):(3*(kj2))] * lambda[6] + S2[1:(kj2),(3*kj2+1):(4*(kj2))] * lambda[7] 
    # Priors on mean presence spline
    for(i in 1:1){b[i] ~ dnorm(0,0.0092)} 
    b[2:kj1] ~ dmnorm(zero[2:kj1],K1) 
	b[(kj1+1):(kj1+kj2)] ~ dmnorm(zero[(kj1+1):(kj1+kj2)],K2) 
    	
		beta2 ~dnorm(0,0.001) #linear effect priors detection
  	beta3 ~dnorm(0,0.001) #
  	beta4 <-0

	for(j in n2s){
 	 logit(mu[j]) <-  X[j,] %*% b ## linear predictor
    z[j] ~ dbern(mu[j]) # latent true presence vector
}
 		for (i in 1:n3) { #observed data
  	y[i] ~ dbern(p[i]*z[key[i]]) ## response
    logit(p[i]) <- eta + beta2*wind[i] + beta3*cloud[i]

# Bayesian P
		ynew[i] ~ dbern(eta2[i]) ## sim response
    d[i]<-  (y[i] - eta2[i])/sqrt((eta2[i]+0.001)*(1-eta2[i]-0.001))
    dnew[i]<- (ynew[i] - eta2[i])/sqrt((eta2[i]+0.001)*(1-eta2[i]-0.001))
    d2[i]<- pow(d[i],2)
    dnew2[i]<- pow(dnew[i],2)
}
p.fit<-sum(d2[1:n3])
p.fitnew<-sum(dnew2[1:n3])
    }
    ",fill=TRUE)
sink()

setwd(data)
sink("CBC_simple_SPATEMP.txt")			 		
cat("
    model {
    #hyper-priors for mean presence spline
		eta ~dnorm(0,0.001) #random intercepts
	  for (i in 1:7) {
    lambda[i] ~ dgamma(.05,.005)
    }
    K1 <- S1[1:(kj1-1),1:(kj1-1)] * lambda[1]  + S1[1:(kj1-1),kj1:(2*(kj1-1))] * lambda[2] + S1[1:(kj1-1),(1+2*(kj1-1)):(3*(kj1-1))] * lambda[3] 
  	K2 <- S2[1:(kj2),1:(kj2)] * lambda[4]  + S2[1:(kj2),(kj2+1):(2*(kj2))] * lambda[5] + S2[1:(kj2),(2*kj2+1):(3*(kj2))] * lambda[6] + S2[1:(kj2),(3*kj2+1):(4*(kj2))] * lambda[7] 
    # Priors on mean presence spline
    for(i in 1:1){b[i] ~ dnorm(0,0.0092)} 
    b[2:kj1] ~ dmnorm(zero[2:kj1],K1) 
	b[(kj1+1):(kj1+kj2)] ~ dmnorm(zero[(kj1+1):(kj1+kj2)],K2) 
   	
 		beta2 ~dnorm(0,0.001) #
  	beta3 ~dnorm(0,0.001) #
   	beta4 ~dnorm(0,0.001) #

	for(j in n2s){
	 	logit(mu[j]) <-  X[j,] %*% b ## linear predictor
    z[j] ~ dbern(mu[j]) # latent true presence vector
		}
 
 		for (i in 1:n3) { 
 		y[i] ~ dbern(p[i]*z[key[i]]) ## response
    logit(p[i]) <- eta + beta2*wind[i] + beta3*cloud[i] +beta4*obs[i]#
# Bayesian P
#		ynew[i] ~ dbern(eta2[i]) ## sim response
#    d[i]<-  (y[i] - eta2[i])/sqrt((eta2[i]+0.001)*(1-eta2[i]-0.001))
#    dnew[i]<- (ynew[i] - eta2[i])/sqrt((eta2[i]+0.001)*(1-eta2[i]-0.001))
#    d2[i]<- pow(d[i],2)
#    dnew2[i]<- pow(dnew[i],2)
}
#p.fit<-sum(d2[1:n3])
#p.fitnew<-sum(dnew2[1:n3])
    }
    
     ",fill=TRUE)
sink()
#Initial values:
b.init <- jags.ready$jags.ini$b
lamb.init <- jags.ready$jags.ini$lambda
inits <- function(){ list(
  z=zobs,
   lambda=lamb.init,
 	eta=rnorm(1,0,0.1),
	b=b.init
)}

n.iter <- 10000
n.burn <-5000
n.thin <- 5
n.chains<-3

if(survey=="BBS"){
params <-c("b","lambda","eta","beta2","beta3","p.fit","p.fitnew")
}else{
params <-c("b","lambda","eta","beta2","beta3","beta4","p.fit","p.fitnew")#
}
# Now fit the model for the species
#- WARNING THIS CAN TAKE SEVERAL WEEKS!!!
if(survey=="BBS"){
  model.fit <- do.call(jags, list(c("y","n2s","n3","X","S1","S2","zero","kj1","kj2","key","wind","cloud"), inits, params,                          
      "BBS_simple_SPATEMP.txt",    n.chains, n.iter, n.burn, n.thin))
}else{
 model.fit <- do.call(jags, list(c("y","n2s","n3","X","S1","S2","zero","kj1","kj2","key","wind","cloud","obs"), inits, params,                          
      "CBC_simple_SPATEMP.txt",    n.chains, n.iter, n.burn, n.thin))
}
#Extract ket outputs
b.samps <- model.fit$BUGSoutput$sims.list$b
eta <- model.fit$BUGSoutput$sims.list$eta
beta2 <- model.fit$BUGSoutput$sims.list$beta2
beta3 <- model.fit$BUGSoutput$sims.list$beta3
lambda <- model.fit$BUGSoutput$sims.list$lambda
if(survey=="CBC"){beta4 <- model.fit$BUGSoutput$sims.list$beta4}else{beta4<-NA}

p.fit <- model.fit$BUGSoutput$sims.list$p.fit
p.fitnew <- model.fit$BUGSoutput$sims.list$p.fitnew
bayesp <- data.frame(p.fit,p.fitnew)


b.samps <- data.frame((b.samps))
names(b.samps) <- paste("b",1:ncol(b.samps),sep="")
det.samps <- data.frame(eta,beta2,beta3,beta4)
lambda.samps <- data.frame(lambda)
outsum <- data.frame(model.fit$BUGSoutput$summary)
#Save model outputs
setwd(output)
write.csv(b.samps,paste(survey,"_",namlist$name[p],"_UPDATEbsamps3.csv",sep=""))
write.csv(det.samps,paste(survey,"_",namlist$name[p],"_UPDATEdetsamps3.csv",sep=""))
write.csv(outsum,paste(survey,"_",namlist$name[p],"_UPDATEsummary3.csv",sep=""))
write.csv(lambda.samps,paste(survey,"_",namlist$name[p],"_lUPDATEambdasamps3.csv",sep=""))
write.csv(bayesp,paste(survey,"_",namlist$name[p],"_UPDATEbayesp3.csv",sep=""))

	
### Now generate  occupancy rasters:
covzone <- raster(paste("Cover_area.grd",sep=""))#study area mao

# occupancy rasters:
occp <- 1/(1+exp(-(X%*%colMeans(b.samps)))) #calculates mean modelled occupancy 
# Create rasters for each window:
ras <- covzone
for(t in 1:nwindows){
r <- data.frame(lon=datx$lon[datx$window==t],lat=datx$lat[datx$window==t],prob=occp[datx$window==t])
r <- subset(r,!is.na(r$lat))
names(r)<-c("lon","lat","prob")
coordinates(r)<-c("lon","lat")
r <- rasterize(r,covzone,field=r$prob,fun=mean)
ras <- stack(ras,r)
}
ras <- ras[[-1]]
	setwd(rasloc)
	writeRaster(ras1,paste(survey,"_",namlist$species[p],octhresh,"_occMean.grd",sep=""),overwrite=T)
##################################################################################	
	
# RANGE-SHIFT ESTIMATION ALGORITHM
	
##################################
setwd(home)
covzone2 <- aggregate(covzone,fact=2) 
f1 <- boundaries(covzone2,directions=4) # find study area edges

# Step 1:
# Find mean range boundaries and calculate their relative orientation
midlev <- 1
sr1 <- proj4string(covzone)
sr <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-110+x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m"
covzone<-projectRaster(covzone,crs=sr)
covzone[!is.na(covzone)]<-1
ras<-projectRaster(ras,crs=sr)
oras1 <-alras<- ras[[midlev]] #create objects for edge data
oras1[oras1<0.1]<-NA #find range patches
oras1 <- clump(oras1) # classify range patches
ncs <- cellStats(oras1,max) 
eds1 <- data.frame(matrix(nrow=0,ncol=3)) #object to store edge cells
names(eds1)<-c("x","y","clump")
octhresh <- 0.1 #occupancy threshold
for(i in 1:ncs){##loop through range patches
  if(sum(getValues(oras1)==i,na.rm=T)>=10){
    cl1<- ras[[midlev]]
    cl1[oras1!=i]<-0
    cl1[is.na(cl1)]<-0
    if(cellStats(cl1,max)<=octhresh){octhresh1<-0.9*cellStats(cl1,max)}else{octhresh1<-octhresh}
    rls <- rasterToContour(cl1,levels=octhresh) #apply threshold
    #lines(rls)
    rls <- as(rls,"SpatialPointsDataFrame")
    rls <- spTransform(rls,CRS(sr1))
    #extract range edge cells
    eds <- data.frame(rls@coords)
    names(eds)<-c("x","y")
    eds$clump<-i
  #  keepe <- rep(c(1:2),times=nrow(eds)+1)[1:nrow(eds)]
  #  eds<-eds[keepe==1,]
    eds1<-rbind(eds1,eds)
  }
}

#eds1 in orig projection latlong here
alras <- projectRaster(alras,crs=sr1) #orig latlong projection
xyz <- data.frame(xyFromCell(alras,1:ncell(alras)),v=getValues(alras))
xyz <- subset(xyz,complete.cases(xyz) & xyz$v>octhresh)
xyz$v<-NULL
#remove very close edge points
dthresh <- 30
dizs <- distm(eds1[,1:2],eds1[,1:2])/1000
rmv <- c()
for(i in 1:nrow(eds1)){
  if(i%in%rmv==F){
    rmv <- c(rmv,which(dizs[i,]>0 & dizs[i,]<dthresh))
  }
}
dizs <- dizs[-rmv,]
dizs <- dizs[,-rmv]
dim(dizs)
eds1 <- eds1[-rmv,]

#remove duplicates
eds1 <- eds1[!duplicated(eds1),]

#range orientation for starting cells
ras <- projectRaster(ras,crs=sr1) #return to latlongs for circular calcs
orient1 <- c() 
edsf <- data.frame(matrix(nrow=0,ncol=3))
names(edsf)<-c("x","y","clump")
for(x in 1:max(eds1$clump)){ #loop throigh each range patch
  edsc <- eds1[eds1$clump==x,1:2]
  if(nrow(edsc)>3){
    orient <- circular(rep(NA,nrow(edsc)),type="angle",units="degrees")
    #Now apply algorithm to get range edge orientation:
    brr <- circular(bearingRhumb(edsc[1,],edsc[2,])+runif(1,-0.01,0.01),units="degrees") 
    brr<-brr+90
    if(brr<0){brr <- brr+360}
    if(brr>360){brr <- brr-360}
    pts <- destPointRhumb(edsc[1,],brr,c(-25000,25000))
    exv <- extract(ras[[midlev]],pts)
    exv[is.na(exv)]<-0
    xxv <- 1:length(exv)
    if(brr<0){brr <- brr+360}
    orient[1]<-brr
    brr <- circular(bearingRhumb(edsc[(nrow(edsc)-1),],edsc[(nrow(edsc)),])+runif(1,-0.01,0.01),units="degrees")
    brr<-brr+90
    if(brr<0){brr <- brr+360}
    if(brr>360){brr <- brr-360}
    pts <- destPointRhumb(edsc[nrow(edsc),],brr,c(-25000,25000)) #search radius
    exv <- extract(ras[[midlev]],pts)
    exv[is.na(exv)]<-0
    xxv <- 1:length(exv)
    if(lm(exv~xxv)$coefficients[2]>0){
      if(brr<180){brr<-brr+180}else{brr<-brr-180}
    } 
    if(brr<0){brr <- brr+360}
    orient[nrow(edsc)]<-brr
    for(i in 2:(nrow(edsc)-1)){
      brr <- circular(bearingRhumb(edsc[i,],edsc[c(i-1,i+1),])+runif(2,-0.01,0.01),units="degrees")
      if(sum(is.na(brr))==1){
        brr <- brr[!is.na(brr)]+90
      }else if(sum(is.na(brr))==0){
        brr<-mean.circular(brr,units="degrees")
      }else{
        brr <- circular(0,units="degrees")
      }
      if(brr<0){brr <- brr+360}
      pts <- destPointRhumb(edsc[i,],brr,seq(-30000,30000,10000))#search radius
      exv <- extract(ras[[midlev]],pts)
      exv[is.na(exv)]<-0
      xxv <- 1:length(exv)
      if(lm(exv~xxv)$coefficients[2]>0){
        if(brr<180){brr<-brr+180}else{brr<-brr-180}
      } 
      if(brr<0){brr <- brr+360}
      orient[i]<-brr
  
    }  
    orient1 <- c(orient1,orient)
    edsf <- rbind(edsf,edsc)
  }
}
#########
eds1 <- edsf[,1:2]
orient1 <- circular(orient1,units="degrees")
orient1[orient1 <0] <- orient1[orient1 <0] +360
orient1[orient1 >360] <- orient1[orient1 >360] -360

# Now check whether edge cells are internal to range:
intern <- rep(0,nrow(eds1))
dizs <- distm(eds1,xyz)/1000
for(i in 1:nrow(eds1)){
  negb <- xyz[dizs[i,]<500 & dizs[i,]>0,]
  if(nrow(negb)>0){
    brs <- bearingRhumb(eds1[i,],negb) #bearings to all neighbours in window
    brs <- sort(brs)
    brs <- c(brs,brs[1]+360)
    print(max(diff(brs)))
    flush.console()
    intern[i] <-max(diff(brs))
  }
}

#= convert to albers
edsa <- eds1[,1:2]
coordinates(edsa)<- c("x","y")
proj4string(edsa)<-sr1
edsa <- spTransform(edsa,CRS(sr))

# Now create lines objects for range orientation projections:
brr <- orient1
lx <- lcells <- list()
edcel <- rep(NA,nrow(eds1))
#window1
for(i in 1:nrow(eds1)){
  edcel[i] <- as.numeric(extract(covzone,edsa[i],cellnumbers=T)[,"cells"])
  pls <- destPointRhumb(cbind(eds1$x[i],eds1$y[i]), brr[i], c(-d,0,d))
  ls2 <- Line(pls)
  ls2 <- Lines(list(ls2),ID=i)
  ls1 <- SpatialLines(list(ls2))
  proj4string(ls1)<-sr1
  ls1 <- spTransform(ls1,crs(sr))
  cc<-as.numeric(extract(covzone,ls1,cellnumbers=T,along=T)[[1]][,"cell"])
  cels <- unique(cc) #c(cc,adj))
  lcells[[i]] <-cels
  print(i)
  flush.console()
}

####################
nreps <- 1000# number of samples from posterior distribution to map
rcells <- cellFromXY(covzone2,cbind(datx$lon[datx$window==1],datx$lat[datx$window==1]))
bs<-covzone2

dthresh <- 250*1000 #distance to boundarty threshold

t1 <- proc.time()[3]
covzoneX<-covzone
# add exterior cells
for(i in 1:length(lcells)){
  covzone[lcells[[i]]]<-1
}
#get latlongs for cell centroids albers
xyz <- xyFromCell(covzone,1:ncell(covzone),spatial=T)
xyzdat <- data.frame(vs=extract(covzone,xyz),cell=1:ncell(covzone))
xyz <- spTransform(xyz,CRS(sr1))
xyz <- as(xyz,"data.frame")
xyz <- cbind(xyz,xyzdat)
xyz <- subset(xyz,!is.na(xyz$vs))

covzoneX1 <- projectRaster(covzoneX,crs=sr1)

#Object to store range change estimate data:
output <- data.frame(x=NA,y=NA,prob.index=NA,dirs=NA,changedist=NA,species=NA,survey=NA,window=NA,bprobchange=NA,intern=NA,iter=NA,point=NA)
output <- output[-1,]

for(t in 2:nwindows){
  ord1 <- sample(1:nrow(b.samps),nreps,replace=F)#randomly pick posterior samples
  ord2 <- sample(1:nrow(b.samps),nreps,replace=F)
  for(ri in 1:nreps){ #posterior samples loop
    print(ri)
    flush.console()
    outd <- data.frame(x=rep(NA,length(lcells)),y=NA,prob.index=NA,dirs=NA,changedist=NA,species=NA,survey=NA,window=NA,bprobchange=NA,intern=NA,iter=NA,point=NA)
    # Baseline first
    occp2 <- 1/(1+exp(-(X%*%as.numeric(b.samps[ord1[ri],]))))
    bs<-covzone2
    bs[rcells] <- occp2[datx$window==1]
    bs <- projectRaster(bs,crs=sr)
    bs <- resample(bs,covzone)*covzoneX
    bs[is.na(bs)]<-0
    bs[bs>=octhresh]<-1 #apply threshold
    bs[bs<octhresh]<-0
    rs<-covzone2
    occp2 <- 1/(1+exp(-(X%*%as.numeric(b.samps[ord2[ri],])))) #next window random sample
    rs[rcells] <- occp2[datx$window==t]
    rs <- projectRaster(rs,crs=sr)
    rs <- resample(rs,covzone)*covzoneX
    rs[is.na(rs)]<-0
    rs[rs>=octhresh]<-1
    rs[rs<octhresh]<-0
    
    #now extract presence absence along each range limit vector:
    for(w in 1:length(lcells)){
      ex1 <-bs[lcells[[w]]] 
      ex2 <-rs[lcells[[w]]] 
      for(e in 3:(length(ex1)-2)){
        if(ex1[e]==0 & sum(ex1[e-c(1,2)])>0 & sum(ex1[e+c(1,2)])>0){ex1[e]<-1}
        if(ex2[e]==0 & sum(ex2[e-c(1,2)])>0 & sum(ex2[e+c(1,2)])>0){ex2[e]<-1}
      }
      mdp<-which(lcells[[w]]==edcel[w])
      #find edge cell in new window:
      if(ex1[mdp]==1){#if inside range search out (expan nsion)
        opt<- which(ex1==0)
        if(sum(opt>mdp)>0){
          clp <- min(opt[opt>mdp])[1]-1
          pick <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }else{
          clp <- max(which(ex1==1))[1]
          pick <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }
      }else{# if outside range search back (extranction)
        opt<- which(ex1==1)
        if(sum(opt<mdp)>0){
          clp <- max(opt[opt<mdp])
          pick <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }else if(sum(opt>mdp)>0){
          clp <- min(opt)
          pick <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }else{
          pick<-eds1[w,1:2]
        }
      }  
      if(ex2[mdp]==1){#if inside range
        opt<- which(ex2==0)
        if(sum(opt>mdp)>0){
          clp <- min(opt[opt>mdp])-1
          pick2 <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }else{
          clp <- max(which(ex2==1))
          pick2 <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }
      }else{# if outside range
        opt<- which(ex2==1)
        if(sum(opt<mdp)>0){
          clp <- max(opt[opt<mdp])
          pick2 <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }else if(sum(opt>mdp)>0){
          clp <- min(opt)
          pick2 <- xyz[xyz$cell==lcells[[w]][clp],1:2]
        }else{
          pick2 <- eds1[w,1:2]
        }
      }  
      if(identical(as.numeric(round(pick,5)),as.numeric(round(pick2,5)))==F){
        nbr <- circular(bearingRhumb(pick,pick2),type="angle",units="degrees")
        if(range.circular(c(orient1[w],nbr))>140){ #determines if rrange expansion or contraction/extirpation
          extr<--1
         else{
          extr <- 1
        }	
      }else{
        extr<-1
      }
      if((orient1[w]>330 | orient1[w]<30)|(orient1[w]>60 & orient1[w]<120)|(orient1[w]>240 & orient1[w]<300)){card<-1.266}else{card<-1}#corrects for non-square cell dimensions
      outd$intern[w] <-intern[w]
      outd$dirs[w] <- orient1[w]
      outd$changedist[w] <- distRhumb(pick,pick2)*extr*card
      outd$x[w] <- pick[1,1]
      outd$y[w] <- pick[1,2]
      outd$window[w] <- t
      outd$survey[w]<-survey
      outd$species[w]<-p
      outd$iter[w]<-ri
      outd$point[w]<-w
    }
    # Now perform checks & cleanup errors:
    outs1 <- outd
    ################################
    #Find implausible shifts
    for(i in 1:nrow(outs1)){
      if(!is.na(outs1$changedist[i])){
        if(outs1$changedist[i]<(-50000)){ #long contraction
          pts <-destPointRhumb(cbind(outs1$x[i],outs1$y[i]),outs1$dirs[i],c(0,outs1$changedist[i]))
          ls2 <- Line(pts)
          ls2 <- Lines(list(ls2),ID=i)
          ls2 <- SpatialLines(list(ls2))
          proj4string(ls2)<-sr1
          ls2 <- spTransform(ls2,crs(sr))
          bvals<-extract(rs,ls2)[[1]]
          bvals2<-extract(bs,ls2)[[1]]
          if((mean(bvals>octhresh)>0.5)|(mean(bvals2<octhresh)>0.5)){ #implausible acontraction
            outs1$changedist[i] <- NA
          }
        }else if(outs1$changedist[i]>(50000)){ #long expansion 
          pts <-destPointRhumb(cbind(outs1$x[i],outs1$y[i]),outs1$dirs[i],c(0,outs1$changedist[i]))
          ls2 <- Line(pts)
          ls2 <- Lines(list(ls2),ID=i)
          ls2 <- SpatialLines(list(ls2))
          proj4string(ls2)<-sr1
          ls2 <- spTransform(ls2,crs(sr))
          bvals<-extract(bs,ls2)[[1]]
          bvals2<-extract(rs,ls2)[[1]]
          if((mean(bvals>octhresh)>0.5)|(mean(bvals2<octhresh)>0.5)){ #implausible for expansion
            outs1$changedist[i] <- NA
          }
        }
      }
    }
    #check longer changes for more plausible origins/destinations
    if(sum(!is.na(outs1$changedist))>0){
      xys1 <- destPointRhumb(cbind(outs1$x[!is.na(outs1$changedist)],outs1$y[!is.na(outs1$changedist)]),outs1$dirs[!is.na(outs1$changedist)],outs1$changedist[!is.na(outs1$changedist)])
      xys2 <- destPointRhumb(cbind(outs1$x[!is.na(outs1$changedist)],outs1$y[!is.na(outs1$changedist)]),outs1$dirs[!is.na(outs1$changedist)],0)
      
      for(i in 1:nrow(outs1)){
        if(!is.na(outs1$changedist[i])){
          if(abs(outs1$changedist[i])>(250000)){ #long change
            ds <- xys1[distRhumb(cbind(outs1$x[i],outs1$y[i]),xys1)<abs(outs1$changedist[i]),]
            if(length(ds)>4){
              mds <- distRhumb(cbind(outs1$x[i],outs1$y[i]),ds)
              outs1$changedist[i]<- ifelse(outs1$changedist[i]<0,-1,1)*mds[which.min(mds)[1]]
            }
          }
        }
      }
    }
    
   output <- rbind(output,outs1)
    
  }
  }
  #Save final range shift estimate posterior samples
  setwd(data)
  write.csv(output,paste(survey,"_win",t,"_",namlist$name[p],"_rangeshift_posterior_samples.csv",sep=""))
}

  proc.time()[3]-t1
	
