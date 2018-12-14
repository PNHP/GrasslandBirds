###################################################################################################
## CUSTOM CODE FOR ANALYZING PENNSYLVANIA WILDS POINT COUNT DATA
## AUTHOR: Nicole Michel, National Audubon Society
## DATE: February 2016
## Updated December 2018 by Christopher Tracey, PNHP
## 
## SOME CODE ORIGINALLY MODIFIED FROM SOLYMOS ET AL. 2013 MEE
###################################################################################################

if (!requireNamespace("plyr", quietly=TRUE)) install.packages("plyr")
require(plyr)
if (!requireNamespace("here", quietly=TRUE)) install.packages("here")
require(here)
if (!requireNamespace("detect", quietly=TRUE)) install.packages("detect")
require(detect)
if (!requireNamespace("StreamMetabolism", quietly=TRUE)) install.packages("StreamMetabolism")
require(StreamMetabolism)
if (!requireNamespace("lme4", quietly=TRUE)) install.packages("lme4")
require(lme4)
if (!requireNamespace("AICcmodavg", quietly=TRUE)) install.packages("AICcmodavg")
require(AICcmodavg)
if (!requireNamespace("multcomp", quietly=TRUE)) install.packages("multcomp")
require(multcomp)
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
require(ggplot2)
if (!requireNamespace("glmmADMB", quietly=TRUE)) install.packages("glmmADMB")
require(glmmADMB)

# read in source code with functions
source(here("WesternPA_Density_Functions.R"))

## @@@@ USER INPUT REQUIRED HERE @@@@
# define the path to where your input and output data files are stored
filedirectory <- "grassland_data2018"

##############################################################################################
## read in data and calculate covariates (Jdate, time since local sunrise) 
## @@@@ ONLY RUN THIS SECTION ON NEW DATASETS. COMPLETE FOR 2015!! @@@@

rawdat <- read.csv(here(filedirectory,"bird_veg_1409_Grasslands_final.csv"), stringsAsFactors=FALSE, colClasses=c("date_start"="Date"))

###CHECK DATE FORMAT, change to M-D-Y###

# convert date to jdate
rawdat$Jdate <- as.POSIXlt(rawdat$date_start, format="%m-%d-%Y")$yday+1

# calculate local sunrise time -
rawdat$Sunrise <- NA
for (i in 1:nrow(rawdat)){
  if (!is.na(rawdat$Y[i])){
    temp <- as.POSIXlt(sunrise.set(lat=rawdat$Y[i], long=rawdat$X[i], date=format(as.Date(rawdat$date_start[i], format="%m-%d-%Y"), "%Y/%m/%d"), timezone = "UTC+4")$sunrise, format="%Y-%m-%d %h:%M:%s")
    rawdat$Sunrise[i] <- as.numeric(format(temp, "%H")) + as.numeric(format(temp, "%M"))/60
  }
}
# ignore the warnings

# convert start time to decimal time
rawdat$starttime <- as.numeric(substr(rawdat$time_start, 0, nchar(rawdat$time_start)-2)) + as.numeric(substr(rawdat$time_start, nchar(rawdat$time_start)-1, nchar(rawdat$time_start)))/60

# calculate time since local sunrise
rawdat$tslr <- rawdat$starttime - rawdat$Sunrise


##############################################################################################
## clean the data, remove unidentified species, keep only the survey visit with the highest count,
##   add zeroes for all species-point combinations where not observed
## @@@@ ONLY RUN THIS SECTION ON NEW DATASETS. @@@@

# get unique list of points and survey visits (unique_id)
UniqPtsVisits <- unique(rawdat$unique_id)

# remove records that have only outside or nosing records. Except keep points where no birds recorded)
rawdat <- subset(rawdat, tot_sing>0 )

# check to see if any records still have missing distances, remove if they still exist
rawdat[which(is.na(rawdat$distance)),]
rawdat <- rawdat[which(!is.na(rawdat$distance)),]

# get species list. check for entry errors (including use of lowercase letters), unidentified species
speclist <- unique(rawdat$elem_name)

# replace any lowercase letters in elem_name field with uppercase
rawdat$elem_name <- toupper(rawdat$elem_name)

# remove records of unidentified bird (UNBI), unidentified woodpecker (UNWO)
rawdat <- rawdat[which(!(rawdat$elem_name=="UNBI")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNWO")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNFL")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNWA")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNCH")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNDU")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNRA")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNHU")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNSP")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNSW")),]
rawdat <- rawdat[which(!(rawdat$elem_name=="UNTH")),]
rawdat$elem_name <- factor(rawdat$elem_name) # reset the factor list
speclist <- unique(rawdat$elem_name) # get HESPsed speclist


# add uniqueID field (unique by pt, spp, survey visit)
rawdat$uniqueID <- paste(rawdat$pt_id, rawdat$survey_vis, rawdat$elem_name, sep=".")

# get total count for each species * point * survey visit combination, summarized across distance bands
ptctTot <- ddply(rawdat, "uniqueID", summarise, sum(sing1 + sing2 + sing3 + sing4 + sing5))
colnames(ptctTot) <- c("uniqueID", "TotSing")
rawdat2 <- merge(ptctTot, rawdat, by.x=c("uniqueID"), by.y=c("uniqueID"))

# write cleaned file to csv
write.csv(rawdat2, file=paste(pathtofiles,"grasslands2018_JDate_TSLR_Cleaned_BothVisits.csv", sep=""))

# scroll through by point, keep only the data from the survey visit with the highest counts
rawdat2$specpt <- paste(rawdat2$elem_name, rawdat2$pt_id, sep=".") # add field with unique combo of species and point

# copy to temp file for identifying survey visit with highest tot_sing. Need to remove repeat lines for distance bands
temprawdat <- rawdat2[!(duplicated(rawdat2$uniqueID)),]

specptlist <- unique(temprawdat$specpt)

uniqueIDlist <- vector()
for (i in 1:length(specptlist)){
  temp <- subset(temprawdat, temprawdat$specpt==specptlist[i])
  maxct <- max(temp$tot_sing)
  if (nrow(temp[which(temp$TotSing==max(temp$TotSing))])==1){
    uniqueIDlist[i] <- temp$uniqueID[which(temp$TotSing==max(temp$TotSing))]
  } else if (nrow(temp[which(temp$TotSing==max(temp$TotSing))])==2){ # if TotSing is the same at both periods, use the first survey visit
    uniqueIDlist[i] <- temp$uniqueID[which(temp$survey_vis==1)]
  }
}

rawdatsv <- rawdat2[which(rawdat2$uniqueID %in% uniqueIDlist),]
rawdatsv$uniqueID <- factor(rawdatsv$uniqueID)


# scroll through by species, add in 0 records for points/visits where they were not detected
zeropts <- data.frame()
speclist <- unique(rawdatsv$elem_name)
for (s in speclist){
  temp <- subset(rawdatsv, rawdatsv$elem_name==s)
  misspt <- unique(rawdatsv$pt_id[which(!(rawdatsv$pt_id %in% temp$pt_id))])
  zeropts <- rbind(zeropts, data.frame(pt_id=misspt, elem_name=rep(s, length(misspt)), tot_sing=rep(0, length(misspt)), distance=rep(0, length(misspt)), sing1=rep(0, length(misspt)), sing2=rep(0, length(misspt)), sing3=rep(0, length(misspt)),sing4=rep(0, length(misspt)),sing5=rep(0, length(misspt)),specpt=paste(s,misspt,sep="."))) ## DY- added Sing1...    sing5=rep(0, length(misspt)), 
}


# merge with rawdat to fill in other columns
# Use c(1,3:15,19,21,23:125)] for 10min surveys, Use c(1,3:15,19,21:125)] for 5min surveys
# Use c(1:13,17,19,21:23,25:32) for peatlands
zerodat <- merge(zeropts, rawdatsv[!duplicated(rawdatsv$pt_id),c("pt_id","uniqueID","TotSing","OBJECTID","spec_pt","site_name","survey_yr","survey_vis","unique_id","date_start","time_start","sky","wind","temp","surveyors","nosing1","nosing2","nosing3","nosing4","nosing5","outside","tot_all","breed_code","comm_bird","All_detections","veg_surveyor","veg_date","grass_type","road","trail","invasive","soil_dist","equip_track","atv","mow","disc","fire","grass_cov","shrub_cov","bare_cov","grass_ht","shrub_ht","grass_pat","shrub_pat","bare_pat","comm_veg","cum_disturb","X","Y","Jdate","Sunrise","starttime","tslr")], by.x="pt_id", by.y="pt_id",all.x=T, all.y=F) #c(1,3:17,19,21,23,25,27:62)

write.csv(rawdatsv, file=paste(pathtofiles,"grasslands2018_rawdatsv.csv", sep=""))
write.csv(zeropts, file=paste(pathtofiles,"grasslands2018_zeropts.csv", sep=""))
write.csv(zerodat, file=paste(pathtofiles,"grasslands2018_zerodat.csv", sep=""))

zerodat <- zerodat[c("uniqueID","TotSing","OBJECTID","spec_pt","site_name","pt_id","survey_yr","survey_vis","unique_id","date_start","time_start","sky","wind","temp","surveyors","elem_name","distance","sing1","nosing1","sing2","nosing2","sing3","nosing3","sing4","nosing4","sing5","nosing5","outside","tot_all","tot_sing","breed_code","comm_bird","All_detections","veg_surveyor","veg_date","grass_type","road","trail","invasive","soil_dist","equip_track","atv","mow","disc","fire","grass_cov","shrub_cov","bare_cov","grass_ht","shrub_ht","grass_pat","shrub_pat","bare_pat","comm_veg","cum_disturb","X","Y","Jdate","Sunrise","starttime","tslr","specpt"  )]

# combine zero data with count data
cleandat <- rbind(rawdatsv, zerodat)
cleandat$uniqueID <- NULL

# write compiled, cleaned file to csv
write.csv(cleandat, file=paste(pathtofiles,"grasslands2018_JDate_TSLR_Cleaned_Zeros_10min.csv", sep=""))






# DENSITY ESTIMATES BEGIN HERE

##############################################################################################
## reshape the data into 3 datasets: a total count dataset for the GLMM, a removal dataset & design matrix, and a distance dataset & design matrix
## @@@@ ONLY RUN THIS SECTION ON NEW DATASETS. COMPLETE FOR 2015!! @@@@

# create total count dataset for GLMM with a single record for each point (remove extra distance bands)
ptctTot <- cleandat[!(duplicated(cleandat$specpt)),]
write.csv(ptctTot, file=paste(pathtofiles,"grassland2018_TotalBySpeciesPointSurvey.csv", sep=""))


# make removal data file (counts by time, collapsed across distance)
rem.sing1 <- ddply(cleandat, "specpt", summarise, sum(sing1))
rem.sing2 <- ddply(cleandat, "specpt", summarise, sum(sing2))
rem.sing3 <- ddply(cleandat, "specpt", summarise, sum(sing3))
rem.sing4 <- ddply(cleandat, "specpt", summarise, sum(sing4))
rem.sing5 <- ddply(cleandat, "specpt", summarise, sum(sing5))

colnames(rem.sing1) <- c("specpt", "TotSing1")
colnames(rem.sing2) <- c("specpt", "TotSing2")
colnames(rem.sing3) <- c("specpt", "TotSing3")
colnames(rem.sing4) <- c("specpt", "TotSing4")
colnames(rem.sing5) <- c("specpt", "TotSing5")

ptctRem <- merge(rem.sing1, cleandat, by.x="specpt", by.y="specpt")
ptctRem <- merge(rem.sing2, ptctRem, by.x="specpt", by.y="specpt")
ptctRem <- merge(rem.sing3, ptctRem, by.x="specpt", by.y="specpt")
ptctRem <- merge(rem.sing4, cleandat, by.x="specpt", by.y="specpt")
ptctRem <- merge(rem.sing5, ptctRem, by.x="specpt", by.y="specpt")
ptctRem <- ptctRem[!(duplicated(ptctRem$specpt)),]
write.csv(ptctRem, file=paste(pathtofiles,"grasslands2018_RemovalData.csv", sep=""))

# sum counts in duplicated distance records
cleandat$distance[which(cleandat$distance==100)] <- Inf # recode distances to detect format = max distance (50, 100, Inf)
cleandat$distance[which(cleandat$distance==50)] <- 1
cleandat$distance[which(cleandat$distance==0)] <- 0.5
cleandat$specptdist <- paste(cleandat$specpt, cleandat$distance, sep=".")
totsing <- ddply(cleandat, "specptdist", summarise, sum(sing1 + sing2 + sing3 + sing4 + sing5)) # use totsing to summarize across species whose data were split into multiple records
colnames(totsing) <- c("specptdist", "TotSingDist")
totsing <- merge(totsing, cleandat, by.x=c("specptdist"), by.y=c("specptdist"))
totsing <- totsing[!(duplicated(totsing$specptdist)),]

# make distance data file (counts summarized by distance, across time)
dist0 <- subset(totsing[which(totsing$distance==0.5),])
dist50 <- subset(totsing[which(totsing$distance==1),])
dist100 <- subset(totsing[which(totsing$distance==Inf),])

# add zero records for unique species records (by point and survey visit) with no observations at each distance category
dist0.add <- unique(totsing$specpt[which(!(totsing$specpt %in% dist0$specpt))])
dist0.add <- data.frame(specpt=dist0.add, TotSing=rep(0, length(dist0.add)),  TotSingDist=rep(0, length(dist0.add)), distance=rep(0.5, length(dist0.add)), sing1=rep(0, length(dist0.add)), sing2=rep(0, length(dist0.add)), sing3=rep(0, length(dist0.add)), sing4=rep(0, length(dist0.add)), sing5=rep(0, length(dist0.add)))
dist0.add$specptdist <- paste(dist0.add$specpt, dist0.add$distance, sep=".")

dist50.add <- unique(totsing$specpt[which(!(totsing$specpt %in% dist50$specpt))])
dist50.add <- data.frame(specpt=dist50.add, TotSing=rep(0, length(dist50.add)),  TotSingDist=rep(0, length(dist50.add)), distance=rep(1, length(dist50.add)), sing1=rep(0, length(dist50.add)), sing2=rep(0, length(dist50.add)), sing3=rep(0, length(dist50.add)), sing4=rep(0, length(dist50.add)), sing5=rep(0, length(dist50.add)))
dist50.add$specptdist <- paste(dist50.add$specpt, dist50.add$distance, sep=".")

dist100.add <- unique(totsing$specpt[which(!(totsing$specpt %in% dist100$specpt))])
dist100.add <- data.frame(specpt=dist100.add, TotSing=rep(0, length(dist100.add)),  TotSingDist=rep(0, length(dist100.add)), distance=rep(Inf, length(dist100.add)), sing1=rep(0, length(dist100.add)), sing2=rep(0, length(dist100.add)), sing3=rep(0, length(dist100.add)), sing4=rep(0, length(dist100.add)), sing5=rep(0, length(dist100.add)))
dist100.add$specptdist <- paste(dist100.add$specpt, dist100.add$distance, sep=".")





####### CODE ADAPTED TO HERE!!!! @@@@@@@




dist0.addm <- merge(dist0.add, totsing[,-c(1,2,3,17,18,20)], by.x=c("specpt"), by.y=c("specpt"), all.x=TRUE, all.y=FALSE)
dist0.addm <- dist0.addm[!duplicated(dist0.addm$specpt),]
dist50.addm <- merge(dist50.add, totsing[,-c(1,2,3,17,18,20)], by.x=c("specpt"), by.y=c("specpt"), all.x=TRUE, all.y=FALSE)
dist50.addm <- dist50.addm[!duplicated(dist50.addm$specpt),]
dist100.addm <- merge(dist100.add, totsing[,-c(1,2,3,17,18,20)], by.x=c("specpt"), by.y=c("specpt"), all.x=TRUE, all.y=FALSE)
dist100.addm <- dist100.addm[!duplicated(dist100.addm$specpt),]
ptctDist <- rbind(totsing, dist0.addm, dist50.addm, dist100.addm)
write.csv(ptctDist, file=paste(pathtofiles,"grasslands2018_DistanceData.csv", sep=""))


##############################################################################################
### FILL IN PA Comm Groups FROM DAVID'S CSV FILE
### @@@@ ONLY RUN THIS ONCE FOR DATA @@@@

ptctTot <- read.csv(file=paste(pathtofiles,"WPAData_TotalBySpeciesPointSurvey.csv", sep=""))
ptctRem <- read.csv(file=paste(pathtofiles,"WPAData_RemovalData.csv", sep=""))
ptctDist <- read.csv(file=paste(pathtofiles,"WPAData_DistanceData.csv", sep=""))

# fill in PA Comm Groups
PACommGroups <- read.csv(file=paste(pathtofiles,"PACommGroups2.csv", sep=""))
ptctTot <- merge(ptctTot, PACommGroups, by.x="pa_comm_code", by.y="pa_comm_code", all.x=T)
ptctRem <- merge(ptctRem, PACommGroups, by.x="pa_comm_code", by.y="pa_comm_code", all.x=T)
ptctDist <- merge(ptctDist, PACommGroups, by.x="pa_comm_code", by.y="pa_comm_code", all.x=T)

# export to CSV
write.csv(ptctTot, file=paste(pathtofiles,"WPAData_TotalBySpeciesPointSurvey.csv", sep=""))
write.csv(ptctRem, file=paste(pathtofiles,"WPAData_RemovalData.csv", sep=""))
write.csv(ptctDist, file=paste(pathtofiles,"WPAData_DistanceData.csv", sep=""))




##############################################################################################
### BUILD REMOVAL-ONLY MODELS TO IDENTIFY BEST COVARIATES FOR SINGING RATE (IF ANY)
### ASSUME HOMOGENEOUS SINGING RATES (INTERCEPT ONLY)

# HYPOTHESIZED COVARIATES AFFECTING SINGING RATE: JDAY, TSLR, JDAY^2, TSLR^2, SKY, WIND, TEMP

# BUILD DATA & DESIGN MATRICES: Y = COUNTS, D = DESIGN, X = COVARIATES
ptctRem <- read.csv(file=paste(pathtofiles,"grasslands2018_RemovalData.csv", sep=""))

### @@@@ USER INPUT REQUIRED HERE @@@@ 
# select study species and subset data to species. HERE: HESP
ptctRem.spp <- subset(ptctRem, ptctRem$elem_name=="HESP")

# Build Y matrix = counts (# rows = # records) by 2 columns (sing03, sing35)
Y.rem <- as.matrix(cbind(ptctRem.spp$TotSing1, ptctRem.spp$TotSing2, ptctRem.spp$TotSing3, ptctRem.spp$TotSing4, ptctRem.spp$TotSing5)) 

# Build D matrix = 2 columns each with the end time of the sampling period (i.e., 3, 5), # rows = # records
D.rem <- as.matrix(cbind(rep(1, nrow(Y.rem)), rep(2, nrow(Y.rem)), rep(3, nrow(Y.rem)), rep(4, nrow(Y.rem)), rep(5, nrow(Y.rem)))) 

# Build X matrix = one column for each covariate, # rows = # records
# Order: Jdate, Jdate^2, tslr, tslr^2, sky, wind, temp
X.remall <- as.matrix(cbind(ptctRem.spp$Jdate, (ptctRem.spp$Jdate^2), ptctRem.spp$tslr, (ptctRem.spp$tslr^2), ptctRem.spp$sky, ptctRem.spp$wind, ptctRem.spp$temp))

remdat <- list("Y"=Y.rem, "D"=D.rem, "X"=X.remall)

### @@@@ USER INPUT REQUIRED HERE @@@@ 
# save to rda file
save(remdat, file=paste(pathtofiles,"remdat_2018_HESP.Rdata", sep=""))

# this line calls up a function in the source code to run all possible combinations of removal models and return a sorted AIC table
# Note: just ignore the string of error messages. The usual tricks to silence the errors aren't working, and the errors are not fatal
# Note #2: all models are saved within the R environment and can be called up as long as your R session is open
# Note #3: please be patient, R can be slow. Time for a coffee break?
rem.aicmodtab <- WPARunRemovalModels(remdat=remdat)


### @@@@ USER INPUT REQUIRED HERE @@@@ 
# print model list, identify and select best model. HERE: SUNRISE + SUNRISE^2 + SKY
# jd = Jdate, ls = local sunrise (tslr), sky = sky, wd = wind, tm = temp
# jd2 = Jdate + Jdate^2; ls2 = tslr + tslr^2
head(rem.aicmodtab)

### @@@@ USER INPUT REQUIRED HERE @@@@ 
### look at best model - if SE are NA, work down the list of models until SE are produced. This should correct any errors in offsets later on.
summary(rem.jd)


##############################################################################################
### BUILD DISTANCEL-ONLY MODELS TO IDENTIFY BEST COVARIATES FOR EFFECTIVE DETECTION RADIUS (EDR; IF ANY)
### ASSUME HOMOGENEOUS EDR (INTERCEPT ONLY)

# HYPOTHESIZED COVARIATES AFFECTING EDR: PA Community Group, Forest Type, surveyors, wind 

# BUILD DATA & DESIGN MATRICES: Y = COUNTS, D = DESIGN, X = COVARIATES
ptctDist <- read.csv(file=paste(pathtofiles,"WPAData_DistanceData.csv", sep=""))

### @@@@ USER INPUT REQUIRED HERE @@@@ 
# select study species and subset data to species. HERE: HESP (forest spp, lots of data)
ptctDist.spp <- subset(ptctDist, ptctDist$elem_name=="HESP")

# Build Y matrix = counts (# rows = # records) by 3 columns (dist 0.5, dist 1, dist Inf)
Y.dist <- as.matrix(cbind(ptctDist.spp$TotSingDist[which(ptctDist.spp$distance==0.5)], ptctDist.spp$TotSingDist[which(ptctDist.spp$distance==1)], ptctDist.spp$TotSingDist[which(ptctDist.spp$distance==Inf)])) 

# Build D matrix = 3 columns each with the far edge of the distance band, measured in 100m increments (i.e., 0.5, 1, Inf), # rows = # records
D.dist <- as.matrix(cbind(rep(0.5, nrow(Y.dist)), rep(1, nrow(Y.dist)), rep(Inf, nrow(Y.dist)))) 

# Build X matrix = one column for each covariate, # rows = # records
# Order: for_type_code, pa_comm_group, surveyors, wind
# call up function from source code that converts categorical covariates to dummy-style coding
X.distall <- DummyCodeXDist(ptctDist.spp)

### @@@@ USER INPUT REQUIRED HERE @@@@ 
# save to rda file
distdat <- list("Y"=Y.dist, "D"=D.dist, "X"=X.distall)
save(distdat, file=paste(pathtofiles,"distdat_2015_HESP.Rdata", sep=""))


# this line calls up a function in the source code to run all possible combinations of distance models and return a sorted AIC table
# Note: just ignore the string of error messages. The usual tricks to silence the errors aren't working, and the errors are not fatal
# Note #2: all models are saved within the R environment and can be called up as long as your R session is open
# Note #3: please be patient, R can be slow. Time for another quick break?
dist.aicmodtab <- WPARunDistanceModels(distdat=distdat)

### @@@@ USER INPUT REQUIRED HERE @@@@ 
# print model list, identify and select best model. HERE: forest type
head(dist.aicmodtab)

### @@@@ USER INPUT REQUIRED HERE @@@@ 
### look at best model - if SE are NA, work down the list of models until SE are produced. This should correct any errors in offsets later on.
summary(dist.ft.sd)


rm(list=ls(pattern="GLMM"))

##############################################################################################
## NOW CALCULATE OFFSETS THEN ENTER PARAMETERS FROM REMOVAL AND DISTANCE MODELS INTO GLMMs TO ESTIMATE DENSITY


## IF YOU HAVE NOT READ IN THE DATA OR RUN THE BEST MODELS IN THE CURRENT R SESSION, DO SO NOW
ptctTot <- read.csv(file=paste(pathtofiles,"WPAData_TotalBySpeciesPointSurvey.csv", sep=""))
ptctRem <- read.csv(file=paste(pathtofiles,"WPAData_RemovalData.csv", sep=""))
ptctDist <- read.csv(file=paste(pathtofiles,"WPAData_DistanceData.csv", sep=""))


## NOTE: IF YOU CLOSED R SINCE YOU RAN THE REM AND DIST MODELS (ABOVE), YOU NEED TO RUN THEM AGAIN BEFORE PROCEEDING!

## @@@@ USER INPUT NEEDED HERE - enter names of best models here @@@@
bestmods <- list("remmod" = rem.sky.wd, "distmod"=dist.fg.sur)

## @@@@ USER INPUT NEEDED HERE - specify species and survey visit being modeled @@@@
ptctTot.spp <- subset(ptctTot, ptctTot$elem_name=="HESP")

# add columns with point count total duration (5 mins) and radius (Inf) to ptctTot
ptctTot.spp$dur <- rep(5, nrow(ptctTot.spp))
ptctTot.spp$dist <- rep(Inf, nrow(ptctTot.spp))

## calculate offsets to enter into GLMMs. 

## get offsets
det.offset <- with(ptctTot.spp, WPAGetOffsets(t=dur, r=dist, mods=bestmods))
summary(det.offset)
# check to see if the numbers look reasonable:
#  A = effective area sampled, in ha. Should be fairly small, depending on how loud the species is: <1 - ~3 ha, 
#  p = probability of a bird singing at least once during the point count. Must be <1, should usually be >0.5 (probably in the .8-.9 range)
#  q will always be 1 for unlimited radius point counts

# calculate offsets, associate with points and write to a csv file for use in habitat modelling
ptctTot.spp$DetOffset <- corrections2offset(det.offset)

# look at the calculated offsets
summary(ptctTot.spp$DetOffset)
# offsets should range from ~ -1 to maybe 1.5-2 but usually <1

write.csv(ptctTot.spp[,c("TotSing", "pt_id","DetOffset")], file=paste(pathtofiles,"DetOffset_HESP_2015.csv", sep=""))



#############################################################################################
### ENTER OFFSET INTO POISSON GLMMs

# If you've closed your R session since you ran the models and calculated offsets, run this code to 
#  read the data and offsets back in and compile the data to run the GLMMs
rm(list=ls(pattern="GLMM"))

### @@@@ USER INPUT REQUIRED HERE @@@@ 
ptctTot <- read.csv(file=paste(pathtofiles,"WPAData_TotalBySpeciesPointSurvey.csv", sep=""))
ptctTot.spp <- subset(ptctTot, ptctTot$elem_name=="HESP")
offsets <- read.csv(file=paste(pathtofiles,"DetOffset_HESP_2015.csv", sep=""))
ptctTot.spp <- merge(ptctTot.spp, offsets[,c("pt_id","DetOffset")], by.x="pt_id", by.y="pt_id")

### Create new fields with short codes for forest groups, forest communities, and site names (better for models)
ptctTot.spp$SiteCode <- factor(as.numeric(ptctTot.spp$site_name))

####### BUILD GLMM/GLM MODEL SET - FOREST GROUP, FOREST COMMUNITY, SITE AND NULL MIXED AND FIXED EFFECT MODELS
GLMM.ForestGroup <- glmer(TotSing ~ for_group_code  + (1|SiteCode), data=ptctTot.spp, family=poisson, offset=DetOffset)
GLMM.ForestType <- glmer(TotSing ~ for_typ_code  + (1|SiteCode), data=ptctTot.spp, family=poisson, offset=DetOffset)
GLMM.PACommType <- glmer(TotSing ~ pa_comm_code  + (1|SiteCode), data=ptctTot.spp, family=poisson, offset=DetOffset)
GLMM.PACommGroup <- glmer(TotSing ~ pa_comm_group  + (1|SiteCode), data=ptctTot.spp, family=poisson, offset=DetOffset)
GLMM.Site <- glmer(TotSing ~ SiteCode + (1|SiteCode), data=ptctTot.spp, family=poisson, offset=DetOffset)
GLMM.Null <- glmer(TotSing ~ 1 + (1|SiteCode), data=ptctTot.spp, family=poisson, offset=DetOffset)

# produce model table
candset.glmm <- list("GLMM.ForestGroup"=GLMM.ForestGroup, "GLMM.ForestType"=GLMM.ForestType, "GLMM.PACommType"=GLMM.PACommType, 
                     "GLMM.PACommGroup"=GLMM.PACommGroup, "GLMM.Site"=GLMM.Site, "GLMM.Null"=GLMM.Null)
GLMM.aictab <- aictab(cand.set=candset.glmm, second.ord=T)
GLMM.aictab

# first get estimated mean density across the entire study area using GLMM.Null
MeanDensity <- mean(predict(GLMM.Null, type="response"))
MeanDensity

# compare this to a very rough estimate based on the total number of birds encountered / area sampled
# don't report these numbers. This is just a check to see if modeled estimates are reasonable
RawDensity <- sum(ptctTot.spp$TotSing) / (length(unique(ptctTot.spp$pt_id)) * mean(det.offset$A))
RawDensity

# MeanDensity should be somewhat higher than RawDensity because it has been corrected for incomplete detection
#   e.g., if mean detection probability is 0.5, MeanDensity should be 2x higher than RawDensity

# test for overdispersion. Use one of the models for a forest class you want to estimate densities for.
#  If p < 0.05 data is overdispersed; consider using a negative binomial glmm
overdisp_fun(GLMM.ForestType)


###### FOREST GROUPS ##########
# do densities differ between forest groups? LRT: Yes if p < 0.05
anova(GLMM.ForestGroup, GLMM.Null)

# run post-hoc pairwise comparisons of forest groups
summary(glht(GLMM.ForestGroup, linfct=mcp(for_group_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest groups
newdat.fg<-data.frame(for_group_code=unique(ptctTot.spp$for_group_code))
mm<-model.matrix(~for_group_code,newdat.fg)
y<-mm%*%fixef(GLMM.ForestGroup)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.ForestGroup),mm))
tvar1 <- pvar1+VarCorr(GLMM.ForestGroup)$SiteCode[1]
newdat.fg <- data.frame(
  for_group_code=newdat.fg$for_group_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.fg$dum <- rep("x",nrow(newdat.fg))
newdat.fg<-newdat.fg[order(newdat.fg$for_group_code),]
newdat.fg$CIlo[which(is.na(newdat.fg$CIlo))] <- newdat.fg$y[which(is.na(newdat.fg$CIlo))]
newdat.fg$CIhi[which(is.na(newdat.fg$CIhi))] <- newdat.fg$y[which(is.na(newdat.fg$CIhi))]
newdat.fg$CIhi[which(newdat.fg$CIhi==Inf)] <- newdat.fg$y[which(newdat.fg$CIhi==Inf)]

# plot density and CIs for forest groups
ymax <- (max(newdat.fg$CIhi)*1.1)
ymin <- (min(newdat.fg$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.fg$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.fg$CIlo[which(newdat.fg$CIlo<0)] <- 0
}
Plotnewdat.fg <- ggplot(newdat.fg, aes(x=for_group_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.fg$for_group_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.fg


###### FOREST TYPES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.ForestType, GLMM.Null)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.ForestType, linfct=mcp(for_typ_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST TYPES
# get prediction intervals for forest types
newdat.ft<-data.frame(for_typ_code=unique(ptctTot.spp$for_typ_code))
mm<-model.matrix(~for_typ_code,newdat.ft)
y<-mm%*%fixef(GLMM.ForestType)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.ForestType),mm))
tvar1 <- pvar1+VarCorr(GLMM.ForestType)$SiteCode[1]
newdat.ft <- data.frame(
  for_typ_code=newdat.ft$for_typ_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.ft$dum <- rep("x",nrow(newdat.ft))
newdat.ft<-newdat.ft[order(newdat.ft$for_typ_code),]
newdat.ft$CIlo[which(is.na(newdat.ft$CIlo))] <- newdat.ft$y[which(is.na(newdat.ft$CIlo))]
newdat.ft$CIhi[which(is.na(newdat.ft$CIhi))] <- newdat.ft$y[which(is.na(newdat.ft$CIhi))]
newdat.ft$CIhi[which(newdat.ft$CIhi==Inf)] <- newdat.ft$y[which(newdat.ft$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.ft$CIhi)*1.1)
ymin <- (min(newdat.ft$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.ft$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.ft$CIlo[which(newdat.ft$CIlo<0)] <- 0
}
Plotnewdat.ft <- ggplot(newdat.ft, aes(x=for_typ_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.ft$for_typ_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.ft


###### PA COMMUNITY TYPES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.PACommType, GLMM.Null)

# run post-hoc pairwise comparisons of forest types - this takes a LONG time to run
summary(glht(GLMM.PACommType, linfct=mcp(pa_comm_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR PA Community Types
# get prediction intervals for forest types
newdat.pa<-data.frame(pa_comm_code=unique(ptctTot.spp$pa_comm_code))
mm<-model.matrix(~pa_comm_code,newdat.pa)
y<-mm%*%fixef(GLMM.PACommType)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.PACommType),mm))
tvar1 <- pvar1+VarCorr(GLMM.PACommType)$SiteCode[1]
newdat.pa <- data.frame(
  pa_comm_code=newdat.pa$pa_comm_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.pa$dum <- rep("x",nrow(newdat.pa))
newdat.pa<-newdat.pa[order(newdat.pa$pa_comm_code),]
newdat.pa$CIlo[which(is.na(newdat.pa$CIlo))] <- newdat.pa$y[which(is.na(newdat.pa$CIlo))]
newdat.pa$CIhi[which(is.na(newdat.pa$CIhi))] <- newdat.pa$y[which(is.na(newdat.pa$CIhi))]
newdat.pa$CIhi[which(newdat.pa$CIhi==Inf)] <- newdat.pa$y[which(newdat.pa$CIhi==Inf)]

# plot density and CIs for forest groups
ymax <- (max(newdat.pa$CIhi)*1.1)
ymin <- (min(newdat.pa$CIlo)*-20)
if (ymax <=0){
  ymax <- signif((max(newdat.pa$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.pa$CIlo[which(newdat.pa$CIlo<0)] <- 0
}
Plotnewdat.pa <- ggplot(newdat.pa, aes(x=pa_comm_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.pa$pa_comm_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.pa


###### PA COMMUNITY GROUPS ##########
# do densities differ between pa community groups? LRT: Yes if p < 0.05
anova(GLMM.PACommGroup, GLMM.Null)

# run post-hoc pairwise comparisons of forest types - this takes a LONG time to run
summary(glht(GLMM.PACommGroup, linfct=mcp(pa_comm_group="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR PA COMM GROUPS
# get prediction intervals for forest types
newdat.pag<-data.frame(pa_comm_group=unique(ptctTot.spp$pa_comm_group))
mm<-model.matrix(~pa_comm_group,newdat.pag)
y<-mm%*%fixef(GLMM.PACommGroup)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.PACommGroup),mm))
tvar1 <- pvar1+VarCorr(GLMM.PACommGroup)$SiteCode[1]
newdat.pag <- data.frame(
  pa_comm_group=newdat.pag$pa_comm_group,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.pag$dum <- rep("x",nrow(newdat.pag))
newdat.pag<-newdat.pag[order(newdat.pag$pa_comm_group),]
newdat.pag$CIlo[which(is.na(newdat.pag$CIlo))] <- newdat.pag$y[which(is.na(newdat.pag$CIlo))]
newdat.pag$CIhi[which(is.na(newdat.pag$CIhi))] <- newdat.pag$y[which(is.na(newdat.pag$CIhi))]
newdat.pag$CIhi[which(newdat.pag$CIhi==Inf)] <- newdat.pag$y[which(newdat.pag$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.pag$CIhi)*1.1)
ymin <- (min(newdat.pag$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.pag$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.pag$CIlo[which(newdat.pag$CIlo<0)] <- 0
}
Plotnewdat.pag <- ggplot(newdat.pag, aes(x=pa_comm_group, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.pag$pa_comm_group)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.pag

### @@@@ USER INPUT REQUIRED HERE @@@@
write.csv(newdat.pag, file = "HESP_density_pag.csv")

###### SITES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.Site, GLMM.Null)

# run post-hoc pairwise comparisons of forest types - this takes a LONG, LONG time to run
summary(glht(GLMM.Site, linfct=mcp(SiteCode="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.site<-data.frame(SiteCode=unique(ptctTot.spp$SiteCode))
mm<-model.matrix(~SiteCode,newdat.site)
y<-mm%*%fixef(GLMM.Site)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.Site),mm))
tvar1 <- pvar1+VarCorr(GLMM.Site)$SiteCode[1]
newdat.site <- data.frame(
  SiteCode=newdat.site$SiteCode,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.site$dum <- rep("x",nrow(newdat.site))
newdat.site<-newdat.site[order(newdat.site$SiteCode),]
newdat.site$CIlo[which(is.na(newdat.site$CIlo))] <- newdat.site$y[which(is.na(newdat.site$CIlo))]
newdat.site$CIhi[which(is.na(newdat.site$CIhi))] <- newdat.site$y[which(is.na(newdat.site$CIhi))]
newdat.site$CIhi[which(newdat.site$CIhi==Inf)] <- newdat.site$y[which(newdat.site$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.site$CIhi)*1.1)
ymin <- (min(newdat.site$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.site$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.site$CIlo[which(newdat.site$CIlo<0)] <- 0
}
Plotnewdat.site <- ggplot(newdat.site, aes(x=SiteCode, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.site$SiteCode)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.site


##############################################################
## NEGATIVE BINOMIAL GLMMs - only use if overdispersion test for Poisson models is significant



####### BUILD GLMM/GLM MODEL SET - FOREST GROUP, FOREST COMMUNITY, SITE AND NULL MIXED AND FIXED EFFECT MODELS
# fit non-zero-inflated NB models
GLMM.NB.ForestGroup <- glmmadmb(TotSing ~ for_group_code + offset(DetOffset) + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=FALSE)
GLMM.NB.ForestType <- glmmadmb(TotSing ~ for_typ_code + offset(DetOffset)  + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=FALSE)
GLMM.NB.PACommType <- glmmadmb(TotSing ~ pa_comm_code + offset(DetOffset)  + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=FALSE)
GLMM.NB.PACommGroup <- glmmadmb(TotSing ~ pa_comm_group + offset(DetOffset)  + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=FALSE)
GLMM.NB.Site <- glmmadmb(TotSing ~ SiteCode + offset(DetOffset) + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=FALSE)
GLMM.NB.Null <- glmmadmb(TotSing ~ 1  + offset(DetOffset)+ (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=FALSE)

# fit zero-inflated NB models
GLMM.NB.ForestGroup.ZI <- glmmadmb(TotSing ~ for_group_code + offset(DetOffset) + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=TRUE)
GLMM.NB.ForestType.ZI <- glmmadmb(TotSing ~ for_typ_code + offset(DetOffset)  + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=TRUE)
GLMM.NB.PACommType.ZI <- glmmadmb(TotSing ~ pa_comm_code + offset(DetOffset)  + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=TRUE)
GLMM.NB.PACommGroup.ZI <- glmmadmb(TotSing ~ pa_comm_group + offset(DetOffset)  + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=TRUE)
GLMM.NB.Site.ZI <- glmmadmb(TotSing ~ SiteCode + offset(DetOffset) + (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=TRUE)
GLMM.NB.Null.ZI <- glmmadmb(TotSing ~ 1  + offset(DetOffset)+ (1|SiteCode), data=ptctTot.spp, family="nbinom", link="log", zeroInflation=TRUE)

# produce model table. aictab doesn't work on glmmadmb objects, so build table manually
candset.glmm.nb <- list("GLMM.NB.ForestGroup"=GLMM.NB.ForestGroup, "GLMM.NB.ForestType"=GLMM.NB.ForestType, "GLMM.NB.PACommType"=GLMM.NB.PACommType,
                        "GLMM.NB.PACommGroup"=GLMM.NB.PACommGroup, "GLMM.NB.Site"=GLMM.NB.Site, "GLMM.NB.Null"=GLMM.NB.Null,
                        "GLMM.NB.ForestGroup.ZI"=GLMM.NB.ForestGroup.ZI, "GLMM.NB.ForestType.ZI"=GLMM.NB.ForestType.ZI, "GLMM.NB.PACommType.ZI"=GLMM.NB.PACommType.ZI,
                        "GLMM.NB.PACommGroup.ZI"=GLMM.NB.PACommGroup.ZI, "GLMM.NB.Site.ZI"=GLMM.NB.Site.ZI, "GLMM.NB.Null.ZI"=GLMM.NB.Null)
modnames.glmm.nb <- list("GLMM.NB.ForestGroup", "GLMM.NB.ForestType", "GLMM.NB.PACommType", 
                         "GLMM.NB.PACommGroup", "GLMM.NB.Site", "GLMM.NB.Null",
                         "GLMM.NB.ForestGroup.ZI", "GLMM.NB.ForestType.ZI", "GLMM.NB.PACommType.ZI",
                         "GLMM.NB.PACommGroup.ZI", "GLMM.NB.Site.ZI", "GLMM.NB.Null.ZI")

glmm.nb.aic <- lapply(candset.glmm.nb, AIC)
GLMM.NB.aictab   <- data.frame(ModName=modnames.glmm.nb, ModAIC=unlist(glmm.nb.aic))
GLMM.NB.aictab  = GLMM.NB.aictab[order(GLMM.NB.aictab$ModAIC),]
GLMM.NB.aictab$delta.AIC = GLMM.NB.aictab$ModAIC - min(GLMM.NB.aictab$ModAIC)
wt = exp(-0.5*GLMM.NB.aictab$delta.AIC)
GLMM.NB.aictab$Ak.wt = wt/sum(wt)
GLMM.NB.aictab$Ak.wt.cum = cumsum(GLMM.NB.aictab$Ak.wt)
GLMM.NB.aictab

# look at summary of best model
summary(GLMM.NB.PACommGroup)

# test for overdispersion. If p < 0.05 data is still overdispersed; 
overdisp_fun(GLMM.NB.PACommGroup)
###### If data still overdispersed, use GLMM in this order: 1)NB.ZI, 2)NB, 3) regular GLMM, 4) Null/mean density


###### FOREST GROUPS ##########
# do densities differ between forest groups? LRT: Yes if p < 0.05
anova(GLMM.NB.ForestGroup, GLMM.NB.Null)

# run post-hoc pairwise comparisons of forest groups
summary(glht(GLMM.NB.ForestGroup, linfct=mcp(for_group_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest groups
newdat.fg<-data.frame(for_group_code=unique(ptctTot.spp$for_group_code))
mm<-model.matrix(~for_group_code,newdat.fg)
y<-mm%*%fixef(GLMM.NB.ForestGroup)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.ForestGroup),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.ForestGroup)$SiteCode[1]
newdat.fg <- data.frame(
  for_group_code=newdat.fg$for_group_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.fg$dum <- rep("x",nrow(newdat.fg))
newdat.fg<-newdat.fg[order(newdat.fg$for_group_code),]
newdat.fg$CIlo[which(is.na(newdat.fg$CIlo))] <- newdat.fg$y[which(is.na(newdat.fg$CIlo))]
newdat.fg$CIhi[which(is.na(newdat.fg$CIhi))] <- newdat.fg$y[which(is.na(newdat.fg$CIhi))]
newdat.fg$CIhi[which(newdat.fg$CIhi==Inf)] <- newdat.fg$y[which(newdat.fg$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.fg$CIhi)*1.1)
ymin <- (min(newdat.fg$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.fg$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.fg$CIlo[which(newdat.fg$CIlo<0)] <- 0
}
Plotnewdat.fg <- ggplot(newdat.fg, aes(x=for_group_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.fg$for_group_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.fg


###### FOREST TYPES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.NB.ForestType, GLMM.NB.Null)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.NB.ForestType, linfct=mcp(for_typ_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.ft<-data.frame(for_typ_code=unique(ptctTot.spp$for_typ_code))
mm<-model.matrix(~for_typ_code,newdat.ft)
y<-mm%*%fixef(GLMM.NB.ForestType)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.ForestType),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.ForestType)$SiteCode[1]
newdat.ft <- data.frame(
  for_typ_code=newdat.ft$for_typ_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.ft$dum <- rep("x",nrow(newdat.ft))
newdat.ft<-newdat.ft[order(newdat.ft$for_typ_code),]
newdat.ft$CIlo[which(is.na(newdat.ft$CIlo))] <- newdat.ft$y[which(is.na(newdat.ft$CIlo))]
newdat.ft$CIhi[which(is.na(newdat.ft$CIhi))] <- newdat.ft$y[which(is.na(newdat.ft$CIhi))]
newdat.ft$CIhi[which(newdat.ft$CIhi==Inf)] <- newdat.ft$y[which(newdat.ft$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.ft$CIhi)*1.1)
ymin <- (min(newdat.ft$CIlo)*-0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.ft$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.ft$CIlo[which(newdat.ft$CIlo<0)] <- 0
}
Plotnewdat.ft <- ggplot(newdat.ft, aes(x=for_typ_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.ft$for_typ_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.ft


###### PA COMMUNITY TYPES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.NB.PACommType, GLMM.NB.Null)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.NB.PACommType, linfct=mcp(pa_comm_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.pa<-data.frame(pa_comm_code=unique(ptctTot.spp$pa_comm_code))
mm<-model.matrix(~pa_comm_code,newdat.pa)
y<-mm%*%fixef(GLMM.NB.PACommType)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.PACommType),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.PACommType)$SiteCode[1]
newdat.pa <- data.frame(
  pa_comm_code=newdat.pa$pa_comm_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.pa$dum <- rep("x",nrow(newdat.pa))
newdat.pa<-newdat.pa[order(newdat.pa$pa_comm_code),]
newdat.pa$CIlo[which(is.na(newdat.pa$CIlo))] <- newdat.pa$y[which(is.na(newdat.pa$CIlo))]
newdat.pa$CIhi[which(is.na(newdat.pa$CIhi))] <- newdat.pa$y[which(is.na(newdat.pa$CIhi))]
newdat.pa$CIhi[which(newdat.pa$CIhi==Inf)] <- newdat.pa$y[which(newdat.pa$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.pa$CIhi)*1.1)
ymin <- (min(newdat.pa$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.pa$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.pa$CIlo[which(newdat.pa$CIlo<0)] <- 0
}
Plotnewdat.pa <- ggplot(newdat.pa, aes(x=pa_comm_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.pa$pa_comm_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.pa



###### PA COMMUNITY GROUPS ##########
# do densities differ between pa community groups? LRT: Yes if p < 0.05
anova(GLMM.NB.PACommGroup, GLMM.NB.Null)

# run post-hoc pairwise comparisons of forest types - this takes a LONG time to run
summary(glht(GLMM.NB.PACommGroup, linfct=mcp(pa_comm_group="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.pag<-data.frame(pa_comm_group=unique(ptctTot.spp$pa_comm_group))
mm<-model.matrix(~pa_comm_group,newdat.pag)
y<-mm%*%fixef(GLMM.NB.PACommGroup)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.PACommGroup),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.PACommGroup)$SiteCode[1]
newdat.pag <- data.frame(
  pa_comm_group=newdat.pag$pa_comm_group,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.pag$dum <- rep("x",nrow(newdat.pag))
newdat.pag<-newdat.pag[order(newdat.pag$pa_comm_group),]
newdat.pag$CIlo[which(is.na(newdat.pag$CIlo))] <- newdat.pag$y[which(is.na(newdat.pag$CIlo))]
newdat.pag$CIhi[which(is.na(newdat.pag$CIhi))] <- newdat.pag$y[which(is.na(newdat.pag$CIhi))]
newdat.pag$CIhi[which(newdat.pag$CIhi==Inf)] <- newdat.pag$y[which(newdat.pag$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.pag$CIhi)*1.1)
ymin <- (min(newdat.pag$CIlo)*-0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.pag$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.pag$CIlo[which(newdat.pag$CIlo<0)] <- 0
}
Plotnewdat.pag <- ggplot(newdat.pag, aes(x=pa_comm_group, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.pag$pa_comm_group)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.pag

### @@@@ USER INPUT REQUIRED HERE @@@@
write.csv(newdat.pag, file = "HESP_density_pag.csv")

###### Sites ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.NB.Site, GLMM.NB.Null)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.NB.Site, linfct=mcp(SiteCode="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.site<-data.frame(SiteCode=unique(ptctTot.spp$SiteCode))
mm<-model.matrix(~SiteCode,newdat.site)
y<-mm%*%fixef(GLMM.NB.Site)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.Site),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.Site)$SiteCode[1]
newdat.site <- data.frame(
  SiteCode=newdat.site$SiteCode,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.site$dum <- rep("x",nrow(newdat.site))
newdat.site<-newdat.site[order(newdat.site$SiteCode),]
newdat.site$CIlo[which(is.na(newdat.site$CIlo))] <- newdat.site$y[which(is.na(newdat.site$CIlo))]
newdat.site$CIhi[which(is.na(newdat.site$CIhi))] <- newdat.site$y[which(is.na(newdat.site$CIhi))]
newdat.site$CIhi[which(newdat.site$CIhi==Inf)] <- newdat.site$y[which(newdat.site$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.site$CIhi)*1.1)
ymin <- (min(newdat.site$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.site$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.site$CIlo[which(newdat.site$CIlo<0)] <- 0
}
Plotnewdat.site <- ggplot(newdat.site, aes(x=SiteCode, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.site$SiteCode)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.site

##########Zero-inflated Negative Binomial GLMM###############

###### FOREST GROUPS ##########
# do densities differ between forest groups? LRT: Yes if p < 0.05
anova(GLMM.NB.ForestGroup.ZI, GLMM.NB.Null.ZI)

# run post-hoc pairwise comparisons of forest groups
summary(glht(GLMM.NB.ForestGroup.ZI, linfct=mcp(for_group_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest groups
newdat.fg<-data.frame(for_group_code=unique(ptctTot.spp$for_group_code))
mm<-model.matrix(~for_group_code,newdat.fg)
y<-mm%*%fixef(GLMM.NB.ForestGroup.ZI)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.ForestGroup.ZI),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.ForestGroup)$SiteCode[1]
newdat.fg <- data.frame(
  for_group_code=newdat.fg$for_group_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.fg$dum <- rep("x",nrow(newdat.fg))
newdat.fg<-newdat.fg[order(newdat.fg$for_group_code),]
newdat.fg$CIlo[which(is.na(newdat.fg$CIlo))] <- newdat.fg$y[which(is.na(newdat.fg$CIlo))]
newdat.fg$CIhi[which(is.na(newdat.fg$CIhi))] <- newdat.fg$y[which(is.na(newdat.fg$CIhi))]
newdat.fg$CIhi[which(newdat.fg$CIhi==Inf)] <- newdat.fg$y[which(newdat.fg$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.fg$CIhi)*1.1)
ymin <- (min(newdat.fg$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.fg$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.fg$CIlo[which(newdat.fg$CIlo<0)] <- 0
}
Plotnewdat.fg <- ggplot(newdat.fg, aes(x=for_group_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.fg$for_group_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.fg


###### FOREST TYPES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.NB.ForestType.ZI, GLMM.NB.Null.ZI)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.NB.ForestType.ZI, linfct=mcp(for_typ_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.ft<-data.frame(for_typ_code=unique(ptctTot.spp$for_typ_code))
mm<-model.matrix(~for_typ_code,newdat.ft)
y<-mm%*%fixef(GLMM.NB.ForestType.ZI)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.ForestType.ZI),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.ForestType.ZI)$SiteCode[1]
newdat.ft <- data.frame(
  for_typ_code=newdat.ft$for_typ_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.ft$dum <- rep("x",nrow(newdat.ft))
newdat.ft<-newdat.ft[order(newdat.ft$for_typ_code),]
newdat.ft$CIlo[which(is.na(newdat.ft$CIlo))] <- newdat.ft$y[which(is.na(newdat.ft$CIlo))]
newdat.ft$CIhi[which(is.na(newdat.ft$CIhi))] <- newdat.ft$y[which(is.na(newdat.ft$CIhi))]
newdat.ft$CIhi[which(newdat.ft$CIhi==Inf)] <- newdat.ft$y[which(newdat.ft$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.ft$CIhi)*1.1)
ymin <- (min(newdat.ft$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.ft$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.ft$CIlo[which(newdat.ft$CIlo<0)] <- 0
}
Plotnewdat.ft <- ggplot(newdat.ft, aes(x=for_typ_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.ft$for_typ_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.ft


###### PA COMMUNITY TYPES ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.NB.PACommType.ZI, GLMM.NB.Null.ZI)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.NB.PACommType.ZI, linfct=mcp(pa_comm_code="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.pa<-data.frame(pa_comm_code=unique(ptctTot.spp$pa_comm_code))
mm<-model.matrix(~pa_comm_code,newdat.pa)
y<-mm%*%fixef(GLMM.NB.PACommType.ZI)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.PACommType.ZI),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.PACommType.ZI)$SiteCode[1]
newdat.pa <- data.frame(
  pa_comm_code=newdat.pa$pa_comm_code,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.pa$dum <- rep("x",nrow(newdat.pa))
newdat.pa<-newdat.pa[order(newdat.pa$pa_comm_code),]
newdat.pa$CIlo[which(is.na(newdat.pa$CIlo))] <- newdat.pa$y[which(is.na(newdat.pa$CIlo))]
newdat.pa$CIhi[which(is.na(newdat.pa$CIhi))] <- newdat.pa$y[which(is.na(newdat.pa$CIhi))]
newdat.pa$CIhi[which(newdat.pa$CIhi==Inf)] <- newdat.pa$y[which(newdat.pa$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.pa$CIhi)*1.1)
ymin <- (min(newdat.pa$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.pa$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.pa$CIlo[which(newdat.pa$CIlo<0)] <- 0
}
Plotnewdat.pa <- ggplot(newdat.pa, aes(x=pa_comm_code, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.pa$pa_comm_code)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.pa



###### PA COMMUNITY GROUPS ##########
# do densities differ between pa community groups? LRT: Yes if p < 0.05
anova(GLMM.NB.PACommGroup.ZI, GLMM.NB.Null.ZI)

# run post-hoc pairwise comparisons of forest types - this takes a LONG time to run
summary(glht(GLMM.NB.PACommGroup.ZI, linfct=mcp(pa_comm_group="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.pag<-data.frame(pa_comm_group=unique(ptctTot.spp$pa_comm_group))
mm<-model.matrix(~pa_comm_group,newdat.pag)
y<-mm%*%fixef(GLMM.NB.PACommGroup.ZI)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.PACommGroup.ZI),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.PACommGroup.ZI)$SiteCode[1]
newdat.pag <- data.frame(
  pa_comm_group=newdat.pag$pa_comm_group,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.pag$dum <- rep("x",nrow(newdat.pag))
newdat.pag<-newdat.pag[order(newdat.pag$pa_comm_group),]
newdat.pag$CIlo[which(is.na(newdat.pag$CIlo))] <- newdat.pag$y[which(is.na(newdat.pag$CIlo))]
newdat.pag$CIhi[which(is.na(newdat.pag$CIhi))] <- newdat.pag$y[which(is.na(newdat.pag$CIhi))]
newdat.pag$CIhi[which(newdat.pag$CIhi==Inf)] <- newdat.pag$y[which(newdat.pag$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.pag$CIhi)*1.1)
ymin <- (min(newdat.pag$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.pag$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.pag$CIlo[which(newdat.pag$CIlo<0)] <- 0
}
Plotnewdat.pag <- ggplot(newdat.pag, aes(x=pa_comm_group, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.pag$pa_comm_group)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.pag

### @@@@ USER INPUT REQUIRED HERE @@@@
write.csv(newdat.pag, file = "HESP_density_pag.csv")

###### Sites ##########
# do densities differ between forest types? LRT: Yes if p < 0.05
anova(GLMM.NB.Site, GLMM.NB.Null)

# run post-hoc pairwise comparisons of forest types
summary(glht(GLMM.NB.Site, linfct=mcp(SiteCode="Tukey")))

### PLOT PREDICTED DENSITIES AND CIs FOR FOREST GROUPS
# get prediction intervals for forest types
newdat.site<-data.frame(SiteCode=unique(ptctTot.spp$SiteCode))
mm<-model.matrix(~SiteCode,newdat.site)
y<-mm%*%fixef(GLMM.NB.Site)
pvar1 <- diag(mm %*% tcrossprod(vcov(GLMM.NB.Site),mm))
tvar1 <- pvar1+VarCorr(GLMM.NB.Site)$SiteCode[1]
newdat.site <- data.frame(
  SiteCode=newdat.site$SiteCode,
  y=exp(y),
  CIlo = exp(y-1.96*sqrt(pvar1))   # Lower 95%  CI without random effects 
  , CIhi = exp(y+1.96*sqrt(pvar1))  # Upper 95%  CI without random effects 
  , tlo = exp(y-1.96*sqrt(tvar1))  # Lower 95%  CI with random effects 
  , thi = exp(y+1.96*sqrt(tvar1))  # Upper 95%  CI with random effects 
)
newdat.site$dum <- rep("x",nrow(newdat.site))
newdat.site<-newdat.site[order(newdat.site$SiteCode),]
newdat.site$CIlo[which(is.na(newdat.site$CIlo))] <- newdat.site$y[which(is.na(newdat.site$CIlo))]
newdat.site$CIhi[which(is.na(newdat.site$CIhi))] <- newdat.site$y[which(is.na(newdat.site$CIhi))]
newdat.site$CIhi[which(newdat.site$CIhi==Inf)] <- newdat.site$y[which(newdat.site$CIhi==Inf)]

# plot density and CIs for forest types
ymax <- (max(newdat.site$CIhi)*1.1)
ymin <- (min(newdat.site$CIlo)*0.9)
if (ymax <=0){
  ymax <- signif((max(newdat.site$CIhi)*1.1), digits=1)
} 
if (ymin <0) {
  ymin <- -ymax
  breakvec <- c(signif(ymin, digits=1), signif((ymin/2), digits=1), 0, signif((ymax/2), digits=1), signif(ymax, digits=1))
} else {
  ymin <- 0
  breakvec <- c(0, signif((ymax/2), digits=1), signif(ymax, digits=1))
  newdat.site$CIlo[which(newdat.site$CIlo<0)] <- 0
}
Plotnewdat.site <- ggplot(newdat.site, aes(x=SiteCode, y=y)) + 
  geom_point(aes(shape=dum, size=10)) +
  scale_shape_manual(values=c(1,0,15)) +
  geom_errorbar(aes(ymin=(CIlo), ymax=(CIhi)), size=1, width=0) +
  scale_fill_identity()+
  scale_x_discrete(labels=newdat.site$SiteCode)+
  xlab("")+
  ylab("Singing male density per ha") +
  scale_y_continuous(expand=c(0,0), limits=c(ymin, ymax), breaks=breakvec) +   
  theme(axis.text.x=element_text(size=10, colour="#000000"), axis.text.y=element_text(size=10, colour="#000000"),
        axis.title.x=element_text(size=10, colour="#000000"), axis.title.y=element_text(size=10, vjust=1.0),
        panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), plot.title=element_text(face="bold", size=10, vjust=0.9),
        legend.position="none")
Plotnewdat.site