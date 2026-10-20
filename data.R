### ============================================================================
### library loading and setting up paths
### ============================================================================

rm(list=ls())

library(icesTAF)
library(FLCore)
library(tidyverse)
library(FSA)

mkdir("data")

dataPath     <- file.path(".", "data")
dataSource   <- file.path(".","bootstrap", "data")

### ============================================================================
### read data
### ============================================================================

caton           <- read.csv(file.path(dataSource,'caton.csv'),sep = ',',check.names = F)
canum.trach     <- read.csv(file.path(dataSource,'canum_trachurus.csv'),sep = ',',check.names = F)
canum.treca     <- read.csv(file.path(dataSource,'canum_trecae.csv'),sep = ',',check.names = F)
survey_rec      <- read.csv(file.path(dataSource,'survey_rec.csv'),sep = ',',check.names = F)
survey_cpue     <- read.csv(file.path(dataSource,'cpue_russe.csv'),sep = ',',check.names = F)
survey_ssb      <- read.csv(file.path(dataSource,'survey_nansen_alamir_biomass.csv'),sep = ',',check.names = F)
matprop         <- read.csv(file.path(dataSource,'matprop.csv'),sep = ';')
west            <- read.csv(file.path(dataSource,'west.csv'),sep = ';',check.names = F)
natmort         <- read.csv(file.path(dataSource,'natmort.csv'),sep = ';',check.names = F)

### ============================================================================
### prepare the FLStock object
### ============================================================================
minYear <- 1990
maxYear <- 2024
yearVec <- minYear:maxYear
nYears  <- length(yearVec)
minAge  <- 0
maxAge  <- 9
ageVec  <- minAge:maxAge
nAges   <- length(ageVec)
pg <- 8

#### stock object
FLR_temp <- FLQuant( array(0,dim = c(nAges,nYears,1,1,1,1)),
                     dimnames=list( age=as.character(ageVec),
                                    year=as.character(yearVec),
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'),
                     units='1e6')

hom.trach <- FLStock(name = 'COPACE_hom.trach', 
                       desc = 'COPACE Pel-Nord hom.trach',
                       FLR_temp)

hom.treca <- FLStock(name = 'COPACE_hom.treca', 
                             desc = 'COPACE Pel-Nord hom.treca',
                             FLR_temp)

# caton
hom.trach@catch[]  <- an(caton$trachurus)
hom.trach@landings <- hom.trach@catch
hom.trach@discards[] <- 0

hom.treca@catch[]  <- an(caton$trecae)
hom.treca@landings <- hom.treca@catch
hom.treca@discards[] <- 0

# weca
hom.trach@catch.wt[] <- t(as.matrix(west[,2:dim(west)[2]]))
hom.treca@catch.wt[] <- t(as.matrix(west[,2:dim(west)[2]]))

hom.trach@landings.wt <- hom.trach@catch.wt
hom.trach@discards.wt <- hom.trach@catch.wt

hom.treca@landings.wt <- hom.treca@catch.wt
hom.treca@discards.wt <- hom.treca@catch.wt

# west
hom.trach@stock.wt[] <- t(as.matrix(west[,2:dim(west)[2]]))
hom.treca@stock.wt[] <- t(as.matrix(west[,2:dim(west)[2]]))

# canum
hom.trach@catch.n[ac(0:8)] <- as.matrix(canum.trach[,2:dim(canum.trach)[2]])
hom.treca@catch.n[ac(0:8)] <- as.matrix(canum.treca[,2:dim(canum.treca)[2]])

hom.trach@landings.n <- hom.trach@catch.n
hom.trach@discards.n[] <- 0

hom.treca@landings.n <- hom.treca@catch.n
hom.treca@discards.n[] <- 0

# mat
hom.trach@mat[] <- as.matrix(matprop[,2])
hom.treca@mat[] <- as.matrix(matprop[,2])

# spawning time
hom.trach@m.spwn[]        <- 0.1
hom.trach@harvest.spwn[]  <- 0.1

hom.treca@m.spwn[]        <- 0.1
hom.treca@harvest.spwn[]  <- 0.1

hom.trach@m[] <- natmort$m
hom.treca@m[] <- natmort$m

range(hom.trach)[c('minfbar', 'maxfbar')] <- c(3,9)
range(hom.treca)[c('minfbar', 'maxfbar')] <- c(3,9)

# set plus group
hom.trach <- setPlusGroup(hom.trach,pg,na.rm=T)
hom.treca <- setPlusGroup(hom.treca,pg,na.rm=T)

# start ages to 1
# hom.trach <- hom.trach[ac(1:pg),]
# hom.treca     <- hom.treca[ac(1:pg),]

### ============================================================================
### prepare the FLIndices
### ============================================================================

### rec survey
FLR_temp <- FLQuant(  array(0,dim = c(2,length(unique(unique(survey_rec$Year))),1,1,1,1)),
                              dimnames=list(age=as.character(c(0,1)),
                                            year=unique(unique(survey_rec$Year)),
                                            unit=c("1e6"),
                                            season='all',
                                            area='all'))
rec_russ.trach <- FLIndex(name = 'russian_rec', 
                          desc = 'Russian recruitment survey trachurus',
                          index=FLR_temp)
survey_rec.trach <- t(subset(survey_rec,Species == 'trachurus') %>% select(-c(Species)))
rec_russ.trach@index[] <- as.matrix(survey_rec.trach[2:3,])
range(rec_russ.trach)[c('startf', 'endf')] <- c(0.8,0.9) # Oct-Nov
#rec_russ@index@.Data[which(is.na(rec_russ@index))] <- -1
rec_russ.trach <- rec_russ.trach[ac(0:1),] # reduce to age 1

rec_russ.treca <- FLIndex(name = 'russian_rec', 
                          desc = 'Russian recruitment survey trecae',
                          index=FLR_temp)
survey_rec.treca <- t(subset(survey_rec,Species == 'trecae') %>% select(-c(Species)))
rec_russ.treca@index[] <- as.matrix(survey_rec.treca[2:3,])
range(rec_russ.treca)[c('startf', 'endf')] <- c(0.8,0.9) # Oct-Nov
#rec_russ@index@.Data[which(is.na(rec_russ@index))] <- -1
rec_russ.treca <- rec_russ.treca[ac(0:1),] # reduce to age 1

### ssb survey combining nansen and marocan vessel
FLR_temp <- FLQuant(  array(0,dim = c(1,length(unique(survey_ssb$year)),1,1,1,1)),
                      dimnames=list(age='all',
                                    year=ac(unique(survey_ssb$year)),
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))
# trachurus
idx.ssb.trach <- FLIndex(name = 'nansen+MA_ssb_trachurus',
                         desc = 'nansen+MA SSB index trachurus',
                         index=FLR_temp)
idx.ssb.trach@index[] <- survey_ssb$idx_trachurus
range(idx.ssb.trach)[c('startf', 'endf')] <- c(0.75,0.9) # Sept-Nov

# trecae
idx.ssb.treca <- FLIndex(name = 'nansen+MA_ssb_trecae',
                  desc = 'nansen+MA SSB index trecae',
                  index=FLR_temp)
idx.ssb.treca@index[] <- survey_ssb$idx_trecae
range(idx.ssb.treca)[c('startf', 'endf')] <- c(0.75,0.9) # Sept-Nov

### CPUE indices
FLR_temp <- FLQuant(  array(0,dim = c(1,length(unique(survey_cpue$Year)),1,1,1,1)),
                      dimnames=list(age='all',
                                    year=ac(unique(survey_cpue$Year)),
                                    unit=c("1e6"),
                                    season='all',
                                    area='all'))
# trachurus
cpue.ssb.trach <- FLIndex(name = 'cpue_trachurus',
                         desc = 'cpue index trachurus',
                         index=FLR_temp)
cpue.ssb.trach@index[] <- survey_cpue$`Trachurus trachurus`
range(cpue.ssb.trach)[c('startf', 'endf')] <- c(0.75,0.9) # Sept-Nov

# trecae
cpue.ssb.treca <- FLIndex(name = 'cpue_trecaes',
                          desc = 'cpue index trecae',
                          index=FLR_temp)
cpue.ssb.treca@index[] <- survey_cpue$`Trachurus trecae`
range(cpue.ssb.treca)[c('startf', 'endf')] <- c(0.75,0.9) # Sept-Nov

### compile indices trachurus
hom.trach.tun <- FLIndices(idx.ssb=idx.ssb.trach,
                           idx.cpue=cpue.ssb.trach,
                           rec_russ=rec_russ.trach)
hom.trach.tun   <- lapply(hom.trach.tun,function(x) {x@type <- "number"; return(x)})

### compile indices trecae
hom.treca.tun <- FLIndices(idx.ssb=idx.ssb.treca,
                           idx.cpue=cpue.ssb.treca,
                           rec_russ=rec_russ.treca)
hom.treca.tun   <- lapply(hom.treca.tun,function(x) {x@type <- "number"; return(x)})


### ============================================================================
### save objects
### ============================================================================

save(hom.treca,
     hom.treca.tun, 
     file=file.path(dataPath,'data_hom.pelnord.treca.Rdata'))

save(hom.trach,
     hom.trach.tun, 
     file=file.path(dataPath,'data_hom.pelnord.trach.Rdata'))

