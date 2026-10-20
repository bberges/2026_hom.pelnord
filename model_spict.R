## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

rm(list=ls())

library(icesTAF)

library(FLCore)
library(spict)

# taf.session()

mkdir("model")

dataPath               <- file.path(".", "data")
modelPath              <- file.path(".",'model')

runName <- 'hom.treca_pelnord2025_sf'

source('utilities_spict.R')

retroFlag <- T

### ============================================================================
### load data
### ============================================================================

load(file.path(dataPath,'data_hom.pelnord.treca.Rdata'))

# hom.treca.tun <- hom.treca.tun[c("idx.ssb")]#,"rec_russ",'idx.cpue'

cab <- catch(hom.treca) # + min(c(catch(stk))[c(catch(stk)) > 0]) * 1
ind <- lapply(hom.treca.tun[c("idx.ssb","idx.cpue")], index)

# BUILD input list
inp <- buildTMBinputSPICT(cab, ind)

# ADD settings to input
inp <- spict.settings(inp)

fit <- fit.spict(inp,verbose = F)

windows()
spict::plotspict.biomass(fit)