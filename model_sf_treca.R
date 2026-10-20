## Run analysis, write model results

## Before: data_sf.RData and data_mf.RData (data)
## After:  config_sf.RData, results_sf.RData, 
##         config_mf.RData and results_mf.RData (model)

rm(list=ls())

library(icesTAF)

library(FLCore)
library(stockassessment)
library(FLSAM)

# taf.session()

library(methods)


mkdir("model")

dataPath               <- file.path(".", "data")
modelPath              <- file.path(".",'model')

runName <- 'hom.treca_pelnord2025_sf'

source('utilities_model_config.R')

retroFlag <- T

### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_hom.pelnord.treca.Rdata'))

pg <- range(hom.treca)[2]

# stk <- hom.treca
# tun <- hom.treca

hom.treca.tun <- hom.treca.tun[c("idx.ssb")]#,"rec_russ",'idx.cpue'

hom.treca.ctrl <- config_sf_hom.treca_cecaf2025(hom.treca,
                                                hom.treca.tun,
                                                pg)
hom.treca.ctrl@residuals <- F

### ============================================================================
### Run single fleet model
### ============================================================================

# Run model
hom.treca.sam     <- FLSAM(hom.treca,
                           hom.treca.tun,
                           hom.treca.ctrl)

df.plot <- rec(hom.treca.sam)

ggplot(df.plot,aes(x=year,y=value))+
  geom_line()

# AIC(hom.treca.sam)

hom.treca@stock.n[] <- hom.treca.sam@stock.n[,ac(range(hom.treca)["minyear"]:range(hom.treca)["maxyear"])]
hom.treca@harvest[] <- hom.treca.sam@harvest[,ac(range(hom.treca)["minyear"]:range(hom.treca)["maxyear"])]

save(hom.treca,
     hom.treca.tun, 
     hom.treca.ctrl,
     hom.treca.sam, 
     file=file.path(modelPath,paste0(runName,'.RData')))

# dat <- subset(residuals(pil.sam),fleet=="catch unique")
# print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
#         geom_point(shape=21)+
#         theme_bw()+
#         theme(legend.position = "none")+
#         scale_color_manual(values = c("red", "blue")))

pdf(file.path(paste0(runName,"_fit.pdf")))
residual.diagnostics(hom.treca.sam)
dev.off()


### ============================================================================
### run single fleet retro
### ============================================================================

if(retroFlag){
  # retro
  n.retro.years <- 8  # Number of years for which to run the retrospective
  pil.ctrl@residuals <- FALSE
  pil.retro <- retro(pil.all,
                     pil.all.tun,
                     pil.ctrl, retro=n.retro.years)
  
  save(pil.retro,
       file=file.path(modelPath,paste0(runName,'_retro.RData')))
}
