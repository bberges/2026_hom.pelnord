
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

runName <- 'hom.trach_pelnord2025_sf'

source('utilities_model_config.R')

retroFlag <- T

### ============================================================================
### Construct single fleet control object
### ============================================================================

load(file.path(dataPath,'data_hom.pelnord.trach.Rdata'))

pg <- range(hom.trach)[2]

# stk <- hom.trach.pelnord
# tun <- hom.trach.tun

hom.trach.tun <- hom.trach.tun[c("idx.ssb","rec_russ")]#,'idx.cpue'

hom.trach.ctrl <- config_sf_hom.trach_cecaf2025(hom.trach,
                                                hom.trach.tun,
                                                pg)
hom.trach.ctrl@residuals <- F

### ============================================================================
### Run single fleet model
### ============================================================================

# Run model
hom.trach.sam     <- FLSAM(hom.trach,
                           hom.trach.tun,
                           hom.trach.ctrl)

df.plot <- as.data.frame(hom.trach.sam@harvest)


ggplot(df.plot,aes(x=year,y=data))+
  geom_line()+
  facet_wrap(~age,scales='free')

df.plot <- as.data.frame(quantMeans(hom.trach.sam@harvest[ac(2:5)]))

ggplot(subset(df.plot,year >= 2004),aes(x=year,y=data))+
  geom_line()+
  ylab('F2-5')

df.plot <- ssb(hom.trach.sam)

ggplot(subset(df.plot,year >= 2004),aes(x=year,y=value))+
  geom_line()

# AIC(pil.sam)

hom.trach@stock.n[] <- hom.trach.sam@stock.n[,ac(range(hom.trach)["minyear"]:range(hom.trach)["maxyear"])]
hom.trach@harvest[] <- hom.trach.sam@harvest[,ac(range(hom.trach)["minyear"]:range(hom.trach)["maxyear"])]

save(hom.trach,
     hom.trach.tun, 
     hom.trach.ctrl,
     hom.trach.sam, 
     file=file.path(modelPath,paste0(runName,'.RData')))

# dat <- subset(residuals(pil.sam),fleet=="catch unique")
# print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
#         geom_point(shape=21)+
#         theme_bw()+
#         theme(legend.position = "none")+
#         scale_color_manual(values = c("red", "blue")))

pdf(file.path(paste0(runName,"_fit.pdf")))
residual.diagnostics(hom.trach.sam)
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
