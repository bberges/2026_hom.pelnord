
rm(list=ls())

library(icesTAF)
library(FLCore)
library(stockassessment)
library(FLSAM)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(scales)

mkdir("report")

dataSource     <- file.path(".", "bootstrap",'data')
dataPath       <- file.path(".", "data")
modelPath      <- file.path(".", "model")
reportPath     <- file.path(".", "report")

source('utilities_plots.R')

runName <- 'pil_pelnord2025_sf'

load(file.path(modelPath,paste0(runName,'.RData')))
# load(file.path(modelPath,paste0(runName,'_sf_retro.RData')))

#------------------------------------------
### input data

prefix <- 'input'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))


## prop catches at age since 2000
taf.png(paste0(prefix,"_prop_catches"))
prop_catches <- as.data.frame(pil.all@catch.n)
prop_catches$data[prop_catches$data == -1] <- NA

ggplot(subset(prop_catches,year>=2000),aes(x=year,y=data))+
  geom_bar(aes(fill = as.factor(age)),stat="identity",position = "fill")+
  theme_bw()+
  labs(x = 'Year',
       y = 'Proportion of Catch numbers at age',
       fill = "Age")+
  scale_y_continuous(expand = c(0, 0))
dev.off()

## weight in the stock by cohort
taf.png(paste0(prefix,"_wstock_cohort"))
west.by.cohort      <- as.data.frame(FLCohort(window(pil.all@stock.wt,2000,range(pil.all)["maxyear"])))
west.by.cohort      <- subset(west.by.cohort,!is.na(west.by.cohort$data))
west.by.cohort$year <- west.by.cohort$age + west.by.cohort$cohort
zoom(xyplot(data~year,data=west.by.cohort,
            groups=cohort,
            auto.key=list(space="right",points=FALSE,lines=TRUE,type="b"),
            type="b",
            xlab="Year",ylab="Weight in the stock (kg)",
            main=paste(pil.all@name,"Weight in the stock by cohort"),
            par.settings=list(superpose.symbol=list(pch=as.character(unique(west.by.cohort$cohort)%%10),cex=1.25)),
            panel=function(...) {
              panel.grid(h=-1,v=-1)
              panel.xyplot(...)
            }))
dev.off()

## Internal consistency MA_idx
taf.png(paste0(prefix,"_MA_idx_internal_consistency"))
zoom(plot(pil.all.tun[["MA_idx"]],type="internal"))
dev.off()

setwd('../..')

################################################################################
### ============================================================================
### Single-fleet
### ============================================================================
################################################################################

prefix <- 'assessment'
mkdir(file.path("report",prefix))
setwd(file.path("report",prefix))

##################################################
## assessment fit
##################################################
pdf(file.path(paste0(prefix,"_fit.pdf")))
residual.diagnostics(pil.sam)
dev.off()

##################################################
## stock trajectory
##################################################
taf.png(paste0(prefix,"_stock_trajectory"))
zoom(plot(pil.sam))
dev.off()

##################################################
## Catchability at age
##################################################
taf.png(paste0(prefix,"_survey_catchability"))
q <- catchabilities(pil.sam)
print(ggplot(q,aes(x=age,y=value))+
  geom_point()+
  geom_line()+
  geom_ribbon(aes(ymin=lbnd, ymax=ubnd),alpha=0.2)+
  facet_wrap(~fleet,scales = 'free_y'))+
  theme_bw()
dev.off()

##################################################
## observation variance by data source
##################################################
taf.png(paste0(prefix,"_observation_var_by_source"))
obv <- obs.var(pil.sam)
obv$str <- paste(obv$fleet,ifelse(is.na(obv$age),"",obv$age))
obv <- obv[order(obv$value),]
bp <- barplot(obv$value,ylab="Observation Variance",
              main="Observation variances by data source",col=factor(obv$fleet))
axis(1,at=bp,labels=obv$str,las=3,lty=0,mgp=c(0,0,0))
legend("topleft",levels(obv$fleet),pch=15,col=1:nlevels(obv$fleet),pt.cex=1.5)
dev.off()

##################################################
## variance vs uncertainty for each data source
#################################################
taf.png(paste0(prefix,"_variance_vs_uncertainty"))
plot(obv$value,obv$CV,xlab="Observation variance",ylab="CV of estimate",log="x",
     pch=16,col=obv$fleet,main="Observation variance vs uncertainty")
text(obv$value,obv$CV,obv$str,pos=4,cex=0.75,xpd=NA)
dev.off()

##################################################
## Selectivity pattern per 5 years
#################################################
taf.png(paste0(prefix,"_selectivity"))
sel.pat <- merge(f(pil.sam),fbar(pil.sam),
                 by="year",suffixes=c(".f",".fbar"))
sel.pat$sel <- sel.pat$value.f/sel.pat$value.fbar
sel.pat$age <- as.numeric(as.character(sel.pat$age))
zoom(xyplot(sel ~ age|sprintf("%i's",floor((year)/5)*5),sel.pat,
       groups=year,type="l",as.table=TRUE,
       scale=list(alternating=FALSE),
       main="Selectivity of the Fishery by Pentad",xlab="Age",ylab="F/Fbar"))
dev.off()

##################################################
## correlation matrix model params
#################################################
taf.png(paste0(prefix,"_cor_params"))
cor.plot(pil.sam)
dev.off()

##################################################
## catch residuals per year per age
#################################################
taf.png(paste0(prefix,"_catage_residuals"))
dat <- subset(residuals(pil.sam),fleet=="catch unique")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## nansen residuals per year per age
#################################################
taf.png(paste0(prefix,"_nansen_residuals"))
dat <- subset(residuals(pil.sam),fleet=="nansen")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## MA_idx residuals per year per age
#################################################
taf.png(paste0(prefix,"_MA_idx_residuals"))
dat <- subset(residuals(pil.sam),fleet=="MA_idx")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## IBTS0 residuals per year per age
#################################################
taf.png(paste0(prefix,"_MA_ssb_residuals"))
dat <- subset(residuals(pil.sam),fleet=="MA_ssb")
print(ggplot(dat,aes(x=year,y=age,size=10*abs(std.res),col=std.res>0))+
        geom_point(shape=21)+
        theme_bw()+
        theme(legend.position = "none")+
        scale_color_manual(values = c("red", "blue")))
dev.off()

##################################################
## F at age
#################################################
taf.png(paste0(prefix,"_fatage"))
zoom(timeseries(window(pil.all,2000,range(pil.all)["maxyear"]),slot="harvest"))
dev.off()

##################################################
## F/fbar at age
#################################################
pil.temp <- pil.all
pil.temp@harvest <- sweep(pil.all@harvest,2:6,fbar(pil.all),'/')

taf.png(paste0(prefix,"_ffbaratage"))
zoom(timeseries(window(pil.temp,2000,range(pil.all)["maxyear"]),slot="harvest"))
dev.off()

##################################################
## Recruitment vs SSB
#################################################

df.rec <- rec(pil.sam)
df.rec <- df.rec %>% select(c('year','value'))
colnames(df.rec) <- c('year','rec')
df.ssb <- ssb(pil.sam)
df.ssb <- df.ssb %>% select(c('year','value'))
colnames(df.ssb) <- c('year','ssb')

df.all <- left_join(df.rec,df.ssb,by=('year'))
df.all$yearclass <- df.all$year-1
df.all$ssbR <- df.all$ssb*100/df.all$rec
df.all$log.ssbR <- log10(df.all$ssbR)
df.all$RperS <- df.all$rec/df.all$ssb

taf.png(paste0(prefix,"_SR_periods"))
ggplot(data=df.all,aes(x=ssb,y=rec),aes(x=ssb,y=rec))+
  geom_line()+
  ylab('Recruitment 0-wr (thousands)')+
  xlab('Spawning Stock Biomass (t)')+
  theme_bw()
dev.off()

taf.png(paste0(prefix,"_recruit_spawners"))
ggplot()+
  geom_smooth(data=df.all,aes(x=yearclass,y=RperS),span = 0.1,col='black', linewidth=2)+
  geom_point(data=df.all,aes(x=yearclass,y=RperS), size=2)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_bw()+
  annotation_logticks()+
  xlab('Year class')+
  ylab('Recruits per spawner')
dev.off()

##################################################
## retro stock trajectory
#################################################
# taf.png(paste0(prefix,"_retrospective_stock"))
# zoom(plot(NSH.retro))
# dev.off()

##################################################
## model parameters retrospective
#################################################
# taf.png(paste0(prefix,"_retrospective_SAM_parameters"))
# zoom(retroParams(NSH.retro))
# dev.off()

##################################################
## fishing selectivity retrospective
#################################################
# taf.png(paste0(prefix,"_retrospective_SAM_parameters"))
# zoom(retroSelectivity(NSH.retro,(range(NSH)["maxyear"]-10):range(NSH)["maxyear"]))
# dev.off()

setwd('../..')

