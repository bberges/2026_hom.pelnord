config_sf_hom.treca_cecaf2025 <- function(  stk,
                                            tun,
                                            pg){
  sam.ctrl                                    <- FLSAM.control(stk,
                                                               tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  ssbRow     <- grep("idx.ssb",rownames(sam.ctrl@obs.vars))
  cpueRow     <- grep("idx.cpue",rownames(sam.ctrl@obs.vars))
  # recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(ssbRow,cpueRow)] <- 5
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(ssbRow,cpueRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  # sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # #All fishing mortality states are free except
  # #oldest ages to ensure stablity
  # sam.ctrl@states[catchRow,]                  <- c(0,0,0,0,1,1)#seq(dim(stk)[1])
  # 
  # # Group observation variances of catches to ensure stability
  # sam.ctrl@obs.vars[catchRow,]        <- c(0,1,2,2,3,3)
  # sam.ctrl@obs.vars[idxRow,ac(0:4)]   <- c(100,101,102,103,104)
  # # sam.ctrl@obs.vars[recRow,ac(0:1)]   <- c(200,201)
  # sam.ctrl@obs.vars[nansenRow,ac(0)]       <- c(300)
  # sam.ctrl@obs.vars[MA.ssbRow,ac(0)]       <- c(400)
  # 
  # #Catchabilities of the surveys. Set LAI all to 1 value, rest can be varied
  # sam.ctrl@catchabilities[idxRow,ac(0:4)]   <- c(100,101,102,102,102)
  # # sam.ctrl@catchabilities[recRow,ac(0:1)]   <- c(200,201)
  # sam.ctrl@catchabilities[nansenRow,ac(0)]       <- c(300)
  # sam.ctrl@catchabilities[MA.ssbRow,ac(0)]       <- c(400)
  # 
  # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,1:4]                   <- c(rep(101,4))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  # 
  # #Variance of F random walks
  # #sam.ctrl@f.vars[catchRow,]    <- c(0,0,1,1,2,2)
  # #correlation structure among ages in F random walks
  # sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "horse mackerel trecae cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}

config_sf_hom.trach_cecaf2025 <- function(  stk,
                                            tun,
                                            pg){
  sam.ctrl                                    <- FLSAM.control(stk,
                                                               tun)
  
  # catchRow                                    <- grep("catch",rownames(sam.ctrl@f.vars))
  catchRow      <- grep("catch",rownames(sam.ctrl@f.vars))
  ssbRow     <- grep("idx.ssb",rownames(sam.ctrl@obs.vars))
  cpueRow     <- grep("idx.cpue",rownames(sam.ctrl@obs.vars))
  # recRow        <- grep("rec_russ",rownames(sam.ctrl@obs.vars))
  
  # set up the ssb surveys
  biomassTreat <- sam.ctrl@biomassTreat[]
  biomassTreat[c(ssbRow,cpueRow)] <- 5
  sam.ctrl@biomassTreat[] <- biomassTreat
  
  fleets <- sam.ctrl@fleets[]
  fleets[c(ssbRow,cpueRow)] <- 3
  sam.ctrl@fleets[] <- fleets
  
  # #Variance in N random walks (set 1st one free is usually best)
  # sam.ctrl@logN.vars[]                        <- c(1,rep(2,length(1:pg)))
  # 
  # #All fishing mortality states are free except
  # #oldest ages to ensure stablity
  # sam.ctrl@states[catchRow,]                  <- c(0,0,0,0,1,1)#seq(dim(stk)[1])
  # 
  # # Group observation variances of catches to ensure stability
  # sam.ctrl@obs.vars[catchRow,]        <- c(0,1,2,2,3,3)
  # sam.ctrl@obs.vars[idxRow,ac(0:4)]   <- c(100,101,102,103,104)
  # # sam.ctrl@obs.vars[recRow,ac(0:1)]   <- c(200,201)
  # sam.ctrl@obs.vars[nansenRow,ac(0)]       <- c(300)
  # sam.ctrl@obs.vars[MA.ssbRow,ac(0)]       <- c(400)
  # 
  # #Catchabilities of the surveys. Set LAI all to 1 value, rest can be varied
  # sam.ctrl@catchabilities[idxRow,ac(0:4)]   <- c(100,101,102,102,102)
  # # sam.ctrl@catchabilities[recRow,ac(0:1)]   <- c(200,201)
  # sam.ctrl@catchabilities[nansenRow,ac(0)]       <- c(300)
  # sam.ctrl@catchabilities[MA.ssbRow,ac(0)]       <- c(400)
  # 
  # #Add correlation correction for Q3 survey, not for HERAS
  # idx                                       <- which(rownames(sam.ctrl@cor.obs)=="MA_idx")
  # sam.ctrl@cor.obs[idx,1:4]                   <- c(rep(101,4))
  # sam.ctrl@cor.obs.Flag[idx]                  <- as.factor("AR")
  # 
  # #Variance of F random walks
  # #sam.ctrl@f.vars[catchRow,]    <- c(0,0,1,1,2,2)
  # #correlation structure among ages in F random walks
  # sam.ctrl@cor.F                              <- 2
  #Finalise
  sam.ctrl@name                               <- "horse mackerel trachurus cecaf pelnord"
  sam.ctrl                                    <- update(sam.ctrl)
  
  sam.ctrl@residuals <- FALSE
  
  return(sam.ctrl)
}
