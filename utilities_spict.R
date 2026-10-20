
# buildTMBinputSPICT {{{

#' Build SPiCT Input List from FLR Objects
#'
#' Constructs the input list required by \code{spict::fit.spict()} from an
#' \code{FLQuant} of catches and a list of survey index \code{FLQuant}s.
#'
#' @param catch An \code{FLQuant} containing observed catch time series.
#' @param indices A named list of \code{FLQuant} objects, each containing a
#'   survey abundance index time series.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{\code{obsC}}{Numeric vector of catch observations.}
#'     \item{\code{timeC}}{Numeric vector of years for catch observations.}
#'     \item{\code{obsI}}{List of numeric vectors of index observations,
#'       one element per survey.}
#'     \item{\code{timeI}}{List of numeric vectors of years for index
#'       observations, one element per survey.}
#'   }
#'
#' @seealso \code{\link{spict.settings}}, \code{\link{spict.sa}}
#'
#' @keywords internal

buildTMBinputSPICT <- function(catch, indices) {
  
  # PARSE and CHECK arguments
  ny <- dim(catch)[2]
  
  # C_hist
  
  dataC <- as.data.frame(catch)
  obsC <- dataC$data
  timeC <- dataC$year
  
  # obsI
  
  dataI <- lapply(indices, function(x) as.data.frame(x))
  obsI <- lapply(dataI, '[[', 'data')
  timeI <- lapply(dataI, '[[', 'year')
  
  # BUILD input list
  inp <- list(obsC=obsC, timeC=timeC, obsI=obsI, timeI=timeI)
  
  return(inp)
}
# }}}

# spict2ind {{{

#' Extract SPiCT Indicators as FLQuants
#'
#' Assembles a set of stock status indicators from a fitted SPiCT model
#' into a named \code{FLQuants} object for use within FLR workflows.
#'
#' @param fit A fitted SPiCT model object, as returned by
#'   \code{spict::fit.spict()}.
#'
#' @return A named \code{FLQuants} with the following elements:
#'   \describe{
#'     \item{\code{F}}{Fishing mortality (\code{logFnotS}).}
#'     \item{\code{Fstatus}}{F relative to Fmsy (\code{logFFmsynotS}).}
#'     \item{\code{B}}{Biomass (\code{logB}).}
#'     \item{\code{Bstatus}}{B relative to Bmsy (\code{logBBmsy}).}
#'     \item{\code{Bdepletion}}{Biomass depletion, computed as B / K.}
#'     \item{\code{C}}{Predicted catch (\code{logCpred}).}
#'   }
#'
#' @seealso \code{\link{spict2FLQuant}}, \code{\link{spict.sa}}
#'
#' @keywords internal

spict2ind <- function(fit) {
  
  out <- list(
    F = spict2FLQuant(fit, val="logFnotS"),
    Fstatus = spict2FLQuant(fit, val="logFFmsynotS"),
    B = spict2FLQuant(fit, val="logB"),
    Bstatus = spict2FLQuant(fit, val="logBBmsy"),
    Bdepletion = spict2FLQuant(fit, val="logB") / fit$value[['K']],
    C = spict2FLQuant(fit, val="logCpred"))
  
  return(FLQuants(out))
}
# }}}

# spict2FLQuant {{{

#' Extract SPiCT Model Output as FLQuant
#'
#' Converts a named quantity from a fitted SPiCT model object into an
#' \code{FLQuant}, handling seasonal structure and time indexing automatically.
#'
#' @param x A fitted SPiCT model object, as returned by
#'   \code{spict::fit.spict()}.
#' @param val Character string specifying the quantity to extract. One of
#'   \code{"logB"} (biomass), \code{"logFnotS"} (fishing mortality),
#'   \code{"logCpred"} (predicted catch), \code{"logBBmsy"} (B/Bmsy status),
#'   or \code{"logFFmsynotS"} (F/Fmsy status). Defaults to \code{"logB"}.
#'
#' @return An \code{FLQuant} on the natural scale (i.e. exp-transformed),
#'   with years spanning \code{x$inp$timerange} and seasonal structure
#'   determined by \code{x$inp$dteuler}. For \code{"logCpred"}, seasonal
#'   values are summed to annual totals via \code{seasonSums()}; for all
#'   other quantities the final season only is returned, with
#'   \code{dimnames$season} set to \code{"all"}.
#'
#' @details
#' Biomass (\code{"logB"}) is extracted from \code{x$par.random}; all other
#' quantities are taken from \code{x$value}. The number of seasons per year
#' is derived from \code{round(1 / x$inp$dteuler)}, except for
#' \code{"logCpred"} where \code{x$inp$dtc[1]} is used.
#'
#' @seealso \code{\link{spict2ind}}, \code{\link{spict.sa}}
#'
#' @keywords internal

spict2FLQuant <- function(x,
                          val=c("logB", "logFnotS", "logCpred", "logBBmsy", "logFFmsynotS")) {
  
  val <- match.arg(val)
  
  # logB
  if(val == "logB") {
    vec = x$par.random
  } else {
    vec = x$value
  }
  
  quant <- an(exp(vec[which(names(vec)==val)]))
  
  if(val == "logCpred") {
    season = x$inp$dtc[1]
  } else {
    season = round(1/x$inp$dteuler)
  }
  
  year <- c(rep(seq(x$inp$timerange[1], x$inp$timerange[2]), each=season))
  quant <- quant[1:length(year)]
  
  dat <- data.frame(age=1,
                    year=c(rep(seq(x$inp$timerange[1], x$inp$timerange[2]),each=season)),
                    season=c(rep(seq(season), x$inp$timerange[2] - x$inp$timerange[1]+1)),
                    data=quant)
  
  if(val == "logCpred") {
    out <- seasonSums(as.FLQuant(dat))
  } else {
    out <- as.FLQuant(dat)[,,,season] 
    dimnames(out)$season <- "all"
  }
  return(out)
}
# }}}

# spict.settings {{{

#' Configure Priors and Settings for a SPiCT Model Input
#'
#' Applies a standard set of prior distributions and model settings to a
#' SPiCT input list, as returned by \code{\link{buildTMBinputSPICT}}. All
#' priors follow the SPiCT convention of \code{c(log-mean, log-sd, active)},
#' where \code{active = 1} enables the prior and \code{active = 0} disables it.
#' Setting any argument to \code{NULL} skips that prior.
#'
#' @param inp A SPiCT input list, as returned by
#'   \code{\link{buildTMBinputSPICT}}.
#' @param r.pr Numeric vector \code{c(mean, sd, active)} for the prior on the
#'   intrinsic growth rate \emph{r}. Defaults to \code{c(0.4583, 0.2553, 1)}.
#' @param bk.pr Numeric vector \code{c(mean, sd, active)} for the prior on the
#'   initial depletion (B0/K). Defaults to \code{c(0.5, 0.3, 1)}.
#' @param shape.pr Numeric vector \code{c(mean, sd, active)} for the prior on
#'   the Pella-Tomlinson shape parameter \emph{n}. Defaults to
#'   \code{c(2, 0.001, 1)}, effectively fixing \emph{n} near the Schaefer
#'   value.
#' @param oe Numeric vector \code{c(mean, sd, active)} for the prior on
#'   observation error (survey index SD). Applied identically to all surveys.
#'   Defaults to \code{c(0.3, 0.5, 1)}.
#' @param pe Numeric vector \code{c(mean, sd, active)} for the prior on
#'   process error (biomass SD). Defaults to \code{c(0.2, 0.5, 1)}.
#' @param fdevs Numeric vector \code{c(mean, sd, active)} for the prior on
#'   F deviations SD. Defaults to \code{c(4, 0.5, 1)}.
#' @param ce Numeric vector \code{c(mean, sd, active)} for the prior on catch
#'   observation error SD. Defaults to \code{c(0.05, 0.001, 1)}.
#' @param q.pr Numeric vector \code{c(mean, sd, active)} for the prior on the
#'   survey catchability \emph{q}. Applied identically to all surveys.
#'   Defaults to \code{c(1, 0.5, 1)}.
#' @param timing Numeric scalar (or vector, one per survey) giving the fraction
#'   of the year at which surveys are conducted, added to the integer survey
#'   year. Defaults to \code{0.625} (mid-August).
#' @param dteuler Numeric scalar controlling the Euler discretisation time
#'   step. Defaults to \code{0.25} (quarterly).
#'
#' @return The modified SPiCT input list \code{inp} with priors and settings
#'   applied.
#'
#' @details
#' The \code{logalpha} (ratio of process to observation error) and
#' \code{logbeta} (ratio of catch error to F deviations) priors are
#' unconditionally deactivated, as their estimation is known to be
#' unreliable without strong external information.
#'
#' @seealso \code{\link{buildTMBinputSPICT}}, \code{\link{spict.sa}}
#'
#' @keywords internal

spict.settings = function(inp, r.pr=c(0.4583,  0.2553, 1), bk.pr=c(0.5,0.3,1),
                          shape.pr=c(2, 0.001, 1), oe=c(0.3, 0.5, 1), pe= c(0.2,0.5,1),
                          fdevs=c(4, 0.5, 1), ce=c(0.05, 0.001, 1),q.pr= c(1,0.5,1),timing=0.625, dteuler=0.25) {
  
  # ASSIGN arguments to input
  inp$dteuler=dteuler
  
  #ASSIGN survey timing
  inp$timeI = Map(function(x,y){x+y}, y=as.list(timing),x=inp$timeI)
  
  # DEACTIVATE ratio of proc to obs
  inp$priors$logalpha <- c(0,0,0) 
  
  # DEACTIVATE catch to f devs
  inp$priors$logbeta <- c(0,0,0) 
  
  # SET bk prior'
  if(!is.null(bk.pr))
    inp$priors$logbkfrac <- c(log(bk.pr[1])-bk.pr[2]^2/2,bk.pr[2],bk.pr[3])
  
  # SET obs error
  if(!is.null(oe)){
    inp$priors$logsdi = list()
    for(i in 1:length(inp$obsI)){
      inp$priors$logsdi[[i]] = c(log(oe[1]), oe[2], oe[3])  
    }
  }
  # SET process error   
  if(!is.null(pe))
    inp$priors$logsdb <- c(log(pe[1])-0.5*pe[2]^2, pe[2], pe[3]) 
  
  # SET F devs: should be higher
  if(!is.null(fdevs))
    inp$priors$logsdf <- c(log(fdevs[1])-0.5*fdevs[2]^2, fdevs[[2]], fdevs[3])
  
  # SET catch error
  if(!is.null(ce))
    inp$priors$logsdc <- c(log(ce)[1]-0.5*ce[2]^2, ce[2], ce[3])
  
  # REDUCE shape CV
  if(!is.null(shape.pr))
    inp$priors$logn <- c(log(shape.pr[1]),shape.pr[2],shape.pr[3])
  
  # SET r prior
  if(!is.null(r.pr))
    inp$priors$logr <- c(log(r.pr[1])-0.5*r.pr[2]^2,r.pr[2],r.pr[3])
  
  if(!is.null(q.pr)){
    inp$priors$logq = list()
    for(i in 1:length(inp$obsI)){
      inp$priors$logq[[i]] = c(log(q.pr[1])-0.5*q.pr[2]^2,q.pr[2],q.pr[3]) 
    }}
  
  
  return(inp)
}
# }}}

# spict.settings.CJM {{{ ----

spict.settings.CJM = function(inp){
  # ASSIGN arguments to input
  inp$dteuler=0.25
  
  # DEACTIVATE ratio of proc to obs
  inp$priors$logalpha <- c(0,0,0)
  
  # DEACTIVATE catch to f devs
  inp$priors$logbeta <- c(0,0,0)
  
  # SET bk prior'
  bk.pr=c(1,0.35,1)
  inp$priors$logbkfrac <- c(log(bk.pr[1])-bk.pr[2]^2/2,bk.pr[2],bk.pr[3])
  
  # SET obs error [#><>these error are way too high]
  oe=c(0.2, 0.3, 1)
  inp$priors$logsdi = list()
  for(i in 1:length(inp$obsI)){
    inp$priors$logsdi[[i]] = c(log(oe[1]), oe[2], oe[3])  
  }
  
  # c("Chile_AcousCS","Chile_AcousN","Offshore_CPUE", "Chile_CPUE", "Peru_CPUE")
  
  # SET process error # Too high changed to 0.15
  pe= c(0.10,0.3,1)
  inp$priors$logsdb 
  
  # SET F devs: #><> Indeed  should be higher 
  fdevs=c(0.5, 0.1, 1)
  inp$priors$logsdf <- c(log(fdevs[1])-0.5*fdevs[2]^2, fdevs[[2]], fdevs[3])
  
  # SET catch error
  ce=c(0.05, 0.001, 1)
  inp$priors$logsdc <- c(log(ce)[1]-0.5*ce[2]^2, ce[2], ce[3])
  
  # REDUCE shape CV
  
  #><> SCHAEFER
  shape.pr=c(2, 0.001, 1)
  inp$priors$logn <- c(log(shape.pr[1]),shape.pr[2],shape.pr[3])
  
  # SET r prior
  # For Schaefer use Fmsy = r/2
  r.pr=c(0.29,  0.3, 1)
  inp$priors$logr <- c(log(r.pr[1])-0.5*r.pr[2]^2,r.pr[2],r.pr[3])
  
  # SET q prior
  inp$priors$logq = NULL 
  
  return(inp)
}

# }}}

# spict.sa {{{

#' SPiCT Stock Assessment Module for FLR/mse
#'
#' An \code{est} module for use in \code{\link{mpCtrl}} that fits the Surplus
#' Production in Continuous Time (SPiCT) model to catch and survey index data,
#' iterating over Monte Carlo replicates. Returns stock indicators, reference
#' points, and convergence tracking for use in the MSE feedback loop.
#'
#' @param stk An \code{FLStock} object containing the perceived stock, as
#'   passed by the MSE framework.
#' @param idx An \code{FLIndices} object containing the survey indices.
#' @param indices Character vector of index names within \code{idx} to use.
#'   Defaults to all indices in \code{idx}.
#' @param ... Additional named arguments passed directly to the SPiCT input
#'   list, overriding defaults set by \code{\link{spict.settings}}.
#' @param args A list of MSE run arguments supplied by the framework, including
#'   at minimum \code{ay} (assessment year) and \code{it} (number of
#'   iterations).
#' @param tracking An \code{FLQuant}-based tracking object for recording
#'   diagnostic quantities across the MSE run.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{\code{stk}}{The input \code{FLStock} object, unchanged.}
#'     \item{\code{ind}}{A named \code{FLQuants} of stock indicators (B, F,
#'       Bstatus, Fstatus, Bdepletion, C) assembled across all iterations and
#'       windowed to the assessment year.}
#'     \item{\code{tracking}}{Updated tracking object with convergence flags
#'       (\code{conv.est}), biomass estimates (\code{B.est}), and biomass
#'       status estimates (\code{Bstatus.est}) for the assessment year.}
#'     \item{\code{refpts}}{An \code{FLPar} containing per-iteration estimates
#'       of MSY, Fmsy, Bmsy, and K.}
#'   }
#'
#' @details
#' The function loops over iterations (\code{args$it}), fitting SPiCT
#' independently to each replicate using total catch (summed over units and
#' areas) and the corresponding index slices. Model settings and priors are
#' applied via \code{\link{spict.settings}}, with any additional arguments in
#' \code{...} appended afterwards.
#'
#' Convergence is recorded in \code{tracking} using
#' \code{fit$opt$convergence}: \code{0} indicates successful convergence;
#' \code{3} is assigned when \code{fit.spict()} throws an error, in which case
#' reference points are set to zero.
#'
#' @seealso \code{\link{mp}}, \code{\link{mpCtrl}}, \code{\link{mseCtrl}},
#'   \code{\link{spict2ind}}, \code{\link{spict.settings}}
#'
#' @references
#' Pedersen, M.W. and Berg, C.W. (2017). A stochastic surplus production model
#' in continuous time. \emph{Fish and Fisheries}, 18(2): 226--243.
#' \doi{10.1111/faf.12174}
#'
#' @examples
#' \dontrun{
#' data(sol274)
#' control <- mpCtrl(list(
#'   est = mseCtrl(method = spict.sa,
#'    args = list(indices = "BTS")),
#'   hcr = mseCtrl(method = catchSSB.hcr,
#'     args = list(MSY = 14000))
#' ))
#' res <- mp(om, oem = oem, ctrl = control, args = list(iy = 2021))
#' }

spict.sa <- function(stk, idx, indices=names(idx), ..., args, tracking) {
  
  # EXTRACT inputs
  catch <- catch(stk) # + min(c(catch(stk))[c(catch(stk)) > 0]) * 1
  indices <- lapply(idx[indices], index)
  
  # safeguard for convergence when indices and catch get very low
  for(i in 1:length(indices)){
    indices[i][[1]][indices[i][[1]] < 1e-4] <- 0
  }
  catch[catch < 1e-4] <- 0
  
  # OUTPUT
  empty <- catch %=% 0
  out <- vector(mode="list", length=args$it)
  
  # LOOP over iters
  for (i in seq(args$it)) {
    # cat(".")
    
    # EXTRACT single iter
    cab <- unitSums(areaSums(iter(catch, i)))
    cab <- window(cab,start=1984)
    ind <- iter(indices, i)
    #ind$Chile_AcousN <- window(ind$Chile_AcousN,start=2006)
    
    # BUILD input list
    inp <- buildTMBinputSPICT(cab, ind)
    
    # ADD settings to input
    inp <- spict.settings.CJM(inp)
    
    #save(inp,file=file.path('./model/spict.last.RData'))
    # ADD provided settinga
    spargs <- list(...)
    
    if(length(spargs) > 0) {
      inp[names(spargs)] <- spargs
    }
    
    # FIT and catch nmlinb errors
    fit <- fit.spict(inp,verbose = F)
    
    if(fit$opt$convergence != 0){
      inp$dteuler <- 0.5
      track(tracking, "conv.est", ac(args$ay)) <- fit$opt$convergence
      fit <- fit.spict(inp,verbose = F)
    }else{
      track(tracking, "conv.est", ac(args$ay)) <- fit$opt$convergence
    }
    # convergence flag
    conv <- fit$opt$convergence
    
    # GET refpts
    rps <- unlist(c(fit$report[c("MSY", "Fmsy", "Bmsy")],
                    K=fit$value[['K']]))
    
    # EXTRACT output indicators
    ind <- spict2ind(fit)
    ind <- window(ind,end=args$ay)
    
    out[[i]] <- list(ind=ind, conv=conv, rps=rps, fit=fit, inp=inp)
    
  }
  
  # TRACK convergence
  # mse::track(tracking, "conv.est", ac(args$ay)) <- unlist(lapply(out, '[[', 'conv'))
  
  # COMBINE results
  inds <- lapply(out, "[[", "ind")
  
  # ASSSEMBLE ind ...
  ind <- lapply(setNames(nm=names(inds[[1]])), function(i) {
    window(Reduce(FLCore::combine, lapply(inds, "[[", i)), end=args$ay)
  })
  
    # TRACK B estimates
  track(tracking, "B.est.spict", ac(args$ay)) <- ind[['B']][, ac(args$ay)]
  
  # TRACK B status
  track(tracking, "Bstatus.est.spict", ac(args$ay)) <- ind[['Bstatus']][, ac(args$ay)]
  
  # ... and refpts
  rps <- Reduce(FLCore::combine, lapply(out, function(x) FLPar(x$rps)))
  
  # LISTS of fits and inputs
  fit <- lapply(out, '[[', 'fit')
  inp <- lapply(out, '[[', 'inp')
  
  # TODO: ASSEMBLE FLStockR
  
  return(list(stk=stk, ind=ind, tracking=tracking, refpts=rps, fit=fit, inp=inp))
}
# }}}

# vb(FLombf) {{{
