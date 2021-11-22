library(R6)

aggrParts = function(bl_parts) {
  blns = names(bl_parts[[1]])
  aggr = lapply(blns, function(bln) {
    blp = lapply(bl_parts, function(r) r[[bln]] )
    XtX = Reduce("+", lapply(blp, function(b) b$XtX))
    K   = blp[[1]]$K
    df  = blp[[1]]$df
    pen = compboostSplines::demmlerReinsch(XtX, K, df)

    out = list(XtX = XtX, K = K, df = df, pen = pen)
    return(out)
  })
  names(aggr) = blns
  return(aggr)
}
aggrXtu = function(receive) {
  blns = names(receive[[1]])
  Xtu = lapply(blns, function(bln) {
    blp = lapply(receive, function(r) r[[bln]])
    Xtu = Reduce("+", lapply(blp, function(b) b$Xtu))
    return(list(Xtu = Xtu))
  })
  names(Xtu) = blns
  return(Xtu)
}
estimateParam = function(aggr_parts, aggr_xtu) {
  blns = names(aggr_parts)
  params = lapply(blns, function(bln) {
    XtX = aggr_parts[[bln]]$XtX
    K   = aggr_parts[[bln]]$pen * aggr_parts[[bln]]$K
    Xtu = aggr_xtu[[bln]]$Xtu
    return(solve(XtX + K) %*% Xtu)
  })
  names(params) = blns
  return(params)
}


Host = R6Class("Host",
  public = list(
    initialize = function(cwb, bls, sites) {
      private$p_cwb = cwb
      private$p_bls = bls
      private$p_sites = sites
      private$p_blns  = sapply(private$p_bls, function(bl) bl$getName())
    },
    initializeBaselearner = function() {
      init = lapply(sites, function(site) site$communicateInit())
      init_aggr = lapply(private$p_blns, function(bln) {
        sinit = lapply(init, function(i) i[[bln]])
        if (grepl("spline", bln)) {
          smin = min(sapply(sinit, function(s) s$min))
          smax = max(sapply(sinit, function(s) s$max))
          out = compboostSplines::createKnots(c(smin, smax), n_knots = 10L, degree = 2)
        }
        out = list(knots = out)
        return(out)
      })
      names(init_aggr) = private$p_blns
      private$p_init_aggr = init_aggr

      # Initialize base learner:
      nuisance = lapply(private$p_sites, function(s) s$initBaselearner(private$p_init_aggr))
      bl_parts = lapply(private$p_sites, function(s) s$communicateFittingParts())
      private$p_bl_parts = aggrParts(bl_parts)
    },
    initializePrediction = function(constant_aggregator = function(y) mean(y)) {
      ss = sapply(private$p_sites, function(s) s$getAggregatedInitialization(constant_aggregator))
      sw = sapply(private$p_sites, function(s) s$nrow())
      private$p_offset = weighted.mean(ss, sw)

      nuisance = lapply(sites, function(s) s$initPrediction(private$p_offset))
    },
    initializeCWB = function() {
      bl_parts = lapply(private$p_sites, function(s) s$communicateFittingParts())
      private$p_bl_parts = aggrParts(bl_parts)
    },
    updateCWB = function(m = NULL) {

      # Calculate params
      xtus     = lapply(private$p_sites, function(s) s$communicateXtu())
      aggr_xtu = aggrXtu(xtus)
      blparams = estimateParam(private$p_bl_parts, aggr_xtu)

      # Get SSEs:
      sses = sapply(private$p_blns, function(bln) {
        bp = blparams[[bln]]
        sum(sapply(private$p_sites, function(s) {
          s$ssePrUpdate(bln, bp)
        }))
      })
      names(sses) = private$p_blns

      bl_best = which.min(sses)

      # Site models:
      nuisance = lapply(private$p_sites, function(s) {
        s$update(names(bl_best), blparams[[names(bl_best)]])
      })
      # Host model:
      private$p_cwb$update(names(bl_best), blparams[[names(bl_best)]])

      sw = sapply(private$p_sites, function(s) s$nrow())
      private$p_risk = c(private$p_risk, weighted.mean(sapply(private$p_sites, function(s) s$communicateRisk()), sw))
      cat(m, ": ", tail(private$p_risk, 1), "\n", sep = "")
    }
  ),
  private = list(
    p_cwb = NULL,
    p_bls = NULL,
    p_blns = NULL,
    p_sites = NULL,
    p_offset = NULL,
    p_init_aggr = NULL,
    p_bl_parts  = NULL,
    p_risk = NULL
  ),
  active = list(
    cwb = function(x) {
      if (! missing(x)) stop("Cannot set CWB.")
      return(private$p_cwb)
    },
    bls = function(x) {
      if (! missing(x)) stop("Cannot set list of base learners.")
      return(private$p_bls)
    },
    sites = function(x) {
      if (! missing(x)) stop("Cannot set sites.")
      return(private$p_sites)
    },
    offset = function(x) {
      if (! missing(x)) stop("Cannot set offset")
      return(private$p_offset)
    },
    init_aggr = function(x) {
      if (! missing(x)) stop("Cannot set initialization.")
      return(private$p_init_aggr)
    },
    bl_parts = function(x) {
      if (! missing(x)) stop("Cannot set base learner parts.")
      return(private$p_bl_parts)
    },
    risk = function(x) {
      if (! missing(x)) stop("Cannot set risk.")
      return(private$p_risk)
    }
  )
)

