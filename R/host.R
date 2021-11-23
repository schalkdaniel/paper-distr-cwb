library(R6)

#' Aggregate static base learner information
#'
#' @param bl_parts (`list()`)\cr
#'   List containing all base learner parts.
#' @return
#'   Aggregated parts by base learner.
#' @export
aggrParts = function(bl_parts) {
  blns = names(bl_parts[[1]])
  aggr = lapply(blns, function(bln) {
    blp = lapply(bl_parts, function(r) r[[bln]] )
    XtX = Reduce("+", lapply(blp, function(b) b$XtX))
    K   = blp[[1]]$K
    df  = blp[[1]]$df
    pen = compboostSplines::demmlerReinsch(as.matrix(XtX), K, df)

    out = list(XtX = XtX, K = K, df = df, pen = pen)
    return(out)
  })
  names(aggr) = blns
  return(aggr)
}

#' Aggregate X^Tu from each base learner
#'
#' @param receive (`list()`)\cr
#'   List containing all X^Tu matrices.
#' @return
#'   Aggregated X^Tu by base learner.
#' @export
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

#' Calculate parameter extimates
#'
#' @param aggr_parts (`list()`)\cr
#'   Result from `aggrParts()`.
#' @param aggr_xtu (`list()`)\cr
#'   Result from `aggrXtu()`.
#' @return
#'   Estimated parameter vectors for all base learner
#'   in the lists.
#' @export
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

#' @title Host class
#'
#' @description
#' This class simulates the data hold by the host, such as the
#' site connections, the host CWB model, or the  structure of
#' the base learners.
#'
#' @export
Host = R6Class("Host",
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param cwb (`[CWB]`)\cr
    #'   CWB object containing the fitted model.
    #' @param bls (`list([Baselearner])`)\cr
    #'   List of all base learners used to train the model.
    #' @param sites (`list([Site])`)\cr
    #'   List of [Site]s containing site specific data.
    initialize = function(cwb, bls, sites) {
      private$p_cwb   = cwb
      private$p_bls   = bls
      private$p_sites = sites
      private$p_blns  = sapply(private$p_bls, function(bl) bl$getName())
      names(private$p_bls) = private$p_blns
    },

    #' @description
    #' Initialize all base learners in `bls`.
    initializeBaselearner = function() {
      init = lapply(sites, function(site) site$communicateInit())
      init_aggr = lapply(private$p_blns, function(bln) {
        sinit = lapply(init, function(i) i[[bln]])
        out = list(dummy = 0)
        if (grepl("spline", bln)) {
          smin = min(sapply(sinit, function(s) s$min))
          smax = max(sapply(sinit, function(s) s$max))
          knots = private$p_bls[[bln]]$createKnots(c(smin, smax))
          out = list(knots = knots)
        }
        if (grepl("ridge", bln)) {
          out = list(classes = unname(unlist(sinit)))
        }
        return(out)
      })
      names(init_aggr) = private$p_blns
      private$p_init_aggr = init_aggr

      # Initialize base learner:
      nuisance = lapply(private$p_sites, function(s) s$initBaselearner(private$p_init_aggr))
      bl_parts = lapply(private$p_sites, function(s) s$communicateFittingParts())
      private$p_bl_parts = aggrParts(bl_parts)

      return(invisible(NULL))
    },

    #' @description
    #' Initialize the prediction in CWB.
    initializePrediction = function(constant_aggregator = function(y) mean(y)) {
      ss = sapply(private$p_sites, function(s) s$getAggregatedInitialization(constant_aggregator))
      sw = sapply(private$p_sites, function(s) s$nrow())
      private$p_offset = weighted.mean(ss, sw)

      nuisance = lapply(private$p_sites, function(s) s$initPrediction(private$p_offset))
      return(invisible(NULL))
    },

    #' @description
    #' Initialize the static parts of CWB.
    initializeCWB = function() {
      bl_parts = lapply(private$p_sites, function(s) s$communicateFittingParts())
      private$p_bl_parts = aggrParts(bl_parts)
      return(invisible(NULL))
    },

    #' @description
    #' Function to conduct one update step for CWB. These steps
    #' are sharing X^Tu from all sites, calculate parameter estimates,
    #' find the best base learner, and update CWB.
    #'
    #' @param trace (`integer(1L)`)\cr
    #'   Number after how many iterations the trace should be printed.
    updateCWB = function(trace = 1L) {
      checkmate::assertIntegerish(trace, len = 1L)

      m = length(private$p_cwb$getBlTrace())

      if (! m %% trace)
        cat("> ", m, ":\n", sep = "")

      # Calculate params
      if (! m %% trace)
        cat("  >> Receive X^Tu from all sites for all base learners\n")

      xtus     = lapply(private$p_sites, function(s) s$communicateXtu())
      aggr_xtu = aggrXtu(xtus)

      if (! m %% trace)
        cat("  >> Aggregate to get global parameter estimates for all base learners\n")

      blparams = estimateParam(private$p_bl_parts, aggr_xtu)

      # Get SSEs:
      if (! m %% trace)
        cat("  >> Send best parameters to sites to get SSEs\n")

      idx_bls_included = sapply(private$p_bls, function(bl) bl$isIncluded())
      sses = sapply(private$p_blns[idx_bls_included], function(bln) {
        bp = blparams[[bln]]
        sum(sapply(private$p_sites, function(s) {
          s$ssePrUpdate(bln, bp)
        }))
      })
      names(sses) = private$p_blns[idx_bls_included]

      bl_best = which.min(sses)
      if (! m %% trace)
        cat("  >> Found best base learner:", names(bl_best), "\n")

      # Site models:
      if (! m %% trace)
        cat("  >> Send parameter of best base learner to site to update CWB\n")

      nuisance = lapply(private$p_sites, function(s) {
        s$update(names(bl_best), blparams[[names(bl_best)]])
      })
      # Host model:
      private$p_cwb$update(names(bl_best), blparams[[names(bl_best)]])

      sw = sapply(private$p_sites, function(s) s$nrow())
      private$p_risk = c(private$p_risk, weighted.mean(sapply(private$p_sites, function(s) s$communicateRisk()), sw))

      if (! m %% trace)
        cat("  >> Risk:", tail(private$p_risk, 1), "\n")

      return(invisible(NULL))
    }
  ),
  private = list(

    #' @field p_cwb (`[CWB]`)\cr
    #'   The CWB object.
    p_cwb = NULL,

    #' @field p_bls (`list([Baselearner])`)\cr
    #'   List of [Baselearner]s.
    p_bls = NULL,

    #' @field p_blns (`character()`)\cr
    #'   Names of the base learners in `p_bls`.
    p_blns = NULL,

    #' @field p_sites (`list([Site])`)\cr
    #'   List of [Site]s.
    p_sites = NULL,

    #' @field p_offset (`numeric(1L)`)\cr
    #'   Offset value.
    p_offset = NULL,

    #' @field p_init_aggr (`list()`)\cr
    #'   Static initial values (such as knots or classes).
    p_init_aggr = NULL,

    #' @field p_bl_parts (`list()`)\cr
    #'   List containing the static base learner parts (such
    #'   as X^TX, the penalty, or the penalty matrix).
    p_bl_parts  = NULL,

    #' @field p_risk (`numeric()`)\cr
    #'   Vector of the risk after each iteration.
    p_risk = NULL
  ),
  active = list(

    #' @field cwb (`[CWB]`)\cr
    #'   The CWB object.
    cwb = function(x) {
      if (! missing(x)) stop("Cannot set CWB.")
      return(private$p_cwb)
    },

    #' @field bls (`list([Baselearner])`)\cr
    #'   List of [Baselearner]s.
    bls = function(x) {
      if (! missing(x)) stop("Cannot set list of base learners.")
      return(private$p_bls)
    },

    #' @field sites (`list([Site])`)\cr
    #'   List of [Site]s.
    sites = function(x) {
      if (! missing(x)) stop("Cannot set sites.")
      return(private$p_sites)
    },

    #' @field offset (`numeric(1L)`)\cr
    #'   Offset value.
    offset = function(x) {
      if (! missing(x)) stop("Cannot set offset")
      return(private$p_offset)
    },

    #' @field init_aggr (`list()`)\cr
    #'   Static initial values (such as knots or classes).
    init_aggr = function(x) {
      if (! missing(x)) stop("Cannot set initialization.")
      return(private$p_init_aggr)
    },

    #' @field bl_parts (`list()`)\cr
    #'   List containing the static base learner parts (such
    #'   as X^TX, the penalty, or the penalty matrix).
    bl_parts = function(x) {
      if (! missing(x)) stop("Cannot set base learner parts.")
      return(private$p_bl_parts)
    },

    #' @field risk (`numeric()`)\cr
    #'   Vector of the risk after each iteration.
    risk = function(x) {
      if (! missing(x)) stop("Cannot set risk.")
      return(private$p_risk)
    }
  )
)

