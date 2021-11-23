library(R6)

#' @title Site class
#'
#' @description
#' This class simulates the sites or servers with individual data.
#' It defines what gets shared and received by the host.
#'
#' @export
Site = R6Class("Site",
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param name (`character(1L)`)\cr
    #'   Name of the site.
    #' @param data (`data.frame()`)\cr
    #'   Data frame exclusively hold by the server.
    #' @param cwb (`[CWB]`)\cr
    #'   Object of R6 class CWB.
    #' @param bls (`list()`)\cr
    #'   List of object of R6 class base learner.
    initialize = function(name, data, cwb, bls) {
      checkmate::assertCharacter(name, len = 1L)
      checkmate::assertDataFrame(data)
      checkmate::assertR6(cwb, "CWB")
      nuisance = lapply(bls, function(bl) checkmate::assertR6(bl, "Baselearner"))

      private$p_name = name
      private$p_data = data
      private$p_cwb  = cwb$clone()

      private$p_bls = lapply(bls, function(bl) bl$clone(deep = TRUE))
      names(private$p_bls) = vapply(private$p_bls, FUN.VALUE = character(1L),
        FUN = function(bl) bl$getName())
      return(invisible(NULL))
    },

    #' @description
    #' Get the name of the base learners.
    #'
    #' @return
    #'   Returns the base learner name as character value.
    baselearnerNames = function() {
      return(names(private$p_bls))
    },

    #' @description
    #' Initialize all base learners.
    #'
    #' @param ll_init (`list()`)\cr
    #'   List containing the information to initialize the
    #'   base learners (such as knots or classes).
    initBaselearner = function(ll_init) {
      if (private$p_is_initialized_bl)
        stop("Baselearner already initialized. Call `$unlockBaselearner()` and re initialize.")

      ninits = names(ll_init)
      blc = character()
      for (ni in ninits) {
        if (! private$p_bls[[ni]]$isCombination())
          do.call(private$p_bls[[ni]]$initDesign, c(list(dat = private$p_data), ll_init[[ni]]))
        else
          blc = c(blc, private$p_bls[[ni]]$getName())
      }
      for (bln in blc) {
        blc = private$p_bls[[bln]]
        do.call(private$p_bls[[ni]]$initDesign,
          c(list(dat = private$p_data,
            bl1 = private$p_bls[[blc$getBl1Name()]],
            bl2 = private$p_bls[[blc$getBl2Name()]])))
      }
      if (any(! sapply(private$p_bls, function(bl) bl$isLocked())))
        warning("Not all base learners were initialized.")
      private$p_is_initialized_bl = TRUE
      return(invisible(NULL))
    },

    #' @description
    #' Make the base learner "editable", meaning that it can be
    #' initialized.
    unlockBaselearner = function() {
      private$p_is_initialized_bl = FALSE
      return(invisible(NULL))
    },

    #' @description
    #' Communicate (get) the risk of this site.
    communicateRisk = function() {
      return(private$p_cwb$risk(private$p_data[[private$p_cwb$getTarget()]], private$p_pred))
    },

    #' @description
    #' Update the CWB model.
    #'
    #' @param bln (`character(1L)`)\cr
    #'   Name of the base learner to update
    #' @param param (`numeric()`)\cr
    #'   Parameter update.
    update = function(bln, param) {
      checkmate::assertChoice(bln, choices = names(private$p_bls))
      private$p_cwb$update(bln, param)
      private$p_pred = private$p_offset + private$p_cwb$linPred(private$p_bls)
      private$p_pr = private$p_cwb$pseudoResiduals(
        y = private$p_data[[private$p_cwb$getTarget()]], y_pred = private$p_pred)

      return(invisible(NULL))
    },

    #' @description
    #' Communicate (get) the static part of the base learner.
    #' Static parts are X^TX, the degrees of freedom, penalty,
    #' and the penalty matrix.
    #'
    #' @return
    #'   Returns list containing the static parts per base learner.
    communicateFittingParts = function() {
      if ((private$p_is_initialized_bl + private$p_is_initialized_offset) == 0)
        stop("Initialize prediction and base learner first.")

      out = lapply(private$pBlsInclude(), function(bl) {
        list(XtX = bl$getXtX(), pen = bl$getPenalty(),
          K = bl$getPenaltyMat(), df = bl$getDF())
      })
      return(out)
    },

    #' @description
    #' Communicate (get) the dynamic part X^Tu which changes after each
    #' update due to the adjustment of the pseudo residuals.
    #'
    #' @return
    #'   Returns list containing X^Tu of the current state of the model.
    communicateXtu = function() {
      if ((private$p_is_initialized_bl + private$p_is_initialized_offset) == 0)
        stop("Initialize prediction and base learner first.")

      out = lapply(private$pBlsInclude(), function(bl) list(Xtu = bl$getXty(private$p_pr)))
      return(out)
    },

    #' @description
    #' Communicate (get) the parts required to calculate a "global"
    #' initialization. For splines, this is the minimal and maximum of a feature
    #' wihle the classes are communicated for the ridge base learner.
    #'
    #' @return
    #'   Returns list the initialization parts for each base learner.
    communicateInit = function() {
      out = lapply(private$p_bls, function(bl) bl$communicateInit(private$p_data))
      names(out) = names(private$p_bls)
      return(out)
    },

    #' @description
    #' Calculate the constant initialization.
    #'
    #' @param aggregator (`function(y)`)\cr
    #'   Loss optimal initialization depending on y.
    #' @param ... \cr
    #'   Arguments passed to `aggregator`.
    #'
    #' @return
    #'   Aggregated value by calling the `aggregator` on the site data.
    getAggregatedInitialization = function(aggregator, ...) {
      checkmate::assertFunction(aggregator, args = "y")

      target = private$p_cwb$getTarget()
      checkmate::assertChoice(target, choices = colnames(private$p_data))

      out = aggregator(private$p_data[[private$p_cwb$getTarget()]], ...)
      checkmate::assertNumeric(out, len = 1L)
      return(out)
    },

    #' @description
    #' Initialize the prediction and pseudo residuals on the site.
    #'
    #' @param offset (`numeric(1L)`)\cr
    #'   Constant initialization for the prediction.
    initPrediction = function(offset) {
      if (private$p_is_initialized_offset) stop("Prediction is already initialized.")
      checkmate::assertNumeric(offset, len = 1L)
      private$p_offset = offset
      private$p_pred = rep(offset, nrow(private$p_data))
      private$p_is_initialized_offset = TRUE
      private$p_pr = private$p_cwb$pseudoResiduals(
        y = private$p_data[[private$p_cwb$getTarget()]], y_pred = private$p_pred)

      return(invisible(NULL))
    },

    #' @description
    #' Get the sum of squared errors of the current model state.
    #'
    #' @return
    #'   The sum of squared errors as numeric value.
    sse = function() {
      return(sum((private$p_data[[private$p_cwb$getTarget()]] - private$p_pred)^2))
    },

    #' @description
    #' Get the sum of squared errors of a base learner with parameter
    #' vector `param` and the pseudo residuals. This is required to
    #' get the best base learner in one iteration.
    #'
    #' @return
    #'   The sum of squared errors as numeric value.
    ssePrUpdate = function(bln, param) {
      checkmate::assertChoice(bln, choices = names(private$p_bls))
      return(sum((private$p_pr - private$p_bls[[bln]]$linPred(param))^2))
    },

    #' @description
    #' Get the number of columns of the site data.
    #'
    #' @return
    #'   Integer value containing the number of columns.
    ncol = function() {
      return(ncol(private$p_data))
    },

    #' @description
    #' Get the number of rows of the site data.
    #'
    #' @return
    #'   Integer value containing the number of rows
    nrow = function() {
      return(nrow(private$p_data))
    }

    ### JUST FOR TESTING, NEEDS TO BE REMOVED FOR A REALISTIC SETUP
    #getPR = function() {
      #return(private$p_pr)
    #},
    #getPred = function() {
      #return(private$p_pred)
    #}
  ),
  private = list(

    #' @field p_name (`character(1L`)\cr
    #'   Name/id of the site.
    p_name = NULL,

    #' @field p_data (`data.frame()`)\cr
    #'   Site data.
    p_data = NULL,

    #' @field p_bls (`list([Baselearner])`)\cr
    #'   List of [Baselearner]s.
    p_bls = NULL,

    #' @field p_cwb (`[CWB]`)\cr
    #'   Site CWB model (same as the hots CWB model).
    p_cwb = NULL,

    #' @field p_offset (`numeric(1L`)\cr
    #'   Offset of the site (or global offset over multiple sites).
    p_offset = NULL,

    #' @field p_pred (`numeric()`)\cr
    #'   Prediction of the current model state.
    p_pred = NULL,

    #' @field p_pr (`numeric()`)\cr
    #'   Pseudo residuals of the current model state.
    p_pr = NULL,

    #' @field p_is_initialized_offset (`logical(1L)`)\cr
    #'   Flag indicating whether the offset and hence prediction
    #'   is initialized or not.
    p_is_initialized_offset = FALSE,

    #' @field p_is_initialized_bl (`logical(1L)`)\cr
    #'   Flag indicating whether the base learner are initilized or not.
    p_is_initialized_bl = FALSE,

    #' @description
    #' Get the base learners that are used for modelling.
    #' I.e. the base learners with the `include = TRUE` flag.
    #'
    #' @return
    #'   List containing the base learners.
    pBlsInclude = function() {
      bls_out = list()
      k = 1
      for (i in seq_along(private$p_bls)) {
        if (private$p_bls[[i]]$isIncluded()) {
          bls_out[[k]] = private$p_bls[[i]]
          k = k + 1
        }
      }
      names(bls_out) = sapply(bls_out, function(bl) bl$getName())
      return(bls_out)
    }
  )
)





if (FALSE) {

nsim = 100L
x = runif(nsim, 0, 10)
dat = data.frame(x = x, y = sin(x))
knots = compboostSplines::createKnots(x, 10, 2)

loss  = function(y, y_pred) 0.5 * (y - y_pred)^2
dloss = function(y, y_pred) -(y - y_pred)

cwb = CWB$new("y", lr = 0.01, loss = loss, dloss = dloss)
bl = BaselearnerPSpline$new("x", ord = 2)
bl$initDesign(dat, knots)
param = solve(bl$getXtX() + bl$getPenalty() * bl$getPenaltyMat()) %*% bl$getXty(dat$y)

bls = list(BaselearnerPSpline$new("x", ord = 2))
site = Site$new("test-site", dat, cwb, bls)

site$communicateInit()

init = list(x_spline = list(knots = knots))
site$initBaselearner(init)
offset = site$getAggregatedInitialization(function(y) mean(y))
site$initPrediction(offset)
site$sse()
sum((dat$y - mean(dat$y))^2)

a = site$communicate()

site$ssePrUpdate("x_spline", param)
cwb$update("x_spline", param)

site$sse()
site$update("x_spline", param)
site$sse()

all.equal(cwb$bl_trace[[1]]$param, param)
all.equal(cwb$bl_map[[1]], param * 0.01)

all.equal(cwb$pseudoResiduals(dat$y, offset), dat$y - offset)


}
