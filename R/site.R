library(R6)

#source(here::here("R/baselearner.R"))
#source(here::here("R/cwb.R"))

Site = R6Class("Site",
  public = list(
    #' @param name [`character(1L)`] Name of the site.
    #' @param data [`data.frame()`] Data frame exclusively hold by the server.
    #' @param cwb [`CWB`] Object of R6 class CWB.
    #' @param bls [`list()`] List of object of R6 class base learner.
    initialize = function(name, data, cwb, bls) {
      checkmate::assertCharacter(name, len = 1L)
      checkmate::assertDataFrame(data)
      checkmate::assertR6(cwb, "CWB")
      nuisance = lapply(bls, function(bl) checkmate::assertR6(bl, "Baselearner"))

      private$name = name
      private$data = data
      private$cwb  = cwb$clone()

      private$bls = lapply(bls, function(bl) bl$clone(deep = TRUE))
      names(private$bls) = vapply(private$bls, FUN.VALUE = character(1L),
        FUN = function(bl) bl$getName())
    },
    baselearnerNames = function() {
      return(names(private$bls))
    },
    initBaselearner = function(ll_init) {
      if (private$is_initialized_bl)
        stop("Baselearner already initialized. Call `$unlockBaselearner()` and re initialize.")

      ninits = names(ll_init)
      for (ni in ninits) {
        do.call(private$bls[[ni]]$initDesign, c(list(dat = private$data), ll_init[[ni]]))
      }
      if (any(! sapply(private$bls, function(bl) bl$isLocked())))
        warning("Not all base learners were initialized.")
      private$is_initialized_bl = TRUE
      return(invisible(NULL))
    },
    unlockBaselearner = function() {
      private$is_initialized_bl = FALSE
      return(invisible(NULL))
    },
    communicateRisk = function() {
      return(private$cwb$risk(private$data[[private$cwb$getTarget()]], private$pred))
    },
    #' @param bln [`character(1L)`] Name of the base learner to update
    #' @param param [`numeric()`] Parameter update.
    update = function(bln, param) {
      checkmate::assertChoice(bln, choices = names(private$bls))
      private$cwb$update(bln, param)
      private$pred = private$offset + private$cwb$linPred(private$bls)
      private$pr = private$cwb$pseudoResiduals(
        y = private$data[[private$cwb$getTarget()]], y_pred = private$pred)

      return(invisible(NULL))
    },
    communicate = function() {
      if ((private$is_initialized_bl + private$is_initialized_offset) == 0)
        stop("Initialize prediction and base learner first.")

      out = lapply(private$bls, function(bl) {
        list(XtX = bl$getXtX(), Xtu = bl$getXty(private$pr), pen = bl$getPenalty(),
          K = bl$getPenaltyMat())
      })
      return(out)
    },
    communicateInit = function() {
      out = lapply(bls, function(bl) bl$communicateInit(private$data))
      names(out) = names(private$bls)
      return(out)
    },
    #' @param aggregator [`function(y)`] Loss optimal initialization depending on y.
    #' @param ... Arguments passed to `aggregator`.
    getAggregatedInitialization = function(aggregator, ...) {
      checkmate::assertFunction(aggregator, args = "y")
      out = aggregator(private$data[[private$cwb$getTarget()]], ...)
      checkmate::assertNumeric(out, len = 1L)
      return(out)
    },
    #' @param offset [`numeric(1L)`] Constant initialization for the prediction
    initPrediction = function(offset) {
      if (private$is_initialized_offset) stop("Prediction is already initialized.")
      checkmate::assertNumeric(offset, len = 1L)
      private$offset = offset
      private$pred = rep(offset, nrow(private$data))
      private$is_initialized_offset = TRUE
      private$pr = private$cwb$pseudoResiduals(
        y = private$data[[private$cwb$getTarget()]], y_pred = private$pred)

      return(invisible(NULL))
    },
    sse = function() {
      return(sum((private$data[[private$cwb$getTarget()]] - private$pred)^2))
    },
    ssePrUpdate = function(bln, param) {
      checkmate::assertChoice(bln, choices = names(private$bls))
      return(sum((private$pr - private$bls[[bln]]$linPred(param))^2))
    },
    ncol = function() {
      return(ncol(private$data))
    },
    nrow = function() {
      return(nrow(private$data))
    },
    ### ATM JUST FOR TESTING, NEEDS TO BE REMOVED FOR A REALISTIC SETUP
    getPR = function() {
      return(private$pr)
    },
    getPred = function() {
      return(private$pred)
    }
  ),
  private = list(
    name = NULL,
    data = NULL,
    bls = NULL,
    cwb = NULL,
    offset = NULL,
    pred = NULL,
    pr = NULL,
    is_initialized_offset = FALSE,
    is_initialized_bl = FALSE
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
