library(R6)

CWB = R6Class("CWB",
  public = list(
    #' @param target [`character(1L)`] Target value. Must be contained by data.
    #' @param lr [`numeric(1L)`] Learning rate.
    #' @param loss [`function(y, y_pred)`] Loss with formal arguments `y` and `y_pred`.
    #' @param dloss [`function(y, y_pred)`] Derivative of the loss w.r.t. y_pred.
    initialize = function(target, lr = 0.01, loss, dloss) {
      checkmate::assertCharacter(target, len = 1L)
      checkmate::assertNumeric(lr, len = 1L, lower = 0)
      checkmate::assertFunction(loss, args = c("y", "y_pred"))
      checkmate::assertFunction(dloss, args = c("y", "y_pred"))

      private$target = target
      private$lr = lr
      private$loss = loss
      private$dloss = dloss
    },
    getLR = function() {
      return(private$lr)
    },
    update = function(bln, param) {
      private$bl_trace[[length(private$bl_trace) + 1L]] = list(name = bln, param = param)
      if (is.null(private$bl_map[[bln]])) private$bl_map[[bln]] = 0
      private$bl_map[[bln]] = private$bl_map[[bln]] + private$lr * param
    },
    linPred = function(bls) {
      blLinPred = function(bl) {
        if (! is.null(private$bl_map[[bl$getName()]]))
          return(bl$linPred(private$bl_map[[bl$getName()]]))
        else
          return(0)
      }
      pred = rowSums(do.call(cbind, lapply(bls, blLinPred)))
      return(pred)
    },
    pseudoResiduals = function(y, y_pred) {
      return(-private$dloss(y, y_pred))
    },
    getTarget = function() {
      return(private$target)
    },
    risk = function(y, y_pred) {
      return(mean(private$loss(y, y_pred)))
    },
    getBlMap = function() {
      return(private$bl_map)
    },
    getBlTrace = function(just_names = FALSE) {
      if (just_names)
        return(vapply(private$bl_trace, function(blt) blt$name, character(1L)))

      return(private$bl_trace)
    },
    blTable = function() {
      return(table(self$getBlTrace(just_names = TRUE)))
    }
  ),
  private = list(
    target = NULL,
    lr = NULL,
    loss = NULL,
    dloss = NULL,
    bl_map = list(),
    bl_trace = list()
  )
)

