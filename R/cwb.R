library(R6)

#' @title CWB class
#'
#' @description
#' This holds all information about a CWB model. It also serves
#' as a container for the included base learner as well as the
#' estimated parameter vector.
#'
#' @export
CWB = R6Class("CWB",
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param target (`character(1L)`)\cr
    #'   Target value. Must be contained by data.
    #' @param lr (`numeric(1L)`)\cr
    #'   Learning rate.
    #' @param loss (`function(y, y_pred)`)\cr
    #'   Loss with formal arguments `y` and `y_pred`.
    #' @param dloss (`function(y, y_pred)`)\cr
    #'   Derivative of the loss w.r.t. y_pred.
    initialize = function(target, lr = 0.01, loss, dloss) {
      checkmate::assertCharacter(target, len = 1L)
      checkmate::assertNumeric(lr, len = 1L, lower = 0)
      checkmate::assertFunction(loss, args = c("y", "y_pred"))
      checkmate::assertFunction(dloss, args = c("y", "y_pred"))

      private$p_target = target
      private$p_lr = lr
      private$p_loss = loss
      private$p_dloss = dloss
    },

    #' @description
    #' Getter for the learning rate.
    #'
    #' @return
    #'   Returns the learning rate as numerical value.
    getLR = function() {
      return(private$p_lr)
    },

    #' @description
    #' Updates the CWB model by including a new base learner and
    #' updating the parameter map.
    #'
    #' @param bln (`character(1L)`)\cr
    #'   Name of the base learner which is added to the ensemble.
    #' @param param (`numeric()`)\cr
    #'   Parameter vector of this base learner.
    update = function(bln, param) {
      private$p_bl_trace[[length(private$p_bl_trace) + 1L]] = list(name = bln, param = param)
      if (is.null(private$p_bl_map[[bln]])) private$p_bl_map[[bln]] = 0
      private$p_bl_map[[bln]] = private$p_bl_map[[bln]] + private$p_lr * param
      return(invisible(NULL))
    },

    #' @description
    #' Calculate linear prediction for a given list of base learners.
    #'
    #' @param bls (`list([Baselearner])`)\cr
    #'   List of base learners for which the linear predictor is calculated.
    #'
    #' @return
    #'   Returns linear predictor.
    linPred = function(bls) {
      blLinPred = function(bl) {
        if (! is.null(private$p_bl_map[[bl$getName()]]))
          return(bl$linPred(private$p_bl_map[[bl$getName()]]))
        else
          return(0)
      }
      pred = rowSums(do.call(cbind, lapply(bls, blLinPred)))
      return(pred)
    },

    #' @description
    #' Update the pseudo residuals.
    #'
    #' @param y (`numeric()`)\cr
    #'   Response.
    #' @param y_pred (`numeric()`)\cr
    #'   Predicted values.
    #'
    #' @return
    #'   Pseudo residuals based on the loss.
    pseudoResiduals = function(y, y_pred) {
      return(-private$p_dloss(y, y_pred))
    },

    #' @description
    #' Getter for the target name.
    #'
    #' @return
    #'   Returns the target name as character value.
    getTarget = function() {
      return(private$p_target)
    },

    #' @description
    #' Get the empirical risk based on the response and prediction.
    #'
    #' @param y (`numeric()`)\cr
    #'   Response.
    #' @param y_pred (`numeric()`)\cr
    #'   Predicted values.
    #'
    #' @return
    #'   Risk value based on the loss.
    risk = function(y, y_pred) {
      return(mean(private$p_loss(y, y_pred)))
    },

    #' @description
    #' Get the base learner map containing the estimates
    #' of the paramters.
    #'
    #' @return
    #'   List with name and paramter vector.
    getBlMap = function() {
      return(private$p_bl_map)
    },

    #' @description
    #' Get the base learner trace.
    #'
    #' @param just_names (`logical(1L)`)\cr
    #'   Should the function just return the names without
    #'   parameters?.
    #'
    #' @return
    #'   List containing the names and (depnding on `just_names`) the
    #'   added parameter for each iteration.
    getBlTrace = function(just_names = FALSE) {
      if (just_names)
        return(vapply(private$p_bl_trace, function(blt) blt$name, character(1L)))

      return(private$p_bl_trace)
    },

    #' @description
    #' Get a table of all included base learners.
    #'
    #' @return
    #'   The result of table on the base learner trace with
    #'   `just_names = TRUE`.
    blTable = function() {
      return(table(self$getBlTrace(just_names = TRUE)))
    }
  ),
  private = list(

    #' @field p_target (`character(1L)`)\cr
    #'   Name of the target variable.
    p_target = NULL,

    #' @field p_lr (`numeric(1L)`)\cr
    #'   Learning rate.
    p_lr = NULL,

    #' @field p_loss (`function(y, y_pred)`)\cr
    #'   Loss function.
    p_loss = NULL,

    #' @field p_dloss (`function(y, y_pred)`)\cr
    #'   Derivative of the loss function.
    p_dloss = NULL,

    #' @field p_bl_map (`list()`)\cr
    #'   List containing the estimated parameters.
    p_bl_map = list(),

    #' @field p_bl_trace (`list()`)\cr
    #'   List containing all names and parameters of the
    #'   included base learners.
    p_bl_trace = list()
  )
)

