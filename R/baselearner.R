library(R6)

Baselearner = R6Class("Baselearner",
  public = list(
    initialize = function(feat, name,  df) {
      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertCharacter(name, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)

      private$feat = feat
      private$name = name
      private$df = df
    },
    getFeat = function() {
      return(private$feat)
    },
    getName = function() {
      return(private$name)
    },
    linPred = function(param) {
      return(private$design %*% param)
    },
    getXtX = function() {
      return(private$xtx)
    },
    getXty = function(y) {
      return(t(private$design) %*% y)
    },
    getFeature = function() {
      return(private$feat)
    },
    getDF = function() {
      return(private$df)
    },
    setPenalty = function(penalty) {
      checkmate::assertNumeric(penalty, len = 1L, lower = 0)
      if (! is.null(private$penalty))
        stop("Penalty was already set.")
      private$penalty = penalty
    },
    getPenalty = function() {
      return(private$penalty)
    },
    getPenaltyMat = function() {
      return(private$penalty_mat)
    },
    isLocked = function() {
      return(private$locked)
    }
  ),
  private = list(
    feat = NULL,
    name = NULL,
    design = NULL,
    xtx = NULL,
    df = NULL,
    penalty = NULL,
    penalty_mat = NULL,
    locked = FALSE,
    setDesign = function(design) {
      private$design = design
      private$xtx = t(design) %*% design
    }
  ))

BaselearnerPSpline = R6Class("BaselearnerPSpline",
  inherit = Baselearner,
  public = list(
    initialize = function(feat, df = 4, ord = 3, derivs = 2L) {

      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)
      checkmate::assertIntegerish(ord, len = 1L, lower = 1)
      checkmate::assertIntegerish(derivs, len = 1L, lower = 1)

      super$initialize(feat, paste0(feat, "_spline"), df)

      private$ord = ord
      private$derivs = derivs
    },
    initDesign = function(dat, knots, set_penalty = FALSE) {
      checkmate::assertNumeric(knots)
      private$knots = knots

      if (! private$locked) {
        private$setDesign(compboostSplines::createSplineBasis(dat[[private$feat]],
          private$ord, knots))

        private$penalty_mat = compboostSplines::penaltyMat(ncol(super$getXtX()), private$derivs)
        if (set_penalty) {
          private$penalty = compboostSplines::demmlerReinsch(private$xtx,
            private$penalty_mat, private$df)
        }

        private$locked = TRUE
      }
    },
    design = function(dat) {
      return(compboostSplines::createSplineBasis(dat[[private$feat]],
          private$ord, private$knots))
    },
    communicateInit = function(dat) {
      x = dat[[private$feat]]
      return(list(min = min(x), max = max(x)))
    }
  ),
  private = list(
    knots = NULL,
    ord = NULL,
    derivs = NULL
  )
)
