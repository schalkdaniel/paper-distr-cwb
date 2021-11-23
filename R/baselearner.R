library(R6)
library(Matrix)

Baselearner = R6Class("Baselearner",
  public = list(
    initialize = function(feat, name,  df, include) {
      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertCharacter(name, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)
      checkmate::assertLogical(include, len = 1L)

      private$feat = feat
      private$name = name
      private$df = df
      private$include = include

      if (grepl("<<>>", feat)) private$is_combined = TRUE
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
    },
    isCombination = function() {
      return(private$is_combined)
    },
    isIncluded = function () {
      return(private$include)
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
    is_combined = FALSE,
    include = TRUE,
    setDesign = function(design) {
      private$design = design
      private$xtx = t(design) %*% design
    }
  ))

BaselearnerPSpline = R6Class("BaselearnerPSpline",
  inherit = Baselearner,
  public = list(
    initialize = function(feat, df = 4, n_knots = 20L, ord = 3, derivs = 2L, include = TRUE) {

      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)
      checkmate::assertIntegerish(n_knots, len = 1L, lower = 1)
      checkmate::assertIntegerish(ord, len = 1L, lower = 1)
      checkmate::assertIntegerish(derivs, len = 1L, lower = 1)

      super$initialize(feat, paste0(feat, "_spline"), df, include)

      private$n_knots = n_knots
      private$ord     = ord
      private$derivs  = derivs
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
    },
    createKnots = function(x) {
      return(compboostSplines::createKnots(x, private$n_knots, private$ord))
    }
  ),
  private = list(
    n_knots = NULL,
    knots   = NULL,
    ord     = NULL,
    derivs  = NULL
  )
)

BaselearnerRidge = R6Class("BaselearnerRidge",
  inherit = Baselearner,
  public = list(
    initialize = function(feat, df = 4, include = TRUE) {

      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)

      super$initialize(feat, paste0(feat, "_ridge"), df, include)
    },
    initDesign = function(dat, classes, set_penalty = FALSE) {

      checkmate::assertCharacter(classes)
      private$classes = classes

      if (! private$locked) {
        private$setDesign(self$design(dat))

        private$penalty_mat = diag(rep(1, length(classes)))
        if (set_penalty) {
          private$penalty = compboostSplines::demmlerReinsch(private$xtx,
            private$penalty_mat, private$df)
        }
        private$locked = TRUE
      }
    },
    design = function(dat) {
      mat = do.call(rbind, lapply(dat[[private$feat]], function(x) {
        rout = rep(0, length(private$classes))
        rout[which(private$classes == x)] = 1
        return(rout)
      }))
      return(mat)
    },
    communicateInit = function(dat) {
      x = dat[[private$feat]]
      return(list(classes = unique(x)))
    }
  ),
  private = list(
    classes = NULL
  )
)


BaselearnerCombine = R6Class("BaselearnerCombine",
  inherit = Baselearner,
  public = list(
    initialize = function(bl1, bl2, anisotrop = TRUE, include = TRUE) {

      checkmate::assertR6(bl1, "Baselearner")
      checkmate::assertR6(bl2, "Baselearner")

      feat = paste0(bl1$getFeat(), "<<>>", bl2$getFeat())
      df   = bl1$getDF() * bl2$getDF()
      super$initialize(feat, paste0(feat, "_combine"), df, include)

      private$bl1_name = bl1$getName()
      private$bl2_name = bl2$getName()
    },
    initDesign = function(dat, bl1, bl2, set_penalty = FALSE) {

      private$bl1 = bl1
      private$bl2 = bl2

      design1 = private$bl1$design(dat)
      design2 = private$bl2$design(dat)

      design = private$kronDesign(design1, design2)

      K1 = private$bl1$getPenaltyMat()
      K2 = private$bl2$getPenaltyMat()

      K = private$kronPenalty(K1, K2)

      if (! private$locked) {
        private$setDesign(design)
        private$penalty_mat = K
        if (set_penalty) {
          private$penalty = compboostSplines::demmlerReinsch(as.matrix(private$xtx),
            private$penalty_mat, private$df)
        }
        private$locked = TRUE
      }
    },
    design = function(dat) {
      design1 = private$bl1$design(dat)
      design2 = private$bl2$design(dat)

      return(private$kronDesign(design1, design2))
    },
    communicateInit = function(dat) {
      return(list())
    },
    getBl1Name = function() return(private$bl1_name),
    getBl2Name = function() return(private$bl2_name)
  ),
  private = list(
    bl1 = NULL,
    bl2 = NULL,
    bl1_name = NULL,
    bl2_name = NULL,
    kronDesign = function(design1, design2) {
      kr1 = kronecker(design1, Matrix::Matrix(1, ncol = ncol(design2)))
      kr2 = kronecker(Matrix::Matrix(1, ncol = ncol(design1)), design2)
      kr1 * kr2
    },
    kronPenalty = function(K1, K2) {
      kronecker(K1, diag(ncol(K2))) + kronecker(diag(ncol(K1)), K2)
    }
  )
)


if (FALSE) {
  dat = data.frame(x = rnorm(20), site = sample(LETTERS[1:5], 20, TRUE))

  bx = BaselearnerPSpline$new("x", ord = 2)
  knots = compboostSplines::createKnots(unlist(bx$communicateInit(dat)), 10, 2)
  bx$initDesign(dat, knots, TRUE)


  bl = BaselearnerRidge$new("site", include = FALSE)
  classes = bl$communicateInit(dat)
  bl$initDesign(dat, classes[[1]], TRUE)

  ## Check if penalty is calculated correctly:
  mboost:::df2lambda(X = bl$design(dat), df = 4, dmat = bl$getPenaltyMat(), XtX = bl$getXtX())
  bl$getPenalty()
  bl$getPenaltyMat()
  bl$getXtX()


  blc = BaselearnerCombine$new(bx, bl)
  blc$initDesign(dat, bx, bl, TRUE)

  blc$getPenalty()
  blc$getPenaltyMat()
}


