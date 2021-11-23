library(R6)
library(Matrix)

#' @title Abstract base learner class
#'
#' @description
#' This abstract class holds information requires for all base learners.
#' The information are the degrees of freedom, design matrix X, X^TX, as well as
#' the penalty and penalty matrix. The class also defines custom routines
#' such as calculating the linear predictor for a given parameter (`X %*% param`).
#'
#' @export
Baselearner = R6Class("Baselearner",
  public = list(

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param feat (`character(1L)`)\cr
    #'   Name of the feature this base learner applies to.
    #' @param name (`character(1L)`)\cr
    #'   Specific name (or id) of the base learner.
    #' @param df (`numeric(1L)`)\cr
    #'   Degrees of freedom used for this base learner.
    #' @param include (`character(1L)`)\cr
    #'   Flag to specify if the base learner should be included in
    #'   a later fitting process conducted by the [Host] class.
    initialize = function(feat, name,  df, include) {
      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertCharacter(name, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)
      checkmate::assertLogical(include, len = 1L)

      private$p_feat = feat
      private$p_name = name
      private$p_df = df
      private$p_include = include

      if (grepl("<<>>", feat)) private$p_is_combined = TRUE
    },

    #' @description
    #' Getter for the feature name.
    #'
    #' @return
    #'   Returns the feature name as character value.
    getFeat = function() {
      return(private$p_feat)
    },

    #' @description
    #' Getter for the base learner name.
    #'
    #' @return
    #'   Returns the base learner name as character value.
    getName = function() {
      return(private$p_name)
    },

    #' @description
    #' Calculate linear prediction for a given parameter vector.
    #'
    #' @param param (`numeric()`)\cr
    #'   Parameter vector multiplied with the base learner design matrix.
    #'
    #' @return
    #'   Returns the matrix product `$design %*% param`.
    linPred = function(param) {
      return(private$p_design %*% param)
    },

    #' @description
    #' Getter for X^TX.
    #'
    #' @return
    #'   Returns the X^TX matrix.
    getXtX = function() {
      return(private$p_xtx)
    },

    #' @description
    #' Getter for X^Ty with variable y.
    #'
    #' @param y (`numeric()`)\cr
    #'   Vector which is multiplied to X^T.
    #'
    #' @return
    #'   Returns the X^Ty matrix.
    getXty = function(y) {
      return(t(private$p_design) %*% y)
    },

    #' @description
    #' Getter for the degrees of freedom.
    #'
    #' @return
    #'   Returns the degrees of freedom as numeric value.
    getDF = function() {
      return(private$p_df)
    },

    #' @description
    #' Setter for the penalty.
    #'
    #' @param penalty (`numeric(1)`).
    setPenalty = function(penalty) {
      checkmate::assertNumeric(penalty, len = 1L, lower = 0)
      if (! is.null(private$p_penalty))
        stop("Penalty was already set.")
      private$p_penalty = penalty
      return(invisible(NULL))
    },

    #' @description
    #' Getter for the penalty.
    #'
    #' @return
    #'   Returns the penalty as numeric value.
    getPenalty = function() {
      return(private$p_penalty)
    },

    #' @description
    #' Getter for the penalty matrix.
    #'
    #' @return
    #'   Returns the penalty matrix.
    getPenaltyMat = function() {
      return(private$p_penalty_mat)
    },

    #' @description
    #' Check whether the base learner is locked (initialized).
    #'
    #' @return
    #'   Returns a logical value indicating whether the base
    #'   learner is locked or not.
    isLocked = function() {
      return(private$p_locked)
    },

    #' @description
    #' Check whether the base learner is a combination or not.
    #'
    #' @return
    #'   Returns a logical value indicating whether the base
    #'   learner is a combination or not.
    isCombination = function() {
      return(private$p_is_combined)
    },

    #' @description
    #' Check whether the base learner should be included in the host fitting process.
    #'
    #' @return
    #'   Returns a logical value indicating whether the base
    #'   learner should be included in the fitting or not.
    isIncluded = function() {
      return(private$p_include)
    }
  ),
  private = list(
    #' @field p_feat (`character(1L)`)\cr
    #'   Feature name.
    p_feat = NULL,

    #' @field p_name (`character(1L)`)\cr
    #'   Base learner name.
    p_name = NULL,

    #' @field p_design (`matrix()`)\cr
    #'   Model matrix.
    p_design = NULL,

    #' @field p_xtx (`matrix()`)\cr
    #'   Cross product of the design matrix.
    p_xtx = NULL,

    #' @field p_df (`numeric(1L)`)\cr
    #'   Degrees of freedom.
    p_df = NULL,

    #' @field p_penalty (`numeric(1L)`)\cr
    #'   Penalty.
    p_penalty = NULL,

    #' @field p_penalty_mat (`matrix()`)\cr
    #'   Penalty matrix.
    p_penalty_mat = NULL,

    #' @field p_locked (`logical(1L)`)\cr
    #'   Flag indicating whether the base learner is
    #'   locked (initialized) or not.
    p_locked = FALSE,

    #' @field p_is_combined (`logical(1L)`)\cr
    #'   Flag indicating whether the base learner is
    #'   a combination of other base learners or not..
    p_is_combined = FALSE,

    #' @field p_include (`logical(1L)`)\cr
    #'   Flag indicating whether the base learner should be
    #'   included into the fitting process or not.
    p_include = TRUE,

    #' @description
    #' Set the design matrix and the cross product of the
    #' design matrix.
    #'
    #' @param design (`matrix()`)\cr
    #'   Model/design matrix.
    pSetDesign = function(design) {
      private$p_design = design
      private$p_xtx = t(design) %*% design
      return(invisible(NULL))
    }
  )
)

#' @title P-spline base learner class
#'
#' @description
#' This class defines a P-spline base learner.
#'
#' @export
BaselearnerPSpline = R6Class("BaselearnerPSpline",
  inherit = Baselearner,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param feat (`character(1L)`)\cr
    #'   Name of the feature this base learner applies to.
    #' @param df (`numeric(1L)`)\cr
    #'   Degrees of freedom used for this base learner.
    #' @param n_knots (`integer(1L)`)\cr
    #'   Number of knots used to build the model matrix.
    #' @param ord (`integer(1L)`)\cr
    #'   Polynomial degree/order of the base functions.
    #' @param derivs (`integer(1L)`)\cr
    #'   Number of differences that are penalized.
    #' @param include (`character(1L)`)\cr
    #'   Flag to specify if the base learner should be included in
    #'   a later fitting process conducted by the [Host] class.
    initialize = function(feat, df = 4, n_knots = 20L, ord = 3, derivs = 2L, include = TRUE) {

      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)
      checkmate::assertIntegerish(n_knots, len = 1L, lower = 1)
      checkmate::assertIntegerish(ord, len = 1L, lower = 1)
      checkmate::assertIntegerish(derivs, len = 1L, lower = 1)

      super$initialize(feat, paste0(feat, "_spline"), df, include)

      private$p_n_knots = n_knots
      private$p_ord     = ord
      private$p_derivs  = derivs
    },

    #' @description
    #' Initializes the P-spline design matrix.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   Data holding the feature for the base learner.
    #' @param knots (`numeric()`)\cr
    #'   Vector of knots (created by `$createKnots()`).
    #' @param set_penalty (`logical(1L)`)\cr
    #'   Should the penalty set just by calculating it from
    #'   the data of this base learner? The penalty can be
    #'   set later. This is e.g. done by the [Host] which
    #'   calculates a "global" penalty based on the global
    #'   X^TX matrix.
    initDesign = function(dat, knots, set_penalty = FALSE) {
      checkmate::assertNumeric(knots)
      private$p_knots = knots

      if (! private$p_locked) {
        private$pSetDesign(compboostSplines::createSplineBasis(dat[[private$p_feat]],
          private$p_ord, knots))

        private$p_penalty_mat = compboostSplines::penaltyMat(ncol(super$getXtX()), private$p_derivs)
        if (set_penalty) {
          private$p_penalty = compboostSplines::demmlerReinsch(private$p_xtx,
            private$p_penalty_mat, private$p_df)
        }

        private$p_locked = TRUE
      }
    },

    #' @description
    #' Calculates the P-spline design matrix specifically for
    #' an initialized base learner.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   New data holding the feature for the base learner.
    design = function(dat) {
      return(compboostSplines::createSplineBasis(dat[[private$p_feat]],
          private$p_ord, private$p_knots))
    },

    #' @description
    #' Communication information required to calculate static
    #' information (such as knots or classes) to build the
    #' design matrix for each base learner equally. For this
    #' base learner, the minimum and maximum of the feature is
    #' communicated.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   New data holding the feature for the base learner.
    communicateInit = function(dat) {
      x = dat[[private$p_feat]]
      return(list(min = min(x), max = max(x)))
    },

    #' @description
    #' Calculate knot values for a given numeri vector x.
    #'
    #' @param x (`numeric()`)\cr
    #'   Numeric vector for which the knots are created.
    createKnots = function(x) {
      return(compboostSplines::createKnots(x, private$p_n_knots, private$p_ord))
    }
  ),
  private = list(

    #' @field p_n_knots (`integer(1L)`)\cr
    #'   Number of knots.
    p_n_knots = NULL,

    #' @field p_knots (`numeric()`)\cr
    #'   Knots as numeric vector.
    p_knots   = NULL,

    #' @field p_ord (`integer(1L)`)\cr
    #'   Order of the polynomial base functions.
    p_ord     = NULL,

    #' @field p_derivs (`integer(1L)`)\cr
    #'   Number of differences that are penalized.
    p_derivs  = NULL
  )
)

#' @title Ridge base learner class
#'
#' @description
#' This class defines a ridge base learner.
#'
#' @export
BaselearnerRidge = R6Class("BaselearnerRidge",
  inherit = Baselearner,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param feat (`character(1L)`)\cr
    #'   Name of the feature this base learner applies to.
    #' @param df (`numeric(1L)`)\cr
    #'   Degrees of freedom used for this base learner.
    #' @param include (`character(1L)`)\cr
    #'   Flag to specify if the base learner should be included in
    #'   a later fitting process conducted by the [Host] class.
    initialize = function(feat, df = 4, include = TRUE) {

      checkmate::assertCharacter(feat, len = 1L)
      checkmate::assertNumeric(df, len = 1L, lower = 1)

      super$initialize(feat, paste0(feat, "_ridge"), df, include)
    },

    #' @description
    #' Initializes the ridge design matrix.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   Data holding the feature for the base learner.
    #' @param classes (`character()`)\cr
    #'   Vector all classes.
    #' @param set_penalty (`logical(1L)`)\cr
    #'   Should the penalty set just by calculating it from
    #'   the data of this base learner? The penalty can be
    #'   set later. This is e.g. done by the [Host] which
    #'   calculates a "global" penalty based on the global
    #'   X^TX matrix.
    initDesign = function(dat, classes, set_penalty = FALSE) {

      checkmate::assertCharacter(classes)
      private$p_classes = classes

      if (! private$p_locked) {
        private$pSetDesign(self$design(dat))

        private$p_penalty_mat = diag(rep(1, length(classes)))
        if (set_penalty) {
          private$p_penalty = compboostSplines::demmlerReinsch(private$p_xtx,
            private$p_penalty_mat, private$p_df)
        }
        private$p_locked = TRUE
      }
    },

    #' @description
    #' Calculates the ridge design matrix specifically for
    #' an initialized base learner.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   New data holding the feature for the base learner.
    design = function(dat) {
      mat = do.call(rbind, lapply(dat[[private$p_feat]], function(x) {
        rout = rep(0, length(private$p_classes))
        rout[which(private$p_classes == x)] = 1
        return(rout)
      }))
      return(mat)
    },

    #' @description
    #' Communication information required to calculate static
    #' information (such as knots or classes) to build the
    #' design matrix for each base learner equally. For this
    #' base learner, the classes of the feature are communicated.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   New data holding the feature for the base learner.
    communicateInit = function(dat) {
      x = dat[[private$p_feat]]
      return(list(classes = unique(x)))
    }
  ),
  private = list(
    #' @field p_classes (`character()`)\cr
    #'   Classes of the categorical feature.
    p_classes = NULL
  )
)

#' @title Conbine two base learners
#'
#' @description
#' This class combines two base learners based on
#' a row-wise Kronecker operation on the design matrices.
#'
#' @export
BaselearnerCombine = R6Class("BaselearnerCombine",
  inherit = Baselearner,
  public = list(
    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param bl1 (`[Baselearner]`)\cr
    #'   First base learner.
    #' @param bl2 (`[Baselearner]`)\cr
    #'   Second base learner.
    #' @param anisotrop (`logical(1L)`)\cr
    #'   Method used to calculate the penalty of the combined base learner.
    #'   If `anisotrop = TRUE` (default), the penalty is calculated equally for both
    #'   base learners. If `anisotrop = FALSE`, the direction of the penalty
    #'   for both base learners is taken into account.
    #' @param include (`character(1L)`)\cr
    #'   Flag to specify if the base learner should be included in
    #'   a later fitting process conducted by the [Host] class.
    initialize = function(bl1, bl2, anisotrop = TRUE, include = TRUE) {

      checkmate::assertR6(bl1, "Baselearner")
      checkmate::assertR6(bl2, "Baselearner")

      feat = paste0(bl1$getFeat(), "<<>>", bl2$getFeat())
      df   = bl1$getDF() * bl2$getDF()
      super$initialize(feat, paste0(feat, "_combine"), df, include)

      private$p_bl1_name = bl1$getName()
      private$p_bl2_name = bl2$getName()
    },

    #' @description
    #' Initializes the ridge design matrix.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   Data holding the feature for the base learner.
    #' @param bl1 (`[Baselearner]`)\cr
    #'   First base learner.
    #' @param bl2 (`[Baselearner]`)\cr
    #'   Second base learner.
    #' @param set_penalty (`logical(1L)`)\cr
    #'   Should the penalty set just by calculating it from
    #'   the data of this base learner? The penalty can be
    #'   set later. This is e.g. done by the [Host] which
    #'   calculates a "global" penalty based on the global
    #'   X^TX matrix.
    initDesign = function(dat, bl1, bl2, set_penalty = FALSE) {

      private$p_bl1 = bl1
      private$p_bl2 = bl2

      design1 = private$p_bl1$design(dat)
      design2 = private$p_bl2$design(dat)

      design = private$pKronDesign(design1, design2)

      K1 = private$p_bl1$getPenaltyMat()
      K2 = private$p_bl2$getPenaltyMat()

      K = private$pKronPenalty(K1, K2)

      if (! private$p_locked) {
        private$pSetDesign(design)
        private$p_penalty_mat = K
        if (set_penalty) {
          private$p_penalty = compboostSplines::demmlerReinsch(as.matrix(private$p_xtx),
            private$p_penalty_mat, private$p_df)
        }
        private$p_locked = TRUE
      }
    },

    #' @description
    #' Calculates the ridge design matrix specifically for
    #' an initialized base learner.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   New data holding the feature for the base learner.
    design = function(dat) {
      design1 = private$p_bl1$design(dat)
      design2 = private$p_bl2$design(dat)

      return(private$p_kronDesign(design1, design2))
    },

    #' @description
    #' Communication information required to calculate static
    #' information (such as knots or classes) to build the
    #' design matrix for each base learner equally. Nothing
    #' is communicated for base learner since it depends
    #' on the initialization of the two base learner it
    #' combines.
    #'
    #' @param dat (`data.frame()`)\cr
    #'   New data holding the feature for the base learner.
    communicateInit = function(dat) {
      return(list())
    },
    getBl1Name = function() return(private$p_bl1_name),
    getBl2Name = function() return(private$p_bl2_name)
  ),
  private = list(
    #' @field p_bl1 (`[Baselearner]`)\cr
    #'   First base learner.
    p_bl1 = NULL,

    #' @field p_bl2 (`[Baselearner]`)\cr
    #'   Second base learner.
    p_bl2 = NULL,

    #' @field p_bl1_name (`character(1L)`)\cr
    #'   Name of the first base learner.
    p_bl1_name = NULL,

    #' @field p_bl2_name (`character(1L)`)\cr
    #'   Name of the second base learner.
    p_bl2_name = NULL,

    #' @description
    #' Function to calculate the design matirx products
    #' via row-wise Kronecker.
    #'
    #' @param design1 (`matrix()`)\cr
    #'   First design matrix.
    #' @param design2 (`matrix()`)\cr
    #'   Second design matrix.
    #'
    #' @return
    #'   The matrix product.
    pKronDesign = function(design1, design2) {
      kr1 = kronecker(design1, Matrix::Matrix(1, ncol = ncol(design2)))
      kr2 = kronecker(Matrix::Matrix(1, ncol = ncol(design1)), design2)
      return(kr1 * kr2)
    },

    #' @description
    #' Function to calculate the penalty matrix
    #' via row-wise Kronecker.
    #'
    #' @param K1 (`matrix()`)\cr
    #'   First penalty matrix.
    #' @param K2 (`matrix()`)\cr
    #'   Second penalty matrix.
    #'
    #' @return
    #'   The matrix product.
    pKronPenalty = function(K1, K2) {
      return(kronecker(K1, diag(ncol(K2))) + kronecker(diag(ncol(K1)), K2))
    }
  )
)
