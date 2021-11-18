splineBL = function(feat, nknots, df = 4, ord = 3, derivs = 2, sparse = TRUE) {
  out = list(feat = feat,
    bln = paste0(feat, "_spline"),
    nknots = nknots,
    df = df)
  class(out) = "Spline"
  return(out)
}

ridgeBL = function(feat, df = 1) {
  out = list(feat = feat,
    bln = paste0(feat, "_ridge"),
    df = df)
  class(out) = "Ridge"
  return(out)
}

combinedBL = function(bl1, bl2, anistrop = FALSE) {
  out = list(feat = feat,
    bln = paste0(feat, "_combined"),
    df = df)
  class(out) = "Combined"
  return(out)
}


bls = list(splineB)

#' Distributed CWB prototype
#'
#' @param bls [`list()`]
#' @param iters [`integer(1L)`] Number of boosting iterations.
#' @param learning_rate [`numeric(1L)`] Learning rate for the fitting.
#' @return List containing the trace of fitted base learners.
distCWB = function (site_data, target, bls, iters, learning_rate, loss) {

  K = length(site_data)
  bl_classes = c("Spline", "Ridge", "Combined")
  bl_names = vapply(X = bls, FUN.VALUE = character(1L), FUN = function(bl) bl$bln)

  # Checks:
  checkmate::assertList(site_data)
  nuisance = lapply(site_data, checkmate::assertDataFrame)
  checkmate::assertChoice(target, choices = colnames(sites_data[[1]]))
  checkmate::assertList(bls)
  nuisance = lapply(bls, function(bl) {
    checkmate::assertChoice(bl$feat, choices = colnames(site_data[[1]]))
    checkmate::assertChoice(class(bl), bl_classes)
  })
  checkmate::assertIntegerish(iters, len = 1L, lower = 1L)
  checkmate::assertNumeric(learning_rate, len = 1L, lower = 0)

  # The objects have the following terminology:
  # - lls: Data at a sites
  # - llh: Data at the host
  lls = list()
  llh = list()

  # Rules:
  # - lls has a slot called 'share' which is allowed to get shared with the host

  # Initialize

  ## Generate site "packages" and initialize
  for (k in seq_len(K)) {
    init = lapply(bls, function(bl) do.call(paste0("initMeta", class(bl)), list(site_data[[k]])))
    lls[[k]] = list(data = site_data[[k]], target = target,
      share = list(bls = bls, init = init))
  }
  for (bln in bl_names) {
    binit = lapply(lls, )
  }
  lapply(lls, )



  ## Get init data per base learner
  iinit = lapply(bls, function(bl) do.call(paste0("init", class(bl), list(site_data))))

  ## Aggregate init data
  iaggr = lapply(iinit, function(init) )

  # Fitting
}
