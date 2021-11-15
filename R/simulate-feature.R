#' Get linear predictor from B-spline.
#'
#' @param x [numeric] Vector of x values
#' @param bs_dim [integer(1)] Number of base functions for the spline (default = 10L). (Corresponds to number of inner knots for the spline).
#' @param sigma [numeric(1)] Standard deviation for the normally distributed random variable from which the parameter are drawn (default = 3).
#' @param offset [numeric(1)] Shift on the y-axis of the linear predictor (default = 0).
#' @return The sum of \code{x} and \code{y}.
simSpline = function(x, bs_dim = 10L, sigma = 3, offset = 0, ...) {
  checkmate::assertNumeric(x = sigma, len = 1L)
  checkmate::assertNumeric(x = offset, len = 1L)
  if (bs_dim < 7) stop("Need bs_dim >= 7 !")

  nk = bs_dim - 2

  xu = max(x)
  xl = min(x)

  xr = xu - xl
  xl = xl - xr * 0.001
  xu = xu + xr * 0.001

  dx = (xu - xl)/(nk - 1)
  kn = seq(xl - dx * 3, xu + dx * 3, length = nk + 4 + 2)

  # create the spline basis functions
  X = splines::spline.des(kn, x, 4, x * 0)$design

  # multiply with random coefficients to get random functions
  coefs = rnorm(bs_dim, sd = sigma)

  return (list(y = X %*% coefs + offset, x = x, X = X, offset = offset, coefs = coefs, knots = kn))
}

simulateFeature = function(n = 1000L, n_sites = 3L, from = NULL, to = NULL, ...) {
  checkmate::assertIntegerish(x = n, len = 1L)
  checkmate::assertIntegerish(x = n_sites, len = 1L)

  if (is.null(from)) from = runif(1, -100, 100)
  if (is.null(to)) to = from + runif(1, 0, 100)

  checkmate::assertNumeric(x = from, len = 1L, upper = to)
  checkmate::assertNumeric(x = to, len = 1L, lower = from)


  x = seq(from = from, to = to, length.out = n)
  main_effect = simSpline(x, ...)

  idx_sites = sample(seq_len(n_sites), n, TRUE)
  site_effects = lapply(seq_len(n_sites), function(i) simSpline(x[idx_sites == i], ...))

  y = main_effect$y
  for (i in seq_len(n_sites)) {
    y[idx_sites == i] = site_effects[[i]]$y
  }

  out = list()
  out$splines = list(main = main_effect, site_effects = site_effects)
  out$data = data.frame(y = y, x = x, site = idx_sites)

  return(out)
}


devtools::install_github("schalkdaniel/compboost", ref = "dev")
library(compboost)

dat = simulateFeature()

# Define data:
dsx = InMemoryData$new(cbind(dat$data$x), "x")
dsc = CategoricalDataRaw$new(dat$data$site, "c")

facx = BaselearnerPSpline$new(dsx, "spline", list(df = 10, n_knots = 10))
facc = BaselearnerCategoricalRidge$new(dsc, "category", list(df = 1))

tensor1 = BaselearnerTensor$new(facx, facc, "tensor")

fl = BlearnerFactoryList$new()

fl$registerFactory(facx)
fl$registerFactory(tensor1)

loss = LossQuadratic$new()
optimizer = OptimizerCoordinateDescent$new()

log_iterations = LoggerIteration$new("iter", TRUE, 2000)

# Define new logger list:
logger_list = LoggerList$new()

# Register the logger:
logger_list$registerLogger(log_iterations)

cboost = Compboost_internal$new(
  response      = ResponseRegr$new("y", cbind(dat$data$y)),
  learning_rate = 0.05,
  stop_if_all_stopper_fulfilled = FALSE,
  factory_list = fl,
  loss         = loss,
  logger_list  = logger_list,
  optimizer    = optimizer
)

cboost$train(trace = TRUE)

# Dummy to plot baselearner traces
mod = list(model = cboost)
mod$getSelectedBaselearner = cboost$getSelectedBaselearner

plotBaselearnerTraces(mod)

df = dat$data
df$y_pred_tensor = cboost$predictFactoryTrainData("x_c_tensor")
df$y_pred_main   = cboost$predictFactoryTrainData("x_spline")
df$y_main        = dat$splines$main$y

gg_main = ggplot() +
  geom_line(data = df, aes(x = x, y = y_pred_main, color = "Predicted effect")) +
  geom_line(data = df, aes(x = x, y = y_main, color = "True effect")) +
  ggtitle("Main effect")

ggs = lapply(unique(dat$data$site), function(i) {
  df_tmp  = df[df$site == i, ]
  df_site = data.frame(x = dat$splines$site_effects[[i]]$x,
                       y = dat$splines$site_effects[[i]]$y)

  ggplot() +
    geom_line(data = df_tmp, mapping = aes(x = x, y = y_pred_tensor + y_pred_main, color = "Predicted effect")) +
    geom_line(data = df_site, mapping = aes(x = x, y = y, color = "True effect")) +
    ggtitle(paste0("Site ", i))
})

do.call(gridExtra::grid.arrange, c(list(gg_main), ggs))
