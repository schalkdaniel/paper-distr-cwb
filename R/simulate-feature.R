src_pen = "
arma::mat penaltySumKronecker (const arma::mat& Pa, const arma::mat& Pb)
{
  // Variables
  arma::mat out;

  // Create Diagonal matrices
  arma::mat eyePa = arma::diagmat( arma::vec(Pa.n_cols, arma::fill::ones) );
  arma::mat eyePb = arma::diagmat( arma::vec(Pb.n_cols, arma::fill::ones) );

  // sum of Kroneckers with diagonal marices
  out = arma::kron(Pa, eyePb) + arma::kron(eyePa, Pb);

  return out;
}"

src_center = "
arma::mat centerDesignMatrix (const arma::mat& X1, const arma::mat& X2)
{
  // Cross Product X1 and X2
  arma::mat cross = X1.t() * X2;

  // QR decomp
  // We require and orthogonal matrix Q
  arma::mat R;
  arma::mat Q;
  arma::qr(Q,R,cross);

  // get rank of R and add 1
  int rankR = arma::rank(R);

  // construct Z from rows 0 to last row and column R+1 to last column
  arma::mat Z = Q.cols(rankR,Q.n_cols-1);

  return Z;
}"

Rcpp::cppFunction(src_pen, "RcppArmadillo")
Rcpp::cppFunction(src_center, "RcppArmadillo")

#' Get linear predictor from B-spline.
#'
#' @param x [numeric] Vector of x values
#' @param bs_dim [integer(1)] Number of base functions for the spline
#'   (default = 10L). (Corresponds to number of inner knots for
#'   the spline).
#' @param sigma [numeric(1)] Standard deviation for the normally
#'   distributed random variable from which the parameter are
#'   drawn (default = 3).
#' @param offset [numeric(1)] Shift on the y-axis of the linear
#'   predictor (default = 0).
#' @param stz [logical(1)] Sum to zero constraint for the spline.
#' @return A list with \code{x}, \code{y}, and specific information
#'   about the spline such as parameter.
simSpline = function(x, bs_dim = 10L, sigma = 3, offset = 0, stz = FALSE) {
  checkmate::assertNumeric(x = sigma, len = 1L)
  checkmate::assertNumeric(x = offset, len = 1L)
  #if (bs_dim < 7) stop("Need bs_dim >= 7 !")

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
  if (stz)
    X = X %*% centerDesignMatrix(X, cbind(rep(1, nrow(X))))

  # multiply with random coefficients to get random functions
  coefs = rnorm(ncol(X), sd = sigma)

  return (list(y = X %*% coefs + offset, x = x, X = X, offset = offset, coefs = coefs, knots = kn))
}

#' Simulate feature with site specific effects
#'
#' @param n       [integer(1)] Number of observations.
#' @param n_sites [integer(1)] Number sites (randomly drawn).
#' @param from    [numeric(1)] Lower boundary of the feature.
#' @param up      [numeric(1)] Upper boundary of the feature.
#' @return A list with \code{x}, \code{y}, the site, and specific
#'   information about all splines (for main and site effects).
simulateFeature = function(n = 1000L, n_sites = 3L, from = NULL, to = NULL, offset = 0, ...) {
  checkmate::assertIntegerish(x = n, len = 1L)
  checkmate::assertIntegerish(x = n_sites, len = 1L)

  if (is.null(from)) from = runif(1, -100, 100)
  if (is.null(to)) to = from + runif(1, 0, 100)

  checkmate::assertNumeric(x = from, len = 1L, upper = to)
  checkmate::assertNumeric(x = to, len = 1L, lower = from)

  x = seq(from = from, to = to, length.out = n)
  main_effect = simSpline(x, offset = offset, ...)

  idx_sites = sample(seq_len(n_sites), n, TRUE)
  site_effects = lapply(seq_len(n_sites), function(i) simSpline(x[idx_sites == i],
    offset = 0, stz = TRUE, ...))

  y = main_effect$y
  for (i in seq_len(n_sites)) {
    site_effects[[i]]$y = site_effects[[i]]$y - mean(site_effects[[i]]$y)
    y_site = site_effects[[i]]$y
    y[idx_sites == i] = y[idx_sites == i] + y_site
  }

  out = list()
  out$splines = list(main = main_effect, site_effects = site_effects)
  out$data = data.frame(y = y, x = x, site = idx_sites)

  return(out)
}
