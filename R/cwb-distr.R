library(R6)

source(here::here("R/baselearner.R"))
source(here::here("R/cwb.R"))
source(here::here("R/site.R"))
source(here::here("R/host.R"))

nsim = 1000L
x = runif(nsim, 0, 10)
svec = sample(LETTERS[1:4], nsim, TRUE)
dat = data.frame(x = x, y = sin(x))

# CWB structure:
bl_spline  = BaselearnerPSpline$new("x", ord = 2)
#bl_ridge   =
#bl_combine =

#bls = list(bl_spline, bl_combine)
bls = list(bl_spline)

# Local CWB:
loss  = function(y, y_pred) 0.5 * (y - y_pred)^2
dloss = function(y, y_pred) -(y - y_pred)

cwb = CWB$new("y", lr = 0.1, loss = loss, dloss = dloss)

# Define sites:
sites = lapply(unique(svec), function(s) {
  Site$new(s, dat[svec == s, ], cwb, bls) })


host = Host$new(cwb, bls, sites)

host$initializePrediction()
host$initializeBaselearner()

mstop = 1000L
for (m in seq_len(mstop)) {
  host$updateCWB(m)
}
host$cwb$getBlTrace()
host$cwb$getBlMap()


### Testing:
blst = bls
d0 = dat

lapply(blst, function(bl) bl$initDesign(d0, host$init_aggr$x_spline$knots))
d0$pred = host$offset + blst[[1]]$linPred(host$cwb$getBlMap()[["x_spline"]])

library(ggplot2)

ggplot(d0) +
  geom_point(aes(x = x, y = y, color = "Truth")) +
  geom_line(aes(x = x, y = pred, color = "Predicted"))


## Compare to compboost:
devtools::load_all("~/repos/compboost")
cboost = boostSplines(data = dat, target = "y", loss = LossQuadratic$new(),
  learning_rate = 0.1, iterations = 1000L, n_knots = 10, df = 4, degree = 2,
  differences = 2)

all.equal(cboost$baselearner_list$x_spline$factory$getPenaltyMat(),
  host$bl_parts$x_spline$K)

all.equal(cboost$baselearner_list$x_spline$factory$getPenalty(),
  host$bl_parts$x_spline$pen)

all.equal(check.attributes = FALSE,
  cboost$baselearner_list$x_spline$factory$getData() %*% t(cboost$baselearner_list$x_spline$factory$getData()),
   host$bl_parts$x_spline$XtX)

cf_cboost = cboost$getEstimatedCoef()
cf_dist   = host$cwb$getBlMap()

cf_cboost$offset
host$offset

all.equal(cf_cboost$x_spline, cf_dist$x_spline, check.attributes = FALSE)
all.equal(cboost$model$getRiskVector()[-1], host$risk)

