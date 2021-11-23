
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
library(R6)

source(here::here("R/baselearner.R"))
source(here::here("R/cwb.R"))
source(here::here("R/site.R"))
source(here::here("R/host.R"))
source(here::here("R/simulate-feature.R"))

# SPECS
# =======================================================

nsim    = 2000L
n_sites = 4L

mstop = 2000L
learning_rate = 0.1
df = 5
n_knots = 20L

# SIMULATE DATA
# =======================================================

set.seed(31415)
dat = simulateFeature(n = nsim, n_sites = n_sites, bs_dim = 10L, offset = 10)
SNR = 5

dat$data$ynn = rnorm(nsim, dat$data$y, sd = sd(dat$data$y) / SNR)
dat$data$site = as.character(dat$data$site)


# CWB STRUCTURE
# =======================================================

bl_spline  = BaselearnerPSpline$new("x", df = df)
bl_ridge   = BaselearnerRidge$new("site", df = 1, include = FALSE)
bl_combine = BaselearnerCombine$new(bl_spline, bl_ridge)

# Define base learner list
bls = list(bl_spline, bl_ridge,  bl_combine)

# Local CWB:
loss  = function(y, y_pred) 0.5 * (y - y_pred)^2
dloss = function(y, y_pred) -(y - y_pred)

cwb = CWB$new("ynn", lr = learning_rate, loss = loss, dloss = dloss)

# HOST/SITES
# =======================================================

# Define sites:
sites = lapply(unique(dat$data$site), function(s) {
  Site$new(s, dat$data[dat$data$site == s, c("ynn", "site", "x")], cwb, bls) })

# Define the host:
host = Host$new(cwb, bls, sites)


# DISTR. CWB WORKFLOW
# =======================================================

## Initialize:
## ----------------------------------

host$initializePrediction()
host$initializeBaselearner()

host$offset
#> [1] 10.01

## Fitting:
## ----------------------------------

for (m in seq_len(mstop))
  host$updateCWB(250L)
#> > 0:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x_spline 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 1.769 
#> > 250:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.1426 
#> > 500:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.08479 
#> > 750:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.07653 
#> > 1000:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.07435 
#> > 1250:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.07351 
#> > 1500:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.07309 
#> > 1750:
#>   >> Receive X^Tu from all sites for all base learners
#>   >> Aggregate to get global parameter estimates for all base learners
#>   >> Send best parameters to sites to get SSEs
#>   >> Found best base learner: x<<>>site_combine 
#>   >> Send parameter of best base learner to site to update CWB
#>   >> Risk: 0.07283

host$cwb$blTable()
#> 
#>          x_spline x<<>>site_combine 
#>                28              1972


# COMPBOOST COMPARISON
# =======================================================

devtools::load_all("~/repos/compboost")

dsx = InMemoryData$new(cbind(dat$data$x), "x")
dsc = CategoricalDataRaw$new(dat$data$site, "site")

blx = compboost::BaselearnerPSpline$new(dsx, "spline", list(df = df, n_knots = n_knots))
blc = BaselearnerCategoricalRidge$new(dsc, "category", list(df = 1))

tensor1 = BaselearnerTensor$new(blx, blc, "tensor", TRUE)

fl = BlearnerFactoryList$new()

fl$registerFactory(blx)
fl$registerFactory(tensor1)

loss = LossQuadratic$new()
optimizer = OptimizerCoordinateDescent$new()

newdata_full  = list(dsx, dsc)
response_full = ResponseRegr$new("y", cbind(dat$data$y))

log_iterations = LoggerIteration$new("iter", TRUE, mstop)
log_oob        = LoggerOobRisk$new("oob", TRUE, loss, 0.00, 10, newdata_full, response_full)

# Define new logger list:
logger_list = LoggerList$new()

# Register the logger:
logger_list$registerLogger(log_iterations)
logger_list$registerLogger(log_oob)

cboost = Compboost_internal$new(
  response      = ResponseRegr$new("y", cbind(dat$data$ynn)),
  learning_rate = learning_rate,
  stop_if_all_stopper_fulfilled = FALSE,
  factory_list = fl,
  loss         = loss,
  logger_list  = logger_list,
  optimizer    = optimizer
)

cboost$train(trace = 500)
#>    1/2000   risk = 1.8  oob = 1.7   
#>  500/2000   risk = 0.085  oob = 0.014   
#> 1000/2000   risk = 0.074  oob = 0.0037   
#> 1500/2000   risk = 0.073  oob = 0.0026   
#> 2000/2000   risk = 0.073  oob = 0.0023   
#> 
#> 
#> Train 2000 iterations in 0 Seconds.
#> Final risk based on the train set: 0.073



## Compare estimates:
## --------------------------------------

host$cwb$blTable()
#> 
#>          x_spline x<<>>site_combine 
#>                28              1972
table(cboost$getSelectedBaselearner())
#> 
#> x_site_tensor      x_spline 
#>          1972            28

cfd = host$cwb$getBlMap()
cfc = cboost$getEstimatedParameter()

all.equal(as.matrix(cfd[["x<<>>site_combine"]]), cfc[["x_site_tensor"]], check.attributes = FALSE)
#> [1] TRUE
all.equal(as.matrix(cfd[["x_spline"]]), cfc[["x_spline"]], check.attributes = FALSE)
#> [1] TRUE

all.equal(cboost$getRiskVector()[-1], host$risk)
#> [1] TRUE
```
