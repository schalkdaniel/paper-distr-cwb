library(R6)

source(here::here("R/baselearner.R"))
source(here::here("R/cwb.R"))
source(here::here("R/site.R"))

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

# First Aggregation:
init = lapply(sites, function(site) site$communicateInit())
str(init)
blns = sapply(bls, function(bl) bl$getName())
init_agg = lapply(blns, function(bln) {
  sinit = lapply(init, function(i) i[[bln]])
  if (grepl("spline", bln)) {
    smin = min(sapply(sinit, function(s) s$min))
    smax = max(sapply(sinit, function(s) s$max))
    out = compboostSplines::createKnots(c(smin, smax), n_knots = 10L, degree = 3)
  }
  out = list(knots = out)
  return(out)
})
names(init_agg) = blns

# Initialize base learner:
nuisance = lapply(sites, function(s) s$initBaselearner(init_agg))

# Initialize prediction:
ss = sapply(sites, function(s) s$getAggregatedInitialization(function(y) mean(y)))
sw = sapply(sites, function(s) s$nrow())
offset = weighted.mean(ss, sw)

nuisance = lapply(sites, function(s) s$initPrediction(offset))


mstop = 20L
for (m in seq_len(mstop)) {
  # Calculate params
  receive = lapply(sites, function(s) s$communicate())
  blparams = lapply(blns, function(bln) {
    blp = lapply(receive, function(r) r[[bln]] )
    XtX = Reduce("+", lapply(blp, function(b) b$XtX))
    Xtu = Reduce("+", lapply(blp, function(b) b$Xtu))
    pK = blp[[1]]$K * blp[[1]]$pen
    return(solve(XtX + pK) %*% Xtu)
  })
  names(blparams) = blns

  # Get SSEs:
  sses = sapply(blns, function(bln) {
    bp = blparams[[bln]]
    sum(sapply(sites, function(s) {
      s$ssePrUpdate(bln, bp)
    }))
  })
  names(sses) = blns

  bl_best = which.min(sses)

  # Site models:
  nuisance = lapply(sites, function(s) {
    s$update(names(bl_best), blparams[[names(bl_best)]])
  })
  # Host model:
  cwb$update(names(bl_best), blparams[[names(bl_best)]])

  rr = weighted.mean(sapply(sites, function(s) s$communicateRisk()), sw)
  cat(m, "/", mstop, ": ", rr, "\n", sep = "")
}
