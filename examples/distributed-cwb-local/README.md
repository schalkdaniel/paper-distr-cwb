
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Distributed CWB prototype

## General design

  - Focus on simple base learner:
      - Splines
      - Categorial base learner via dummy coding and ridge penalty
      - Row wise Kronecker product
  - `R6` classes for sites and host to “simulate” communication.
      - Sites contains their data.
      - The host does not contain any data, but information about the
        structure of the base learners, the global CWB model, etc.
  - CWB is also a `R6` class as well as the base learners.

<!-- end list -->

``` r
library(R6)

source(here::here("R/baselearner.R"))
source(here::here("R/cwb.R"))
source(here::here("R/site.R"))
source(here::here("R/host.R"))

nsim = 1000L
x = runif(nsim, 0, 10)
svec = sample(LETTERS[1:4], nsim, TRUE)
dat = data.frame(x = x, y = rnorm(nsim, sin(x), 0.2))

# CWB structure:
bl_spline  = BaselearnerPSpline$new("x", ord = 2)
#bl_ridge   =
#bl_combine =

# Define base learner list
#bls = list(bl_spline, bl_combine)
bls = list(bl_spline)

# Local CWB:
loss  = function(y, y_pred) 0.5 * (y - y_pred)^2
dloss = function(y, y_pred) -(y - y_pred)

cwb = CWB$new("y", lr = 0.1, loss = loss, dloss = dloss)

# Define sites:
sites = lapply(unique(svec), function(s) {
  Site$new(s, dat[svec == s, ], cwb, bls) })

# Define the host:
host = Host$new(cwb, bls, sites)
```

## Initialization stage

Since initializing base learner often requires knowledge about the
pooled feature we have to consider: - Splines: - Min. and max.
calculation for the whole feature to calculate equal knots for all
sites. - **Shared data:** Min. and max. value. - Categorical - Ridge:
All groups of the categorical feature to build the 0-1-design matrix. -
**Shared data:** Groups at the site. - For both base learner: The
penalty depending on the degrees of freedom needs to be calculated on
the “global” \(X^TX\) matrix.

After calculating these properties and sharing them with all sites, we
can calculate the design matrices at each site required for modelling.

``` r
host$initializePrediction()
host$initializeBaselearner()
```

## Data at the host

In order to estimate the parameters, the penalty \(\lambda\) and penalty
matrix \(K_j\) per base learner \(j\) to conduct
\((F_j + \lambda_j K_j)^{-1}s\) with is kept at the host. Additionally,
\(F_j = X_j^TX_j\), which is also initialized once in the beginning, is
kept at the host. \[
F_j = \sum_{k=1}^K F_{k,j}, \qquad F_{k,j} = X_{k,j}^TX_{k,j}
\]

``` r
## Aggregated static parts:
str(host$init_aggr)
#> List of 1
#>  $ x_spline:List of 1
#>   ..$ knots: num [1:16, 1] -1.80855 -0.90055 0.00744 0.91544 1.82343 ...
str(host$bl_parts)
#> List of 1
#>  $ x_spline:List of 4
#>   ..$ XtX: num [1:13, 1:13] 6.537 13.145 0.869 0 0 ...
#>   ..$ K  : num [1:13, 1:13] 1 -2 1 0 0 0 0 0 0 0 ...
#>   ..$ df : num 4
#>   ..$ pen: num 583

host$offset
#> [1] 0.1556
```

## Fitting stage

The fitting stage only requires the calculation of \(X_{k,j}^Tu\) with
pseudo residuals \(u\) in each iteration. For each site \(k\) and base
learner \(j\), \(X^T_{k,j} u_k\) is calculated and used to find the best
base learner. -\[
s_j = \sum_{k=1}^K s_{k,j}, \qquad s_k = X_{k,j}^Tu_k
\]

**Note:** Finding the best base learner requires multiple communication
steps: - Communicate \(X_{k,j}^T u_k\) - Send parameters
\(\hat{\theta}_j = F_j^{-1} u_j\) - Communicate sum of squared errors
\(\sum_k \| u_k - X_{k,j}\hat{\theta}_j\) - Update via sending best
parameter \(\hat{\theta}_{j^\ast}\), \(j^\ast = \argmin_j SSE_j\)

``` r
mstop = 1000L
trace = character()
for (m in seq_len(mstop)) {
  trace = c(trace, capture.output(host$updateCWB(m)))
}
cat(tail(trace), sep = "")
#> 995: 0.01941996: 0.01941997: 0.01941998: 0.01941999: 0.019411000: 0.01941
```

## Check the algorithm

### Estimated effect

The CWB object of the host contains the final model as base learner map
`$getBlMap()` while the whole trace of added base learner is accessible
via `getBlTrace()`:

``` r
str(host$cwb$getBlTrace()[1:3])
#> List of 3
#>  $ :List of 2
#>   ..$ name : chr "x_spline"
#>   ..$ param: num [1:13, 1] 0.6549 0.4755 0.28 0.0332 -0.2373 ...
#>  $ :List of 2
#>   ..$ name : chr "x_spline"
#>   ..$ param: num [1:13, 1] 0.592 0.4346 0.2609 0.0346 -0.2191 ...
#>  $ :List of 2
#>   ..$ name : chr "x_spline"
#>   ..$ param: num [1:13, 1] 0.5344 0.3971 0.2434 0.0359 -0.2024 ...
str(host$cwb$getBlMap())
#> List of 1
#>  $ x_spline: num [1:13, 1] -0.5 0.31 0.88 0.673 -0.195 ...
```

Visualizing the estimated effect reveals an accurate effect estimation:

``` r
## Copy everything:
blst = bls
d0 = dat

lapply(blst, function(bl) bl$initDesign(d0, host$init_aggr$x_spline$knots))
#> [[1]]
#> [1] TRUE
d0$pred = host$offset + blst[[1]]$linPred(host$cwb$getBlMap()[["x_spline"]])

library(ggplot2)

ggplot(d0) +
  geom_point(aes(x = x, y = y, color = "Truth"), alpha = 0.2) +
  geom_line(aes(x = x, y = pred, color = "Predicted"))
```

<img src="Readme_files/unnamed-chunk-7-1.png" width="60%" style="display: block; margin: auto;" />

### Comparison with compboost

A comparison with compboost also shows the correctness of the
distributed algorithm:

``` r
## Compare to compboost:
devtools::load_all("~/repos/compboost")
cboost = boostSplines(data = dat, target = "y", loss = LossQuadratic$new(),
  learning_rate = 0.1, iterations = 1000L, n_knots = 10, df = 4, degree = 2,
  differences = 2, trace = 200L)
#>    1/1000   risk = 0.24  time = 0   
#>  200/1000   risk = 0.02  time = 10793   
#>  400/1000   risk = 0.02  time = 22488   
#>  600/1000   risk = 0.019  time = 35395   
#>  800/1000   risk = 0.019  time = 49378   
#> 1000/1000   risk = 0.019  time = 64591   
#> 
#> 
#> Train 1000 iterations in 0 Seconds.
#> Final risk based on the train set: 0.019

all.equal(cboost$baselearner_list$x_spline$factory$getPenaltyMat(),
  host$bl_parts$x_spline$K)
#> [1] TRUE

all.equal(cboost$baselearner_list$x_spline$factory$getPenalty(),
  host$bl_parts$x_spline$pen)
#> [1] TRUE

all.equal(check.attributes = FALSE,
  cboost$baselearner_list$x_spline$factory$getData() %*%
    t(cboost$baselearner_list$x_spline$factory$getData()),
  host$bl_parts$x_spline$XtX)
#> [1] TRUE

cf_cboost = cboost$getEstimatedCoef()
cf_dist   = host$cwb$getBlMap()

all.equal(cf_cboost$offset, host$offset)
#> [1] TRUE
all.equal(cf_cboost$x_spline, cf_dist$x_spline, check.attributes = FALSE)
#> [1] TRUE
all.equal(cboost$model$getRiskVector()[-1], host$risk)
#> [1] TRUE
```
