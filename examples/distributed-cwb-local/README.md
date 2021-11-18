
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Distributed CWB prototype

## General design

  - Focus on simple base learner:
      - Splines
      - Categorial base learner via dummy coding and ridge penalty
      - Row wise Kronecker product

## Initialization stage

Since initializing base learner often requires knoledge about the pooled
feature we have to consider: - Splines: - Min. and max. calculation for
the whole feature to calculate equal knots for all sites. - **Shared
data:** Min. and max. value. - Categorical - Ridge: All groups of the
categorical feature to build the 0-1-design matrix. - **Shared data:**
Groups at the site. - For both base learner: The penalty depending on
the degrees of freedom needs to be calculated on the “global” \(X^TX\)
matrix.

After calculating these properties and sharing them with all sites, we
can calculate the design matrices at each site required for modelling.

## Data at the host

In order to estimate the parameters, the penalty is kept at the host to
conduct \((X^TX + \lambda K)^{-1}X^Ty\) with \[
F = \sum_{k=1}^K F_k, \qquad F_k = X_k^TX_k
\]
<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=F = \sum_{k=1}^K F_k, \qquad F_k = X_k^TX_k" title="F = \sum_{k=1}^K F_k, \qquad F_k = X_k^TX_k" />
<img src="https://latex.codecogs.com/svg.latex?\Large&space;x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />
<img src="https://render.githubusercontent.com/render/math?math=x_{1,2} = \frac{-b \pm \sqrt{b^2-4ac}}{2b}">

\[
s = \sum_{k=1}^K s_k, \qquad s_k = X_k^Ty_k
\]

## Fitting stage

  -
