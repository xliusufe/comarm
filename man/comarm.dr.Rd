% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comarm.dr.R
\name{comarm.dr}
\alias{comarm.dr}
\title{Fit composed model (CoMARM) with sparsity assumption and unknown ranks.}
\usage{
comarm.dr(
  Y,
  X1,
  X2,
  G1 = NULL,
  group = NULL,
  is.fabs = 1,
  K_index = NULL,
  r_index = NULL,
  method = "BIC",
  ncv = 10,
  penalty = "LASSO",
  lambda = NULL,
  D0_t3 = NULL,
  D0_t4 = NULL,
  intercept = TRUE,
  nlam = 50,
  degr = 3,
  lam_min = 0.01,
  eps = 1e-04,
  max_step = 20,
  eps1 = 1e-04,
  max_step1_t3 = 20,
  max_step1_t4 = 20,
  gamma = 2,
  dfmax1 = NULL,
  dfmax2 = NULL,
  alpha = 1,
  vnorm_ratio = 1
)
}
\arguments{
\item{Y}{A \eqn{n\times q} numeric matrix of responses.}

\item{X1}{A \eqn{n\times p1} numeric design matrix for the model, corresponding to first set of views.}

\item{X2}{A \eqn{n\times p2} numeric design matrix for the model, corresponding to second set of views.}

\item{G1}{The number of views in the first set that we consider without intergroup correlation.}

\item{group}{A vector of the grouping index for predictors, usually set to \code{rep(1, p)} by default.}

\item{is.fabs}{A logical value indicating whether data comes from scenario I.}

\item{K_index}{User-specified sequence of K, usually \code{6} for cubic splines.}

\item{r_index}{A sequence of rank values for the tensors.}

\item{method}{Method for parameter selection, either \code{BIC} or \code{CV}.}

\item{ncv}{Number of cross-validation folds.}

\item{penalty}{Type of penalty applied, options: \code{LASSO}, \code{MCP}, \code{SCAD}.}

\item{lambda}{Sequence of lambda values for regularization.}

\item{D0_t3}{Initial values for third-order tensor decomposition.}

\item{D0_t4}{Initial values for fourth-order tensor decomposition.}

\item{intercept}{Whether to fit an intercept.}

\item{nlam}{Number of lambda values to use.}

\item{degr}{Number of knots in B-spline.}

\item{lam_min}{Minimum lambda as a fraction of the maximum.}

\item{eps}{Convergence threshold.}

\item{max_step}{Maximum number of iterations allowed.}

\item{eps1}{Convergence threshold for coordinate descent.}

\item{max_step1_t3}{Maximum iterations for coordinate descent on third-order tensor.}

\item{max_step1_t4}{Maximum iterations for coordinate descent on fourth-order tensor.}

\item{gamma}{Tuning parameter for MCP/SCAD.}

\item{dfmax1}{Maximum number of non-zero coefficients for third-order tensor.}

\item{dfmax2}{Maximum number of non-zero coefficients for fourth-order tensor.}

\item{alpha}{Balancing parameter between LASSO/MCP/SCAD and ridge.}

\item{vnorm_ratio}{Ratio between lambda values for different tensor types.}
}
\value{
List of model outputs including estimates of tensors, coefficients, and model diagnostics.
}
\description{
Fit a composed model using B-splines with unknown ranks (\eqn{r_{1g}, r_{2g}, r_{3g}} and \eqn{r_{1}, r_{2}, r_{3}, r_{4}}). 
This function estimates multiple third-order and fourth-order coefficient tensors. It utilizes group sparse 
penalties such as LASSO, MCP, or SCAD, and a coordinate descent algorithm to yield a sparsity estimator. 
The BIC or cross-validation method is applied to optimize regularization parameters, ranks, and the 
number of B-spline basis functions simultaneously.
}
\examples{
library(comarm)
n <- 200; q <- 5; p <- 20; s1 <- 5; s2 <- 3; G1 <- 1; ng = 4
group <- rep(1:ng, each = p/ng)
mydata <- comarm.sim.fbs(n, q, p, s1, s2, G1, group)
fit <- with(mydata, comarm.dr(Y = Y, X1 = X1, X2 = X2, G1 = G1,
            group = group, is.fabs = is.fabs, K_index = K,
            r_index = r_index, D0_t3 = D0_t3, D0_t4 = D0_t4, nlam = 5))

}
\seealso{
\code{\link{comarm}}
}
