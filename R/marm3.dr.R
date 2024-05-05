
##--------------Estimation with Penalty by CV----------------------##
#' Fit MARM with sparsity assumption and unknown ranks.
#'
#' Fit a multivariate additive model for multi-view data (MARM) using B-splines with unknown ranks (\eqn{r_{1g}, r_{2g}, r_{3g}}). 
#' Multiple third-order coefficient tensors can be estimated by this function. The group sparse penalty such as LASSO, MCP or SCAD 
#' and the coordinate descent algorithm are used to yield a sparsity estimator. The BIC or cross-validation method are used to 
#' search the optimal regularization parameter, multiple ranks and the number of B-spline basis functions simultaneously.
#'
#' @param Y A \eqn{n \times q} numeric matrix of responses.
#' @param X A \eqn{n \times p} numeric design matrix for the model, where \eqn{p = \sum_{g} p_g}.
#' @param group A \eqn{p} vector of the grouping index of predictors, e.g., \eqn{group = c(1,1,1,2,2,2)} means there are 
#'        \eqn{6} predictors in the model, and the first three predictors are in the same group and the last three predictors are 
#'        in another one. By default, we set \eqn{group = rep(1, p)}.
#' @param K_index The user-specified sequence of \code{K}. Default is a length-1 vector \code{6}.
#' @param r1_index A user-specified sequence of \eqn{r_1} values, where \eqn{r_1} is the first dimension of the tensor. 
#'        Default is \code{r1_index}\eqn{ = 1, \dots, \min(\lceil\log(n)\rceil, p)}.
#' @param r2_index A user-specified sequence of \eqn{r_2} values, where \eqn{r_2} is the second dimension of the tensor. 
#'        Default is \code{r2_index = 1, \dots}, \code{max(K_index)}.
#' @param r3_index A user-specified sequence of \eqn{r_3} values, where \eqn{r_3} is the third dimension of the tensor. 
#'        Default is \code{r3_index}\eqn{ = 1, \dots, \min(\lceil\log(n)\rceil, q)}.
#' @param method The method to be applied to select the number of B-spline basis functions, regularization parameters, and 
#'        multiple ranks simultaneously. Either \code{BIC} (default), or \code{CV}.
#' @param ncv The number of cross-validation folds. Default is \code{10}. If \code{method} is not \code{CV}, \code{ncv} is useless.
#' @param penalty The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length \code{nlam} is computed, 
#'        equally spaced on the log scale.
#' @param D0 A user-specified list of initialized values, including \code{ng} sub-lists where ng is the number of groups. 
#'        For each sub-list, it has four initialized matrix \eqn{S_{(3)}} (called \code{S}), \code{A}, \code{B}, and \code{C}. 
#'        By default, a list of initialization satisfying fixed ranks is computed by random.
#' @param intercept A logical value indicating whether the intercept is fitted. Default is \code{TRUE} or set to zero by \code{FALSE}.
#' @param degr The number of knots of B-spline base function. Default is \code{3}.
#' @param nlam The number of lambda values. Default is \code{50}.
#' @param lam_min The smallest value for lambda, as a fraction of lambda.max. Default is \code{0.01}.
#' @param eps Convergence threshold. The algorithm iterates until the relative change in any coefficient is less than \code{eps}. 
#'        Default is \code{1e-4}.
#' @param max_step Maximum number of iterations. Default is \code{20}.
#' @param eps1 Convergence threshold. The Coordinate descent method algorithm iterates until the relative change in any coefficient 
#'        is less than \code{eps1}. Default is \code{1e-4}.
#' @param max_step1 The maximum number of iterates for the coordinate descent method. Default is \code{20}.
#' @param gamma The tuning parameter of the MCP/SCAD penalty.
#' @param dfmax Upper bound for the number of nonzero coefficients. Default is no upper bound. However, for large data sets, 
#'        computational burden may be heavy for models with a large number of nonzero coefficients.
#' @param alpha Tuning parameter for the Mnet estimator which controls the relative contributions from the LASSO, MCP or SCAD 
#'        penalty and the ridge, or L2 penalty. \code{alpha = 1} is equivalent to LASSO, MCP or SCAD penalty, while \code{alpha = 0} 
#'        would be equivalent to ridge regression. However, \code{alpha = 0} is not supported; \code{alpha} may be arbitrarily small, 
#'        but not exactly 0.
#' 
#' @return A list containing the following components:
#' \itemize{
#'   \item{D}{Estimator of coefficients \eqn{D_{(3)} = (D^1_{(3)},...,D^{ng}_{(3)})}.}
#'   \item{mu}{Estimator of intercept \eqn{\mu}.}
#'   \item{S.opt}{A length-\eqn{ng} list including estimator of the core tensor \eqn{S_{(3)}} of each coefficient tensor.}
#'   \item{A.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{A} of each coefficient tensor.}
#'   \item{B.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{B} of each coefficient tensor.}
#'   \item{C.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{C} of each coefficient tensor.}
#'   \item{rk_opt}{The optimal ranks and the number of B-spline basis functions that selected by \code{BIC}, or \code{CV}. 
#'          It is a vector with length 4, which are selected \eqn{r_1}, \eqn{r_2}, \eqn{r_3}, and \eqn{K}.}
#'   \item{lambda.seq}{The sequence of regularization parameter values in the path.}
#'   \item{lambda_opt}{The value of \code{lambda} with the minimum \code{BIC} or \code{CV} value.}
#'   \item{rss}{Residual sum of squares (RSS).}
#'   \item{df}{Degrees of freedom.}
#'   \item{activeX}{The active set of \eqn{X}. A length-\eqn{p} vector.}
#'   \item{opts}{Other related parameters used in algorithm. Some of them are set by default.}
#'   \item{opts_pen}{Other related parameters used in algorithm (especially parameters in penalty). Some of them are set by default.}
#' }
#' 
#' @examples
#' library(comarm)
#' n <- 200; q <- 5; p <- 100; s <- 3; ng = 4
#' group <- rep(1:ng, each = p/ng)
#' mydata <- marm3.sim.fbs(n, q, p, s, group)
#' fit <- with(mydata, marm3.dr(Y, X, group, K, r1_index, r2_index, r3_index, D0 = D0, nlam = 5))
#' 
#' @seealso \code{\link{marm3}}
#' 
#' @useDynLib comarm, .registration = TRUE
#' @export
marm3.dr <- function(Y, X, group = NULL, K_index = NULL, r1_index = NULL, r2_index = NULL,
                     r3_index = NULL, method = "BIC", ncv = 10, penalty = "LASSO",
                     lambda = NULL, D0 = NULL, intercept = TRUE, nlam = 50, degr = 3,
                     lam_min = 0.01, eps = 1e-4, max_step = 20, eps1 = 1e-4,
                     max_step1 = 20, gamma = 2, dfmax = NULL, alpha = 1) {
    n <- nrow(Y)
    q <- ncol(Y)
    nx <- ncol(X)
    if(is.null(group)) group = rep(1,nx)
    gunique <- unique(group)
    G = length(gunique)
    p = rep(0,G)
    for(g in 1:G) p[g] = sum(group==gunique[g])
    K1 <- 6
    if(degr>min(6,K1-1)-1) stop("K must be larger than degree+1 !")
    if(is.null(K_index)) K_index = min(6,K1-1):max(8,K1+1)
    if(is.null(r1_index)) r1_index = 1:min(floor(log(n)),min(p))
    if(is.null(r2_index)) r2_index = 1:min(K_index)
    if(is.null(r3_index)) r3_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MCP penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = nx + 1
    # initial A,B,C,S
    if(is.null(D0)){
      set.seed(1)
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = max(K_index)
      
      B = rbind(diag(r2_max), matrix(0,K_max-r2_max,r2_max))
      C = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
      S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      D0 = list(); 
      for(j in 1:G){
        A = rbind(diag(r1_max), matrix(0,p[j]-r1_max,r1_max))
        SABC = list(S=S,A=A,B=B,C=C)
        D0[[j]] = SABC
      }
    }
    opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=2,r2=2,r3=2,p=p,q=q,degr=degr,K=max(K_index),G=G,nx=nx)
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      setlam = c(1,lam_min,alpha,nlam)
      Sinit = list()
      Ainit = list()
      Binit = list()
      Cinit = list()
      Z = list()
      for(i in 1:G){
        Sinit[[i]] = D0[[i]]$S
        Ainit[[i]] = D0[[i]]$A
        Binit[[i]] = D0[[i]]$B
        Cinit[[i]] = D0[[i]]$C
        Z[[i]] = bsbasefun(X[,group==gunique[i]],max(K_index),degr)
        Zbar1 = colMeans(Z[[i]])
        Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
      }
      Ybar = colMeans(Y)
      Y1 = Y - matrix(rep(Ybar,each=n),n)
      lambda = setuplambda(Y1,Z,Sinit,Ainit,Binit,Cinit,nx,G,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    opts_pen = list(gamma=gamma,dfmax=dfmax,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=1) 
    #---------------- The selection by CV  ---------------------#  
    if(method=="BIC") fit_dr = marm3.bic(Y,X,group,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen)
    if(method=="CV") fit_dr = marm3.cv(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen)
    
    return(fit_dr)
  }