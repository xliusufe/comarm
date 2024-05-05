#' Fit MARM with sparsity assumption and fixed ranks.
#'
#' Fit a multivariate additive model for multi-view data (MARM) using B-splines with given ranks (\eqn{r_{1g}, r_{2g}, r_{3g}}). 
#' Multiple third-order coefficient tensors can be estimated by this function. The group sparse penalty such as LASSO, MCP or SCAD 
#' and the coordinate descent algorithm are used to yield a sparsity estimator. The BIC or cross-validation method are used to 
#' search the optimal regularization parameter.
#'
#' @param Y A \eqn{n \times q} numeric matrix of responses.
#' @param X A \eqn{n \times p} numeric design matrix for the model, where \eqn{p = \sum_{g} p_g}.
#' @param group A \eqn{p} vector of the grouping index of predictors, e.g., \eqn{group = c(1,1,1,2,2,2)} means there are 
#'        \eqn{6} predictors in the model, and the first three predictors are in the same group and the last three predictors are 
#'        in another one. By default, we set \eqn{group = rep(1, p)}.
#' @param K The number of B-spline basis functions, that is the sum of both degrees of basis functions and the number of knots. 
#'        Default is \code{6}, which means cubic spline.
#' @param r1 The first dimension of the singular value matrix of the tensor. Default is \code{2}.
#' @param r2 The second dimension of the singular value matrix of the tensor. Default is \code{2}.
#' @param r3 The third dimension of the singular value matrix of the tensor. Default is \code{2}.
#' @param method The method to be applied to select regularization parameters. Either \code{BIC} (default), or \code{CV}.
#' @param ncv The number of cross-validation folds. Default is \code{10}. If \code{method} is not \code{CV}, \code{ncv} is useless.
#' @param penalty The penalty to be applied to the model. Either \code{LASSO} (the default), \code{MCP} or \code{SCAD}.
#' @param lambda A user-specified sequence of lambda values. By default, a sequence of values of length \code{nlam} is computed, 
#'        equally spaced on the log scale.
#' @param D0 A user-specified list of initialized values, including \code{ng} sub-lists where ng is the number of groups. 
#'        For each sub-list, it has four initialized matrices \code{S_{(3)}} (called \code{S}), \code{A}, \code{B}, and \code{C}. 
#'        By default, a list of initialization satisfying fixed ranks is computed by random.
#' @param intercept A logical value indicating whether the intercept is fitted. Default is \code{TRUE} or set to zero by \code{FALSE}.
#' @param degr The number of knots of B-spline base function. Default is \code{3}.
#' @param nlam The number of lambda values. Default is \code{20}.
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
#'   \item{D}{Estimator of coefficients \eqn{D_{(3)} = (D^1_{(3)}, ..., D^{ng}_{(3)})} where \eqn{ng} is the number of groups.}
#'   \item{mu}{Estimator of intercept \eqn{\mu}.}
#'   \item{S.opt}{A length-\eqn{ng} list including estimator of the core tensor \eqn{S_{(3)}} of each coefficient tensor.}
#'   \item{A.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{A} of each coefficient tensor.}
#'   \item{B.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{B} of each coefficient tensor.}
#'   \item{C.opt}{A length-\eqn{ng} list including estimator of the factor matrix \eqn{C} of each coefficient tensor.}
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
#' n <- 200; q <- 5; p <- 100; s <- 3; ng <- 4
#' group <- rep(1:ng, each = p/ng)
#' mydata <- marm3.sim.fbs(n, q, p, s, group, isfixedR = 1)
#' fit <- with(mydata, marm3(Y, X, group, K, r10, r20, r30, D0 = D0, nlam = 5))
#' 
#' @seealso \code{\link{marm3.dr}}
#' 
#' @useDynLib comarm, .registration = TRUE
#' @export
marm3 <- function(Y, X, group = NULL, K = 6, r1 = NULL, r2 = NULL, r3 = NULL,
                  method = "BIC", ncv = 10, penalty = "LASSO", lambda = NULL, D0 = NULL,
                  intercept = TRUE, degr = 3, nlam = 20, lam_min = 0.01,
                  eps = 1e-4, max_step = 20, eps1 = 1e-4, max_step1 = 20,
                  gamma = 2, dfmax = NULL, alpha = 1) {
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    nx <- ncol(X)
    if(is.null(group)) group = rep(1,nx)
    gunique <- unique(group)
    G = length(gunique)
    p = rep(0,G)
    for(g in 1:G) p[g] = sum(group==gunique[g])
    if(degr>K-1) stop("K must be larger than degree+1 !")
    if(is.null(r1)) r1 <- 2 
    if(is.null(r2)) r2 <- 2
    if(is.null(r3)) r3 <- 2
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = nx + 1
    # initial A,B,C,S
    if(is.null(D0)){
      set.seed(1)
      B = rbind(diag(r2), matrix(0,K-r2,r2))
      C = rbind(diag(r3), matrix(0,q-r3,r3))
      S = matrix(rnorm(r1*r2*r3),r3,r1*r2)
      D0 = list(); 
      for(j in 1:G){
        A = rbind(diag(r1), matrix(0,p[j]-r1,r1))
        SABC = list(S=S,A=A,B=B,C=C)
        D0[[j]] = SABC
      }
    }
    opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=r1,r2=r2,r3=r3,p=p,q=q,degr=degr,K=K,G=G,nx=nx)
    Sinit = list()
    Ainit = list()
    Binit = list()
    Cinit = list()
    Z = list()
    Zbar = NULL
    for(i in 1:G){
      Sinit[[i]] = D0[[i]]$S[1:r3,1:(r1*r2)]
      Ainit[[i]] = D0[[i]]$A[,1:r1]
      Binit[[i]] = D0[[i]]$B[,1:r2]
      Cinit[[i]] = D0[[i]]$C[,1:r3]
      if(dim(Sinit[[i]])[1]!=dim(Cinit[[i]])[2]) Sinit[[i]] = t(Sinit[[i]])
      Z[[i]] = as.matrix(bsbasefun(X[,group==gunique[i]],K,degr))
      Zbar1 = colMeans(Z[[i]])
      Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
      Zbar = cbind(Zbar, Zbar1)
    }
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      setlam = c(1,lam_min,alpha,nlam)
      lambda = setuplambda(Y1,Z,Sinit,Ainit,Binit,Cinit,nx,G,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    opts_pen = list(gamma=gamma,dfmax=dfmax,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=1) 
    #---------------- The selection by BIC or CV  ---------------------# 
    if(method=="BIC"){
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
      #df = 0; for(g in 1:G) df = df + r1*r2*r3 + fit$df[g,]*r1 + K*r2 + q*r3 - r1^2 - r2^2 - r3^2
      df = NULL
      for (ll in 1:nlam) {
        df1 = 0
        for(g in 1:G) {
          if(fit$df[g,ll]<r1) {
            df1 = df1 + r1*r2*r3 + fit$df[g,ll]*r1 + K*r2 + q*r3 - r2^2 - r3^2
          }else {
            df1 = df1 + r1*r2*r3 + fit$df[g,ll]*r1 + K*r2 + q*r3 - r1^2 - r2^2 - r3^2
          }
        }
        df = c(df, df1)
      }
      bic = log(fit$likhd/(n*q)) + log(n*q)*df/(n*q)
      
      selected = which.min(bic)
      lambda_opt = lambda[selected]

      opts_pen$nlam = length(lambda[1:selected])
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda[1:selected],opts,opts_pen)

      activeX = fit$betapath[,selected]
      Snew = fit$Snew[[selected]]
      Anew = fit$Anew[[selected]]
      Bnew = fit$Bnew[[selected]]
      Cnew = fit$Cnew[[selected]]
      Dn = NULL
      for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
      if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam>1){
      len_cv = ceiling(n/ncv)
      RSS = rep(0,nlam)
      for(jj in 1:ncv){
        cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
        if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
        Ytrain = Y1[-cv.id,]
        Xtrain = X[-cv.id,]
        Ytest = Y1[cv.id,]
        Xtest = X[cv.id,]
        Ztrain = list()
        Ztest = list()
        for(i in 1:G){
          Ztrain[[i]] = bsbasefun(Xtrain[,group==gunique[i]],K,degr)
          Ztest[[i]] = bsbasefun(Xtest[,group==gunique[i]],K,degr)
        } 
        fit = EstPenColumnCV(Ytrain,Ztrain,Ytest,Ztest,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
        RSS = RSS + fit$likhd
      } 
      selected = which.min(RSS)
      lambda_opt = lambda[selected]
      
      opts_pen$nlam = length(lambda[1:selected])
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda[1:selected],opts,opts_pen)
      
      activeX = fit$betapath[,selected]
      Snew = fit$Snew[[selected]]
      Anew = fit$Anew[[selected]]
      Bnew = fit$Bnew[[selected]]
      Cnew = fit$Cnew[[selected]]
      Dn = NULL
      for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
      if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam==1){
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
      selected = 1
      lambda_opt = lambda
      activeX = fit$betapath
      Snew = fit$Snew[[selected]]
      Anew = fit$Anew[[selected]]
      Bnew = fit$Bnew[[selected]]
      Cnew = fit$Cnew[[selected]]
      Dn = NULL
      for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
      if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
      else mu = rep(0,q)
    }

    return(list(D          = Dn,
                mu         = mu,
                S.opt      = Snew,
                A.opt      = Anew,
                B.opt      = Bnew,
                C.opt      = Cnew,
                lambda.seq = lambda,
                lambda.opt = lambda_opt,
                rss        = fit$likhd[selected],
                df         = fit$df,
                activeX    = activeX,
                opts       = opts,
                opts_pen   = opts_pen))
  }