##------------------------------------------------------------------------##
##----------------generate scenario II data in MARM----------------------##
# e.g. mydata=marm3.sim.fsin(200,5,100,10,rep(1:4,each=25))
#' Generate scenario II data from MARM model.
#'
#' Generate scenario II data for a multivariate additive model for multi-view data. This scenario involves
#' generating data based on sinusoidal transformations of the predictors, incorporating multivariate responses
#' and a structured group sparsity across tensor decompositions.
#'
#' @param n Sample size.
#' @param q The number of responses, \eqn{q \geq 1}.
#' @param p The number of covariates, \eqn{p \geq 1}.
#' @param s The true covariates of each view associated with responses, \eqn{s \geq 1}.
#' @param group A length-\eqn{p} vector of the grouping index of predictors, typically set as \code{group = rep(1, p)}.
#'        Default groups \eqn{group = c(1, 1, 1, 2, 2, 2)}, which partitions the predictors into two groups.
#' @param r10 The first dimension of the tensor. Default is \code{2}.
#' @param r20 The second dimension of the tensor. Default is \code{2}.
#' @param r30 The third dimension of the tensor. Default is \code{2}.
#' @param isfixedR Logical value indicating if ranks are fixed.
#' @param D2 Mode of unfolding \eqn{D_{(2)}}, generated randomly by default.
#' @param K The number of B-spline basis functions, combining degrees of basis functions and the number of knots.
#'        Default is \code{6}, indicating cubic splines.
#' @param degr The number of knots in the B-spline base function. Default is \code{3}.
#' @param sigma2 Error variance. Default is \code{0.1}.
#' @param seed_id Seed for random number generation, affecting reproducibility. Default is \code{1000}.
#' @param r1_index User-specified sequence of \eqn{r_1} values. 
#'        Default is \code{r1_index}\eqn{ = 1, \dots, \min(\lceil\log(n)\rceil, q)}.
#'        Ignored if \code{isfixedR = 1}.
#' @param r2_index User-specified sequence of \eqn{r_2} values. Default sequence spans from 1 to the maximum of \code{K}.
#'        Ignored if \code{isfixedR = 1}.
#' @param r3_index User-specified sequence of \eqn{r_3} values. 
#'        Default is \code{r3_index}\eqn{ = 1, \dots, \min(\lceil\log(n)\rceil, q)}.
#'        Ignored if \code{isfixedR = 1}.
#' @param D0 List of initialization values for each group, with each containing initialized matrices \eqn{S_{(3)}}, \code{A},
#'        \code{B}, and \code{C}. Randomly computed by default, satisfying fixed ranks.
#'
#' @return A list with the following components:
#'   \item{Y}{Response matrix, dimensions \eqn{n \times q}.}
#'   \item{X}{Design matrix, dimensions \eqn{n \times p}.}
#'   \item{f0}{Matrix of true functions, dimensions \eqn{n \times p}.}
#'   \item{group}{Vector indicating the grouping index of predictors, length \eqn{p}.}
#'   \item{D0}{Initialized values used in tensor decomposition.}
#'   \item{...}{Other algorithmic options, dynamically specified.}
#'
#' @examples
#' library(comarm)
#' n <- 200; q <- 5; p <- 100; s <- 3; ng <- 4
#' group <- rep(1:ng, each = p/ng)
#' mydata <- marm3.sim.fsin(n, q, p, s, group)
#'
#' @seealso \code{\link{marm3.sim.fbs}}
#' @useDynLib comarm, .registration = TRUE
#' @export
marm3.sim.fsin <- function(n, q, p, s, group = NULL, r10 = 2, r20 = 2, r30 = 2, isfixedR = 0,
                            D2 = NULL, K = 6, degr = 3, sigma2 = NULL, seed_id = NULL,
                            r1_index = NULL, r2_index = NULL, r3_index = NULL, D0 = NULL) {
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D2)) {
    set.seed(2)
    S3 <- matrix(runif(r10*r20*r30,5,10),nrow = r30)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r30),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))
    D2 <- TransferModalUnfoldingsT(D3,3,2,c(s,K,q))
  }
  
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    f0 <- f0 + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + sigma2*eps
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r30),q,r30)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30),r30,r10*r20)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_max),q,r3_max)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D3=D3,D2=D2,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,D0=D0))
}