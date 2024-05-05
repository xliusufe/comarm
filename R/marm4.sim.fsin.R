##------------------------------------------------------------------------##
##-------------generate scenario II data in structural MARM--------------##
# e.g. mydata=marm4.sim.fsin(200,5,100,10,rep(1:4,each=25))
#' Generate scenario II data from structural MARM model.
#'
#' Generate scenario II data for a structural multivariate additive model for multi-view data (MARM).
#' This function simulates data based on sinusoidal transformations of predictors, suitable for 
#' analysis with a structural MARM model using fixed or varying tensor ranks.
#'
#' @param n Sample size.
#' @param q The number of responses, \eqn{q\geq1}.
#' @param p The number of covariates, \eqn{p\geq1}.
#' @param s The true covariates of each view associating with responses, \eqn{s\geq1}.
#' @param group A vector of grouping index for predictors, splitting predictors into groups. By default, it splits 
#'        the predictors into two groups.
#' @param r10 The first dimension of the tensor. Default is \code{2}.
#' @param r20 The second dimension of the tensor. Default is \code{2}.
#' @param r30 The third dimension of the tensor. Default is \code{2}.
#' @param r40 The fourth dimension of the tensor. Default is \code{2}.
#' @param isfixedR Indicates if the tensor ranks are fixed (\code{TRUE} or \code{FALSE}).
#' @param D42 Mode of unfolding \eqn{D_{(2)}}. Generated randomly by default.
#' @param K Number of B-spline basis functions, typically \code{6} for cubic splines.
#' @param degr Number of knots in the B-spline base function. Default is \code{3}.
#' @param sigma2 Error variance. Default is \code{0.1}.
#' @param seed_id Seed for random number generation, for reproducibility. Default is \code{1000}.
#' @param r1_index User-specified sequence of \eqn{r_1} values, ignored if \code{isfixedR} is \code{TRUE}.
#' @param r2_index User-specified sequence of \eqn{r_2} values, ignored if \code{isfixedR} is \code{TRUE}.
#' @param r3_index User-specified sequence of \eqn{r_3} values, ignored if \code{isfixedR} is \code{TRUE}.
#' @param r4_index User-specified sequence of \eqn{r_4} values, ignored if \code{isfixedR} is \code{TRUE}.
#' @param D0 List of initialization values for matrices \eqn{S_{(4)}}, \code{A}, \code{B}, \code{C}, and \code{D}.
#'
#' @return A list containing generated data elements:
#'   \item{Y}{Response matrix, dimensions \eqn{n \times q}.}
#'   \item{X}{Design matrix, dimensions \eqn{n \times p}.}
#'   \item{f0}{Matrix of true functions, dimensions \eqn{n \times p}.}
#'   \item{group}{Vector indicating the grouping index of predictors, length \eqn{p}.}
#'   \item{D0}{Initialized values used in tensor decomposition.}
#'   \item{...}{Additional options for the algorithm.}
#'
#' @examples
#' library(comarm)
#' n <- 200; q <- 5; p <- 100; s <- 3; ng = 4
#' group <- rep(1:ng, each = p/ng)
#' mydata <- marm4.sim.fsin(n, q, p, s, group)
#'
#' @seealso \code{\link{marm4.sim.fbs}}
#' @useDynLib comarm, .registration = TRUE
#' @export
marm4.sim.fsin <- function(n, q, p, s, group = NULL, r10 = 2, r20 = 2, r30 = 2, r40 = 2,
                            isfixedR = 0, D42 = NULL, K = 6, degr = 3, sigma2 = NULL,
                            seed_id = NULL, r1_index = NULL, r2_index = NULL,
                            r3_index = NULL, r4_index = NULL, D0 = NULL) {
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(ng<r30) stop("ng must be not smaller than r30")
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D42)) {
    set.seed(2)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(ng*r30),nrow = ng)
    U3 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r40),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r10*r20*r30*r40,5,10),nrow = r30)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
    D42 = TransferModalUnfoldingsT(D44,4,2,c(s,K,ng,q))
  }
  
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  id = NULL
  for(j in 1:q) id = c(id, c(1:s)+(j-1)*s*ng)
  for(j in 1:ng){
    D2 = D42[,id+(j-1)*s]
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
  if(is.null(r3_index)) r3_index = 1:min(4,ng)
  if(is.null(r4_index)) r4_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r30),ng,r30)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r40),q,r40)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30*r40),r40,r10*r20*r30)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      r4_max = max(r4_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r3_max),ng,r3_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_max),q,r4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max*r4_max),r4_max,r1_max*r2_max*r3_max)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D44=D44,D42=D42,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,r40=r40,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,r4_index=r4_index,D0=D0))
}

