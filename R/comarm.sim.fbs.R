

##------------------------------------------------------------------------##
##--------------generate scenario I data in composed model----------------##
# e.g. mydata=comarm.sim.fbs(200,5,100,15,10,1,rep(1:4,each=25))
#' Generate scenario I data from composed model (CoMARM).
#'
#' This function generates data for scenario I of a composed model (CoMARM), which uses B-splines 
#' to model data with known tensor ranks for multiple third-order and fourth-order coefficient tensors. 
#' It's suitable for simulations where ranks are fixed and known a priori, allowing detailed study 
#' of model behavior under controlled conditions.
#'
#' @param n Sample size.
#' @param q Number of responses.
#' @param p Total number of covariates.
#' @param s1 Number of true covariates associated with the first set of views.
#' @param s2 Number of true covariates associated with the second set of views.
#' @param G1 Number of views without intergroup correlations.
#' @param group Grouping index for predictors.
#' @param r310 First rank dimension for third-order tensors.
#' @param r320 Second rank dimension for third-order tensors.
#' @param r330 Third rank dimension for third-order tensors.
#' @param r410 First rank dimension for fourth-order tensors.
#' @param r420 Second rank dimension for fourth-order tensors.
#' @param r430 Third rank dimension for fourth-order tensors.
#' @param r440 Fourth rank dimension for fourth-order tensors.
#' @param isfixedR Logical indicating if ranks are fixed.
#' @param D3 Initial data for third-order tensor simulation.
#' @param D44 Initial data for fourth-order tensor simulation.
#' @param K Number of B-spline basis functions.
#' @param degr Number of knots in B-spline.
#' @param sigma2 Error variance.
#' @param seed_id Seed for random number generation.
#' @param r1_t3_index Rank indices for the first dimension of third-order tensors.
#' @param r2_t3_index Rank indices for the second dimension of third-order tensors.
#' @param r3_t3_index Rank indices for the third dimension of third-order tensors.
#' @param r1_t4_index Rank indices for the first dimension of fourth-order tensors.
#' @param r2_t4_index Rank indices for the second dimension of fourth-order tensors.
#' @param r3_t4_index Rank indices for the third dimension of fourth-order tensors.
#' @param r4_t4_index Rank indices for the fourth dimension of fourth-order tensors.
#' @param D0_t3 Initial third-order tensor decomposition.
#' @param D0_t4 Initial fourth-order tensor decomposition.
#'
#' @return A list with the following components:
#'   \item{Y}{Response matrix, dimensions \eqn{n \times q}.}
#'   \item{X1}{First set of design matrices, dimensions \eqn{n \times p1}.}
#'   \item{X2}{Second set of design matrices, dimensions \eqn{n \times p2}.}
#'   \item{f01}{True function values for the first set of views.}
#'   \item{f02}{True function values for the second set of views.}
#'   \item{group}{Grouping index of predictors.}
#'   \item{D0_t3}{Initialized values for the third-order tensor.}
#'   \item{D0_t4}{Initialized values for the fourth-order tensor.}
#'   \item{...}{Additional algorithm options.}
#'
#' @examples
#' library(comarm)
#' n <- 200; q <- 5; p <- 100; s1 <- 5; s2 <- 3; G1 <- 1; ng = 4
#' group <- rep(1:ng, each = p/ng)
#' mydata <- comarm.sim.fbs(n, q, p, s1, s2, G1, group)
#'
#' @seealso \code{\link{comarm.sim.fsin}}
#' @useDynLib comarm, .registration = TRUE
#' @export
comarm.sim.fbs <- function(n, q, p, s1, s2, G1 = NULL, group = NULL, r310 = 2, r320 = 2, r330 = 2, 
                            r410 = 2, r420 = 2, r430 = 2, r440 = 2, isfixedR = 0, D3 = NULL, 
                            D44 = NULL, K = 6, degr = 3, sigma2 = NULL, seed_id = NULL, 
                            r1_t3_index = NULL, r2_t3_index = NULL, r3_t3_index = NULL, 
                            r1_t4_index = NULL, r2_t4_index = NULL, r3_t4_index = NULL, r4_t4_index = NULL, D0_t3 = NULL, D0_t4 = NULL) {
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s1<1) stop("s1 must be not smaller than 1")
  if(s2<1) stop("s2 must be not smaller than 1")
  if(is.null(G1)) stop("the parameter G1 must be entered")
  if(is.null(group)) stop("the parameter group must be entered")
  gunique <- unique(group)
  G = length(gunique)
  G2 = G-G1
  pg1 = sum(group==gunique[G1])
  pg2 = sum(group==gunique[G1+1])
  p1 = pg1*G1
  p2 = pg2*G2
  group1 = group[which(group<=G1)]
  group2 = group[which(group>G1)]-G1
  gunique1 <- unique(group1)
  gunique2 <- unique(group2)
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D3) | is.null(D44)) {
    set.seed(2)
    S3 <- matrix(runif(r310*r320*r330,5,10),nrow = r330)
    T1 <- matrix(rnorm(s1*r310),nrow = s1)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r320),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r330),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))

    T1 <- matrix(rnorm(s2*r410),nrow = s2)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r420),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(G2*r430),nrow = G2)
    U3 <- as.matrix(qr.Q(qr(T1)))
    T1 <- matrix(rnorm(q*r440),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r410*r420*r430*r440,5,10),nrow = r430)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
  }
  
  
  set.seed(seed_id)
  X1 <- matrix(runif(n*p1), nrow = n)
  X2 <- matrix(runif(n*p2), nrow = n)

  
  f01 = matrix(0,n,q)
  for(j in 1:G1){
    X1j <- X1[,group1==gunique1[j]]
    X11 <- X1j[,1:s1]
    f01 = f01 + bsbasefun(X11,K,degr)%*%t(D3)
  }
  len = s2
  f02 = matrix(0,n,q)
  for(j in 1:G2){
    D23 = D44[,((j-1)*len*K+1):(j*len*K)]
    X2j <- X2[,group2==gunique2[j]]
    X21 <- X2j[,1:s2]
    f02 = f02 + bsbasefun(X21,K,degr)%*%t(D23)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f01 + f02 + eps*sigma2
  
  
  # initialize
  K_index = K
  if(is.null(r1_t3_index)) r1_t3_index = 1:min(4,s1)
  if(is.null(r2_t3_index)) r2_t3_index = 1:min(4,K)
  if(is.null(r3_t3_index)) r3_t3_index = 1:min(4,q)
  if(is.null(r1_t4_index)) r1_t4_index = 1:min(4,s2)
  if(is.null(r2_t4_index)) r2_t4_index = 1:min(4,K)
  if(is.null(r3_t4_index)) r3_t4_index = 1:min(4,G2) 
  if(is.null(r4_t4_index)) r4_t4_index = 1:min(4,q)
  r_index = list(r1_t3_index = r1_t3_index, r2_t3_index = r2_t3_index, r3_t3_index = r3_t3_index, 
                 r1_t4_index = r1_t4_index, r2_t4_index = r2_t4_index, r3_t4_index = r3_t4_index, r4_t4_index = r4_t4_index)
  if(is.null(D0_t3) | is.null(D0_t4)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg1*r310),pg1,r310)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r320),K,r320)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r330),q,r330) 
      C <- qr.Q(qr(T1))
      S = matrix(runif(r310*r320*r330),r330,r310*r320)
      SABC = list(S=S,A=A,B=B,C=C)
      D0_t3 = list()
      for(j in 1:G1) D0_t3[[j]] = SABC
      
      T1 = matrix(rnorm(pg2*r410),pg2,r410)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r420),K,r420)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(G2*r430),G2,r430)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r440),q,r440)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r410*r420*r430*r440),r440,r410*r420*r430)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_t3_max = max(r1_t3_index) 
      r2_t3_max = max(r2_t3_index) 
      r3_t3_max = max(r3_t3_index) 
      K_max = max(K_index)
      T1 = matrix(rnorm(pg1*r1_t3_max),pg1,r1_t3_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_t3_max),K_max,r2_t3_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_t3_max),q,r3_t3_max) #note! r3_t3_max = r4_t4_max
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_t3_max*r2_t3_max*r3_t3_max),r3_t3_max,r1_t3_max*r2_t3_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0_t3 = list()
      for(j in 1:G1) D0_t3[[j]] = SABC
      
      r1_t4_max = max(r1_t4_index) 
      r2_t4_max = max(r2_t4_index) 
      r3_t4_max = max(r3_t4_index) 
      r4_t4_max = max(r4_t4_index) 
      K_max = max(K_index)
      T1 = matrix(rnorm(pg2*r1_t4_max),pg2,r1_t4_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_t4_max),K_max,r2_t4_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(G2*r3_t4_max),G2,r3_t4_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_t4_max),q,r4_t4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_t4_max*r2_t4_max*r3_t4_max*r4_t4_max),r4_t4_max,r1_t4_max*r2_t4_max*r3_t4_max)
      D0_t4 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X1=X1,X2=X2,f01=f01,f02=f02,group=group,D3=D3,D44=D44,n=n,q=q,p=p,s1=s1,s2=s2,r310=r310,r320=r320,r330=r330,r410=r410,r420=r420,r430=r430,r440=r440,
              K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,is.fabs=1,
              r1_t3_index = r1_t3_index, r2_t3_index = r2_t3_index, r3_t3_index = r3_t3_index, 
              r1_t4_index = r1_t4_index, r2_t4_index = r2_t4_index, r3_t4_index = r3_t4_index, r4_t4_index = r4_t4_index,
              r_index=r_index,D0_t3=D0_t3,D0_t4=D0_t4))
}

