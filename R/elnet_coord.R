#' Elastic Net
#'
#' A function that fits elastic net to the given data using coordinate descent algorithm.
#'
#'@param y The response
#'@param X The covariate
#'@param beta The initialization of the coeffient. The default value is 0.
#'@param alpha Elastic net parameter, 0<=alpha<=1
#'@param lambda Regularization parameter
#'@param maxit Maximum Number of Iteration. The default value is 10000.
#'@param tol stop criterion. If the difference between the current result
#' and result from the last iteration is less than the Tolerance, then we say the algorithm converges
#' The default value is 1e-6
#'@return The solution that fits elastic net to the given data using
#'coordinate descent algorithm.
#'
#'@author Yixiao Lin
#'
#'@export
elnet_coord <- function(y,X, beta, alpha,lambda, maxit= 10000, tol = 1e-6){
  X = as.matrix(X)
  n = dim(X)[1]
  p = dim(X)[2]

  ymean = mean(Y); ysd = sd(Y);  Y = (Y-ymean)/ysd;
  Xmean = apply(X,2,mean);  Xsd = apply(X,2,sd);
  X = (X - matrix(Xmean,nrow(X),ncol(X),byrow=TRUE))%*%diag(1/Xsd)

  if (length(y) != n){
    warning("variable lengths differ")
  }

  if (alpha<0 || alpha>1 || is.na(alpha)){
    warning("alpha is invalid")
  }

  if (maxit%%1 != 0 || maxit<0 || is.na(maxit) ){
    warning("niter is invalid")
  }

  if (lambda<0|| is.na(lambda)){
    warning("lambda is invalid")
  }

  if (tol<0 || is.na(tol) ){
    warning("tol is invalid")
  }

  if(any(is.na(beta))){
    beta = rep(0, p)
  }

  res=CD(X,Y,lambda,alpha,beta,n,p,maxit,tol)
  return(res)
}


#' coordinate descent elastic net
#'@param X The covariate
#'@param Y The response
#'@param beta The initialization of the coeffient.
#'@param alpha Elastic net parameter, 0<=alpha<=1
#'@param lambda Regularization parameter
#'@param n numbers of rows of X
#'@param p numbers of rows of Y
#'@param maxit Maximum Number of Iteration.
#'@param tol stop criterion. If the difference between the current result
#' and result from the last iteration is less than the Tolerance, then we say the algorithm converges
#'@export
CD=function(X,Y,lambda,alpha,beta,n,p,maxit,tol){

  tol.met=F
  beta1=beta
  iter=0

  while(!tol.met){

    beta0=beta1
    iter=iter+1

    for (j in 1:p){
      ej=Y-X[,-j]%*% beta1[-j]
      xj=X[,j]
      z=mean(xj*ej)
      num=sign(z)*max(0,abs(z)-lambda*alpha)
      denom=mean(xj^2)+lambda*(1-alpha)
      beta.j=num/denom
      beta1[j]=beta.j
    }

    if (any(abs(beta1) == Inf))
      stop("The algorithm diverges")

    if (iter == maxit) {
      stop("Max iteration reached")
    }

    if (norm(beta0-beta1,'2')<tol){
      tol.met=T
    }
  }
  return(beta1)
}
