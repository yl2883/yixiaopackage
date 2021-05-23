#' Algorithm leveraging function
#'
#' A function that implements algorithmic leveraging for linear regression using uniform and
#' leverage score based subsampling of rows.
#'
#'@param y The response
#'@param X The covariate
#'@param size subsample size of leveraging
#'@param draws Number of draws/repetition of leveraging. Default value is 1.
#'@param method The method to implement, can be weighted leveraging and uniform
#'leveraging. Default is "uniform".
#'@return the leveraged coefficient that fits the linear regression
#'@author Yixiao Lin
#'@export
algo_leverage <- function(y,X,size,draws=1,method ="uniform"){

  ### Check

  X= as.matrix(X)
  n = dim(X)[1]

  if (length(y) != n){
    warning("variable lengths differ")
  }

  if (size%%1 != 0 || size<=0 || is.na(size) ){
    warning("size is invalid")
  }

  if (method == "uniform"|| method  == "weighted"){
    if (method == "uniform"){
      Pi = rep(1/n, n)
    }
    else {
      XTXinv = solve(t(X)%*%X)
      H =X%*%XTXinv%*%t(X)
      H_ii = diag(H)/sum(diag(H))
      Pi = H_ii
    }
    return(lev_subsample(y, X, size, Pi,draws))
  }

  else {
    warning("method should be either uniform or weighted")
  }

}

#' subsample data with draws repetition, size r, probability Pi
#'@param y The response
#'@param X The covariate
#'@param r subsample size of leveraging
#'@param Pi the probability vector to select r sample out of n
#'@param draws Number of draws/repetition of leveraging
#'@return the leveraged coefficient that fits the linear regression
#'@author Yixiao Lin
#'@export
lev_subsample= function(y,X,r,Pi,draws) {
  n=dim(X)[1]
  p=dim(X)[2]
  # error = rep(NA,draws)
  beta=matrix(NA,draws,p)
  for (i in 1:draws){
    index=sample(1:n, size = r, replace = T,prob=Pi)
    Xstar=X[index,]
    Ystar=Y[index]
    Phi=1/ (Pi[index]^{1/2})
    model <- lm(Ystar ~ 0 + Xstar, weights=Phi)
    beta.hat=as.numeric(model$coefficients)
    beta[i,]=beta.hat
    # error[i]=norm(beta.hat-bhat,'2')  #MSE
  }
  return(beta)
}
