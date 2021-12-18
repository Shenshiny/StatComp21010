#' @title A illustration dataset
#' @name X
#' @description A dataset used for predict matrix.
#' @examples
#' \dontrun{
#' data(X)
#' attach(X)
#' X[,1]
#' }
NULL

#' @title A illustration dataset
#' @name Y
#' @description A dataset used for response vector.
#' @examples
#' \dontrun{
#' data(Y)
#' attach(Y)
#' Y
#' }
NULL

pmEst <-function (x, y, lam0) 
{
  nx = dim(x)[1];px = dim(x)[2]
  if (is.null(lam0)) lam0 = sqrt(2 * log(px)/nx)
  sigmaint = 0.1;sigmanew = 5;flag = 0
  objlasso = lars(x, y, type = "lasso", intercept = FALSE, 
                  normalize = FALSE, use.Gram = FALSE)
  while (abs(sigmaint - sigmanew) > 1e-04 & flag <= 100) {
    flag = flag + 1
    sigmaint = sigmanew
    lam = lam0 * sigmaint
    hy = predict.lars(objlasso, x, s = lam, type = "fit", 
                      mode = "lambda")$fit
    sigmanew = sqrt(mean(y*(y - hy)))
  }
  
  hsigma = sigmanew; hlam = lam;
  hbeta = predict.lars(objlasso, x, s = lam, type = "coefficients", 
                       mode = "lambda")$coef
  y.pmle =predict.lars(objlasso, x, s = lam, type = "fit", 
                       mode = "lambda")$fit
  
  return(list(sigma = hsigma,coef = hbeta,y.pmle = y.pmle))
}

lse<-function (X, y, indexset) 
{
  hbeta = rep(0, dim(X)[2])
  if (length(indexset) > 0) {
    objlm = lm(y ~ X[, indexset] + 0)
    hbeta[indexset] = objlm$coefficients
    hsigma = sqrt(mean(objlm$residuals^2))
    residuals = objlm$residuals
    fitted.values = objlm$fitted.values
  }
  else {
    hsigma = sqrt(mean(y^2))
    residuals = y
    fitted.values = rep(0, dim(X)[2])
  }
  return(list(hsigma = hsigma, coefficients = hbeta, residuals = residuals, 
              fitted.values = fitted.values))
}

#' @title PMLE-penalized maximum likelihood estimate in linear regression model.
#' @description A penalized maximum likelihood estimator using R
#' @param x predictors, an n by p matrix with n > 1 and p > 1
#' @param y	 response, an n-vector with n > 1
#' @param lam0 Initial penalty level. Default is sqrt(2 * log(p)/n)
#' @param LSE	 If TRUE, compute least squares estimates after scaled Lasso selection. Default is FALSE
#' @return \item{sigma}{the estimated noise level}
#' \item{coef}{the estimated coefficients}
#' \item{fitted.value}{the fitted value}
#' \item{residuals}{the residuals}
#' \item{lse}{the least square estimation after the seletion}
#' @importFrom lars lars predict.lars
#' @importFrom stats lm
#' @useDynLib StatComp21010
#' @examples 
#' \dontrun{
#' data(X)
#' attach(X)
#' data(Y)
#' attach(Y)
#' x = X
#' y = Y
#' pmle(x,y,lam0=3)$coef
#' }
#' @export
pmle<-function (x, y, lam0 = NULL, LSE = FALSE) 
{ 
  x <- as.matrix(x)
  y <- as.numeric(y)
  est <- pmEst(x, y, lam0=NULL)
  est$fitted.value <- as.vector(x %*% est$coef)
  est$residuals <- y - est$fitted.values
  if (LSE == TRUE) {
    lse.pmle = lse(x,y,indexset = which(est$coef!= 0))
    est$lse = lse.pmle
  }
  est$call <- match.call()
  class(est) <- "PMLE"
  est
}