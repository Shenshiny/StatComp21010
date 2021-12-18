## -----------------------------------------------------------------------------
RIG <- function(n,mu,lambda)
{
  #As in the case of normal distribution, check lengths
  if(length(n)>1) n<-length(n)
  
  #Check that mu and lambda are positive
  if(any(mu<=0)) stop("mu must be positive")
  if(any(lambda<0)) stop("lambda must be positive")
  
  #Check lengths and adjust them recycling
  if(length(mu)>1 && length(mu)!=n) mu<- rep(mu,length=n)
  if(length(lambda)>1 && length(lambda)!=n) lambda = rep(lambda,length=n)
  
  #Generate random sample from standard normal
  g<-rnorm(n,mean=0,sd=1)
  
  #Transform to  a sample from chi-squared with 1 df
  v<-g*g
  
  w<-mu*v
  cte<-mu/(2*lambda)
  result1<-mu+cte*(w-sqrt(w*(4*lambda+w)))
  result2<-mu*mu/result1
  
  #Uniform random numbers (0,1)
  u<-runif(n)
  
  ifelse(u<mu/(mu+result1),result1,result2)
  
}

## -----------------------------------------------------------------------------
n <- 100
mu <- c(1,1)
lambda <- c(2,4)
RIG(n,mu,lambda)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
library(lars)
library(lasso2)
data("Prostate")
x<- Prostate[,-9]
y<-Prostate[,9]
pmle(x,y,lam0=3)$coef

## -----------------------------------------------------------------------------
pmle(x,y,lam0=3)$sigma

## -----------------------------------------------------------------------------
EMg<-function(e1,X,t_max=100){
  f<-function(e,X){
    sigma1<-as.matrix(e$sigma[[1]])
    sigma2<-as.matrix(e$sigma[[2]])
    M<-matrix(0,dim(X)[1],dim(X)[2])
    for (i in 1:dim(X)[1]){
      X1<-X[i,]
      m1<-dmvnorm(X1,e$mu[[1]],sigma1)
      m2<-dmvnorm(X1,e$mu[[2]],sigma2)
      M[i,1]<-e$lambda[1]*m1
      M[i,2]<-e$lambda[2]*m2
    }
    return(M)
  }
  n<-dim(X)[1]
  for (i in 1:t_max){
    p<-f(e1,X)
    p1<-p/apply(p,1,sum)
    e1$lambda<-c(sum(p1[,1])/n,sum(p1[,2])/n)
    e1$mu<-list(as.vector(p1[,1]%*%X/sum(p1[,1])),as.vector(p1[,2]%*%X/sum(p1[,2])))
    q1<-apply(X,1,function(X) X-e1$mu[[1]])
    sigma1<-q1%*%diag(p1[,1])%*%t(q1)/sum(p1[,1])
    q2<-apply(X,1,function(X) X-e1$mu[[2]])
    sigma2<-q2%*%diag(p1[,2])%*%t(q2)/sum(p1[,2])
    e1$sigma<-list(sigma1,sigma2)
  }
  return(e1)
}

## -----------------------------------------------------------------------------
library(mvtnorm)
e1<-list()
e1$mu<-list(c(5,62),c(6,85))
e1$lambda<-c(.3,.7)
e1$sigma<-list(diag(2),diag(2))
e_hat<-EMg(e1,as.matrix(faithful),100)
e_hat


