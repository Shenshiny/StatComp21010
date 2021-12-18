#' @title RIG for random number generation from Inverse Gaussian Distribution.
#' @description Generate random number from Inverse Gaussian Distribution for Bayesian linear regression.
#' @param n number of observations. If 'length(n) > 1', the length is taken to be the number required.
#' @param mu mean, it can be vector, each element must be bigger than 0
#' @param lambda shape parameter, it can be a vector each element must be bigger than 0
#' @return Random number from Inverse Gaussian Distribution.
#' @examples
#' \dontrun{
#' n <- 100
#' mu <- c(1,1)
#' lambda <- c(2,4)
#'RIG(n,mu,lambda)
#' }
#' @export
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
