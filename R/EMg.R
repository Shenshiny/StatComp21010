#' @title EM solution of binary mixed Gaussian model.
#' @description EM solution of binary mixed Gaussian model.
#' @param e1 a list consisting of the iterative initial value.
#' @param X a data matrix or data.frame.
#' @param t_max  maximum number of iterations.
#' @return a list consisting of mu ,sigma matrix and the The mixing ratio.
#' @importFrom mixtools dmvnorm
#' @examples
#' \dontrun{
#'e1<-list()
#'e1$mu<-list(c(5,62),c(6,85))
#'e1$lambda<-c(.3,.7)
#'e1$sigma<-list(diag(2),diag(2))
#'e_hat<-EMg(e1,as.matrix(faithful),100)
#' }
#' @export
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
