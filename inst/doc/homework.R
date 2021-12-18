## -----------------------------------------------------------------------------
set.seed(1014)
x<-rpois(108,lambda=6)
y<-matrix(rpois(12*9,lambda=6),nrow=12,ncol=9,byrow=TRUE)
x;y

## -----------------------------------------------------------------------------
ts(matrix(rpois(108,6),12,9),frequency=12,start=c(1958,1))

## -----------------------------------------------------------------------------
x<-1:10
y<-runif(10)
z<-rchisq(10,14)
exp1<-expression(x/(y^2+exp(z)))
exp1
eval(exp1)
D(exp1,"x");D(exp1,"y");D(exp1,"z")

## -----------------------------------------------------------------------------
m<-matrix(1:16,4,4)
layout(m,widths=c(1,2,3,4),heights=c(4,3,2,1))
layout.show(16)

## -----------------------------------------------------------------------------
x=seq(-10,10,0.01)
plot(x,exp(((-1/2)*x^2))/sqrt(2*pi),xlim=c(-15,15), ylim=c(0,1), main="标准正态图 ",  
xlab="x", ylab="y")

## -----------------------------------------------------------------------------
set.seed(1014)
n<-1e5
u<-runif(n)
y<--2*log(u)
x1<-sqrt(y);x2<-3*sqrt(y);x3<-5*sqrt(y)  
##Generate Rayleigh(sigma) samples for several choices of sigma>0
hist(x1,prob=TRUE,main = expression(f(x1)==x1*exp(-x1^2/2)))
m1<-seq(0,100,0.1)
lines(m1,m1*exp(-m1^2/2),col="green")
hist(x2,prob=TRUE,main = expression(f(x2)==(x2/9)*exp(-x2^2/18)))
m2<-seq(0,100,0.1)
lines(m2,(m2/9)*exp(-m2^2/18),col="purple")
hist(x3,prob=TRUE,main = expression(f(x3)==(x3/25)*exp(-x3^2/50)))
m3<-seq(0,100,0.1)
lines(m3,(m3/25)*exp(-m3^2/50),col="darkred")
## check that the mode of the generated samples is close to the theoretical mode sigma(check the histogram)

## -----------------------------------------------------------------------------
set.seed(1014)
n<-1000
X1<-rnorm(n,mean=0,sd=1)
X2<-rnorm(n,mean=3,sd=1)
## Generate a random sample of size 1000 from a normal location mixture.
Z1<-0.75*X1+(1-0.75)*X2
hist(Z1,prob=TRUE,main = expression(f(Z1)==0.75*X1+0.25*X2))
x<-seq(-5,5,length.out = 1000)
lines(x,dnorm(x,0.75,0.625),col="darkblue",lty=2)
## Graph the histogram of the sample with density superimposed,for p1=0.75.
p1<-sample(c(0,1),n,replace=TRUE)
Z2<-p1*X1+(1-p1)*X2
hist(Z2)
## The components of the mixture have N(0,1) and N(3,1) distributions with mixing probabilities p1 and p2=1-p1.
Z3<-0.5*X1+0.5*X2;Z4<-0.3*X1+0.7*X2;Z5<-0.1*X1+0.9*X2
hist(Z3);hist(Z4);hist(Z5)
## Repeat with different values for p1 and observe whether the empirical distribution of the mixture appears to be bimodal.

## -----------------------------------------------------------------------------
set.seed(1014)
f<-function(n,t,lambda){
  Nt<-rpois(n,lambda = lambda*t)
  return(Nt)}
Xt<-numeric()
g<-function(alpha, beta){
for(i in 1:length(N)){
Xt[i]<-sum(rgamma(N[i],shape = alpha, scale= 1/beta))}
return(Xt)}
lambda<-seq(1,10,by=1)
alpha<-seq(2,length.out = 10,by=2)
beta<-seq(3,30,by=3)
D<-data.frame()
E_X<-numeric()
Var_X<-numeric()
lambda.t.E_Y<-numeric()
Theoretical_Var<-numeric()
for(i in 1:10){
N<-f(1e4,10,lambda[i])
X<-g(alpha[i],beta[i])
E_X[i]<-mean(X)
Var_X[i]<-(1000-1)*var(X)/1000  
lambda.t.E_Y[i]<-lambda[i]*10*alpha[i]/beta[i] ## theoretical mean :E[X(t)]=lambda*t*E[Yi]
Theoretical_Var[i]<-(lambda[i]*10*alpha[i]*(alpha[i]+1)/(beta[i])^2)}
## theoretical variance  :Var[X(t)]=lambda*t*E[(Ti)^2]
D<-data.frame(lambda,alpha,beta,E_X,Var_X,lambda.t.E_Y,Theoretical_Var)
## putting differeent values of lambda ,alpha ,beta in f() and g() we can simulate compound Poisson process(lambda)-Gauss process(Y has a Gamma distribution).
D

## -----------------------------------------------------------------------------
set.seed(1014)
n<-1e6
f<-function(x)
{
s=0
for(i in 1:n)
{
y=rbeta(1,3,3)
if(y<=x)
s=s+1  
}
n_est<s/n
return(c(n_est,pbeta(x,3,3)))
}


## -----------------------------------------------------------------------------
set.seed(1014)
x<-seq(0.1,0.9,length=9)
m<-1e6
u<-runif(m)
cdf<-numeric(length(x))
for(i in 1:length(x))
 { g<-30*(x[i]^5)*(u^4)-60*(x[i]^4)*(u^3)+30*(x[i]^3)*(u^2)
  cdf[i]<-mean(g)
}
p.value<-pbeta(x,3,3)
print(round(rbind(x,cdf,p.value),5))

## -----------------------------------------------------------------------------
MC.Phi<-function(x,R=1e6,antithetic=TRUE)
{
  u<-runif(R/2)
  sigma<-1
  if(!antithetic)
    v<-runif(R/2)
  else
    v<-1-u
  u<-c(u,v)
  cdf<-numeric(length(x))
  for(i in 1:length(x))
{
  g<-{(u*x[i]^2)/(sigma^2)}*exp(-(u*x[i])^2/(2*sigma^2))
  cdf[i]<-mean(g)
  }
  cdf
}
x<-seq(0.1,3,length=10)
Phi<-pnorm(x)
set.seed(123)
MC1<-MC.Phi(x,antithetic=FALSE)
set.seed(123)
MC2<-MC.Phi(x,antithetic=TRUE)
print(round(rbind(x,MC1,MC2,Phi),5))

## -----------------------------------------------------------------------------
set.seed(1014)
m<-10000
MC1<-MC2<-numeric(m)
x<-4.5
for(i in 1:m)
{
  MC1[i]<-MC.Phi(x,R=1000,antithetic=FALSE)
  MC2[i]<-MC.Phi(x,R=1000,antithetic=TRUE)
}

print(sd(MC1));print(sd(MC2));print((var(MC1)-var(MC2))/var(MC1))*100

## -----------------------------------------------------------------------------
g<-function(x)
{
  (x^2*exp(-x^2/2))/sqrt(2*pi)*(x>1)
}
integrate(g,1,Inf)

## -----------------------------------------------------------------------------
set.seed(1014)
f1<-function(x)
{
  exp(-(x-1))*(x>1)
}
f2<-function(x)
{
  ((1+x^2)^(-1)*(x>1)*4/pi)
}
m<-1e7
u<-runif(m)
x1<-1-log(1-u)
x2<-tan(pi*(1+u)/4)
fg<-cbind(g(x1)/f1(x1),g(x2)/f2(x2))
theta.hat<-se<-numeric(2)
theta.hat<-c(mean(fg[,1]),mean(fg[,2]))
se<-c(sd(fg[,1]),sd(fg[,2]))
rbind(theta.hat,se)

## -----------------------------------------------------------------------------
set.seed(1014)
m<-1e6
theta.hat<-se<-numeric(1)
g<-function(x)
{
  (x^2*exp(-x^2/2))/sqrt(2*pi)*(x>1)
}
integrate(g,1,Inf)
x<-rexp(m,1)
fg<-g(x)/exp(-x)
theta.hat<-mean(fg)
se<-sd(fg)
rbind(theta.hat,se)

## -----------------------------------------------------------------------------
set.seed(1212)
alpha<-0.05
n<-20
m<-10000
UCL<-numeric(m)
LCL<-numeric(m)
for (i in 1:m) 
{
  x<-rchisq(n,2)
  LCL[i]<-mean(x)-qt(alpha/2,n-1,lower.tail=FALSE)*sd(x)/(sqrt(n))
  UCL[i]<-mean(x)+qt(alpha/2,n-1,lower.tail=FALSE)*sd(x)/(sqrt(n))
}
mean(LCL<2&UCL>2) ## coverage probability of the t-interval
UCL1<-numeric(m)
LCL1<-numeric(m)
for (i in 1:m) 
{
  x<-rchisq(n,2)
  LCL1[i]<-mean(x)-qnorm(alpha/2,lower.tail=FALSE)*sd(x)/(sqrt(n))
  UCL1[i]<-mean(x)+qnorm(alpha/2,lower.tail=FALSE)*sd(x)/(sqrt(n))
}
mean(LCL1<2&UCL1>2) ## coverage probability of the normality interval

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
m <- 10000
set.seed(1212)
xbar <- numeric(m)
es <- numeric(m)
for (i in 1:m) {
   xbar[i] <- mean(rchisq(n, df = 2))
   es[i] <- qt((1-alpha/2),n-1)*sd(rchisq(n, df = 2))/sqrt(n)
}
# The symmetric t-interval to estimate a mean is xbar +(-) t(1-α/2)(n-1)s/sqrt(n)
p <- mean((xbar-es)<2 & 2<(xbar+es))
#Judging  whether the mean of the populations falling within the confidence interval generated by the samples subject to the chi-square distribution, thus obtaining the probability that the confidence interval covers the mean
p

## -----------------------------------------------------------------------------
set.seed(1212)
n<-1000
alpha<-0.05
m<-10000
p<-numeric(m)
for (j in 1:m) {
u1=runif(n)
u2=runif(n)
z=sqrt(-2*log(u1))*cos(2*pi*u2)
x=z^2
ttest <- t.test(x, mu = 1)
p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
print(c(p.hat,alpha))

## -----------------------------------------------------------------------------
set.seed(23569)
n<-1000
alpha<-0.05
m<-10000
p<-numeric(m)
for (j in 1:m) {
u=runif(n)
x=2*u
ttest <- t.test(x, mu = 1)
p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
print(c(p.hat,alpha))

## -----------------------------------------------------------------------------
set.seed(17852)
n<-1000
alpha<-0.05
m<-10000
p<-numeric(m)
for (j in 1:m) {
u=runif(1000)
x=-log(1-u)
ttest <- t.test(x, mu = 1)
p[j] <- ttest$p.value
}
p.hat <- mean(p < alpha)
print(c(p.hat,alpha))

## -----------------------------------------------------------------------------
set.seed(1212)
c.v=function(d){
  return(qchisq(0.95,df=d*(d+1)*(d+2)/6))
}
critical.value=c(c.v(1),c.v(2),c.v(3)) ## compute the critiacal values for each d
s<- function(x) {
  statistic=0
  sigma=cov(x,x)
  n=dim(x)[1]
  d=dim(x)[2]
  matrix=(scale(x))%*%t(scale(x))
  statistic=sum(matrix^3)
  return( statistic/n/6 )
} ## computes the sample skewness
n<-c(10,20,30,50,100,500) ## sample sizes
m<-1e4 ## repeat experiment times
Type.I.error=matrix(nrow=length(n),ncol=3)
for(i in 1:length(n)){
  for(d in 1:3){
    result=numeric(m)
    for(j in 1:m){
      x=array(rnorm(n[i]*d),dim=c(n[i],d))
      result[j]=as.integer(abs(s(x)) >=critical.value[d] )
    }
    Type.I.error[i,d]=mean(result)
  }
}
dimnames(Type.I.error)[[1]]=c("size=10","size=20","size=30","size=50","size=100","size=500")
dimnames(Type.I.error)[[2]]=c("d=1","d=2","d=3")
print(Type.I.error)

## -----------------------------------------------------------------------------
set.seed(1212)
alpha<-0.05
n<-100
m<-5000
eps<-c(seq(0, .2, .01), seq(.2, 1, .05))
N<-length(eps)
power1<-numeric(N)
c.v=function(d){
  return(qchisq(0.95,df=d*(d+1)*(d+2)/6))
}
critical.value=c(c.v(1),c.v(2),c.v(3)) ## compute the critiacal values for each d
s<- function(x) {
  statistic=0
  sigma=cov(x,x)
  n=dim(x)[1]
  matrix=(scale(x))%*%t(scale(x))
  statistic=sum(matrix^3)
  return( statistic/n/6 ) 
} ## computes the sample skewness
for (j in 1:N) { 
 e<-eps[j]
 sktests1<-numeric(m)
 for (i in 1:m) { 
 sig1<-sample(c(1, 10), replace = TRUE,size = n*1, prob = c(1-e, e))
 x=array(rnorm(30*1,0,sig1),dim=c(30,1))
 sktests1[i] <- as.integer(s(x) >= c.v(1))
}
power1[j]<-mean(sktests1)
}
plot(eps, power1, type = "b",
xlab = bquote(eps,splice = FALSE), ylim = c(0,1),main="Power of demension one")
abline(a=NULL,b=NULL,h = 0.05, lty = 3)

## -----------------------------------------------------------------------------
set.seed(1212)
alpha<-0.05
n<-100
m<-5000
eps<-c(seq(0, .2, .01), seq(.2, 1, .05))
N<-length(eps)
power2<-numeric(N)
c.v=function(d){
  return(qchisq(0.95,df=d*(d+1)*(d+2)/6))
}
critical.value=c(c.v(1),c.v(2),c.v(3)) ## compute the critiacal values for each d
s<- function(x) {
  statistic=0
  sigma=cov(x,x)
  n=dim(x)[1]
  matrix=(scale(x))%*%t(scale(x))
  statistic=sum(matrix^3)
  return( statistic/n/6 ) 
} ## computes the sample skewness
for (j in 1:N) { 
 e<-eps[j]
 sktests2<-numeric(m)
 for (i in 1:m) { 
  sig2<-sample(c(1, 10), replace = TRUE,size = n*2, prob = c(1-e, e))
 x=array(rnorm(30*2,0,sig2),dim=c(30,2))
 sktests2[i] <- as.integer(s(x) >= c.v(2))
}
power2[j]<-mean(sktests2)
}
plot(eps, power2, type = "b",
xlab = bquote(eps,splice = FALSE), ylim = c(0,1),main="Power of demension two")
abline(a=NULL,b=NULL,h = 0.05, lty = 3)

## -----------------------------------------------------------------------------
set.seed(1212)
alpha<-0.05
n<-100
m<-5000
eps<-c(seq(0, .2, .01), seq(.2, 1, .05))
N<-length(eps)
power3<-numeric(N)
c.v=function(d){
  return(qchisq(0.95,df=d*(d+1)*(d+2)/6))
}
critical.value=c(c.v(1),c.v(2),c.v(3)) ## compute the critiacal values for each d
s<- function(x) {
  statistic=0
  sigma=cov(x,x)
  n=dim(x)[1]
  matrix=(scale(x))%*%t(scale(x))
  statistic=sum(matrix^3)
  return( statistic/n/6 ) 
} ## computes the sample skewness
for (j in 1:N) { 
 e<-eps[j]
 sktests3<-numeric(m)
 for (i in 1:m) { 
  sig3<-sample(c(1, 10), replace = TRUE,size = n*3, prob = c(1-e, e))
 x=array(rnorm(30*3,0,sig3),dim=c(30,3))
 sktests3[i] <- as.integer(s(x) >= c.v(3))
}
power3[j]<-mean(sktests3)
}
plot(eps, power3, type = "b",
xlab = bquote(eps,splice = FALSE), ylim = c(0,1),main="Power of demension 3")
abline(a=NULL,b=NULL,h = 0.05, lty = 3)

## -----------------------------------------------------------------------------
library(bootstrap)
data(scor)
X<-scor
n<-nrow(X)
p<-ncol(X)
set.seed(1212)
B<-1e4
f <- function(X){
m <-cov(X)
lambda <- eigen(m)$values
theta<- lambda[1]/sum(lambda)
}
sigma <- cov(scor)
s <- eigen(sigma)
v <- s$values[order(-(s$values))]
theta <- v[1]/sum(v) ## measures the proportion of variance explained by the first principal component of the original sample.
boot <- 1:B
for(i in 1:B){
index <- sample(1:n,n,replace=TRUE)
X_boot <- X[index,]
boot[i] <- f(X_boot)
}
theta_hat <- f(X)
bias.bootstrap <- mean(boot) - theta_hat ## Obtain the bootstrap estimates of bias of θ.hat
se.bootstrap <- sqrt(var(boot)) ## Obtain the bootstrap estimates of standard error of θ.hat
print(theta,digits = 8)
print(bias.bootstrap,digits = 8)
print(se.bootstrap,digits = 8)

## -----------------------------------------------------------------------------
library(bootstrap)
set.seed(1212)
sigma <- cov(scor)
s <- eigen(sigma)
v <- s$values[order(-(s$values))]
theta <- v[1]/sum(v) ## measures the proportion of variance explained by the first principal component of the original sample.
n <- 88
theta.hat <- numeric(n)
for (i in 1:n) {
  scor.jack <- scor[-i, ]
  sigma.hat <- cov(scor.jack)
  ss <- eigen(sigma.hat)
  vv <- ss$values[order(-(ss$values))]
  theta.hat[i] <- vv[1]/sum(vv)
}
bias.jackknife <- (n-1)*(mean(theta.hat)-theta) ## Obtain the jackknife estimates of bias of θ.hat
se.jackknife <- sqrt((n-1)*mean((theta.hat-theta)^2)) ## Obtain the jackknife estimates of standard error of θ.hat
print(theta,digits = 8)
print(bias.jackknife,digits = 8)
print(se.jackknife,digits = 8)

## -----------------------------------------------------------------------------
library(boot)
set.seed(827)
data(scor, package = "bootstrap")
theta.hat <- numeric(88)
boot.theta<-function(X,i){
  sigma.hat <- cov(X[i,])
  ss <- eigen(sigma.hat)
  vv <- ss$values[order(-(ss$values))]
  theta.hat[i] <- vv[1]/sum(vv)
}
boot.obj <- boot(scor, statistic = boot.theta,R = 2000)
print(boot.obj)
boot.ci(boot.obj,type=c("perc","bca"))

## -----------------------------------------------------------------------------
library(boot)
set.seed(827)
## computes the sample skewness.
skewness <- function(x,i) {
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3)
  m2 <- mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}
s <- 0 ## normal populations (skewness=0)
n <- 20 ## the size of sample
m <- 1000
normal.norm <- normal.basic <- normal.perc <- matrix(0, m, 2)
for (i in 1:m) {
  data.normal <- rnorm(n, 0, 5)
  normal.ske <- boot(data.normal, statistic = skewness, R=1000)
  normal <- boot.ci(normal.ske, type=c("norm","basic","perc"))
  normal.norm[i,] <- normal$norm[2:3]
  normal.basic[i,] <- normal$basic[4:5]
  normal.perc[i,] <- normal$percent[4:5]
}
## Calculate the coverage probability of a normal distribution
    norm <- mean(normal.norm[,1] <= s & normal.norm[,2] >= s)
    basic <- mean(normal.basic[,1] <= s & normal.basic[,2] >= s)
    perc <- mean(normal.perc[,1] <= s & normal.perc[,2] >= s)
## (i)Calculate the probability of the left side of the normal distribution
   norm.left <- mean(normal.norm[,1] >= s )
    basic.left <- mean(normal.basic[,1] >= s )
    perc.left <- mean(normal.perc[,1] >=s )
## (ii)Calculate the probability of the right side of the normal distribution
   norm.right <- mean(normal.norm[,2] <= s )
    basic.right <- mean(normal.basic[,2] <= s )
    perc.right <- mean(normal.perc[,2] <= s)
Distribution<-"N(0,25)"
Type <- c( "norm","basic","perc")
Left <- c(norm.left, basic.left,perc.left)
Right <- c(norm.right, basic.right, perc.right)
P.coverage <- c(norm,basic,perc)
result <- data.frame(Distribution, Type, Left, Right, P.coverage)
knitr::kable(result) 

## -----------------------------------------------------------------------------
library(boot)
set.seed(827)
## computes the sample skewness.
skewness <- function(x,i) {
  xbar <- mean(x[i])
  m3 <- mean((x[i] - xbar)^3)
  m2 <- mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}
s<-sqrt(8/5) ## chisq with 5 freedom distribution (positive skewness=1.6).
n<-20 ## the size of sample
m<-1000
chisq.norm<-chisq.basic<-chisq.perc<-matrix(0, m, 2)
for (i in 1:m) {
  data.chisq<-rchisq(n,5)
  chisq.ske<-boot(data.chisq,statistic=skewness, R=1000)
  chisq<- boot.ci(chisq.ske,type=c("norm","basic","perc"))
  chisq.norm[i,]<-chisq$norm[2:3];
  chisq.basic[i,]<-chisq$basic[4:5];
  chisq.perc[i,]<-chisq$percent[4:5];
}
## Calculate the coverage probability of the chi-square distribution.
    norm <- mean(chisq.norm[,1] <= s & chisq.norm[,2] >= s)
    basic <- mean(chisq.basic[,1] <= s & chisq.basic[,2] >= s)
    perc <- mean(chisq.perc[,1] <= s & chisq.perc[,2] >= s)
## (i) Calculate the probability of the left side of the chi-square distribution.
    norm.left <- mean(chisq.norm[,1] >= s )
    basic.left <-mean(chisq.basic[,1] >= s )
    perc.left <- mean(chisq.perc[,1] >=s )
## (ii) Calculate the probability of the right side of the chi-square distribution.
    norm.right <- mean(chisq.norm[,2] <= s )
    basic.right <- mean(chisq.basic[,2] <= s )
    perc.right <- mean(chisq.perc[,2] <= s)
Distribution<-"chisq(5)"
Type <- c( "norm","basic","perc")
Left <- c(norm.left, basic.left,perc.left)
Right <- c(norm.right, basic.right, perc.right)
P.coverage <- c(norm,basic,perc)
result <- data.frame(Distribution, Type, Left, Right, P.coverage)
knitr::kable(result) 

## -----------------------------------------------------------------------------
set.seed(1212)
data <- as.matrix(iris[1:50, 1:4])
x <- data[ , 3:4]
y <- data[ , 1:2]
R <- 999
z <- c(x,y)
K <- 1:200
n <- length(x)
permutation <- numeric(R)
cor0 <- cor.test(x,y)$statistic
for (i in 1:R) {
  k <- sample(K,size = n,replace = FALSE)
  x1 <- z[k]
  y1 <- z[-k]
  permutation[i] <- cor.test(x1,y1)$statistic
}
permutation.p <- mean(abs(c(cor0,permutation))>=abs(cor0))
permutation.p
cor.test(x, y,alternative = "two.sided", conf.level = 0.95)

## -----------------------------------------------------------------------------
set.seed(1212)
library(boot)
library(RANN)
library(energy)
library(Ball)
library(ggplot2)
## write the function of the NN test.
Tn <- function(z, ix, sizes,k){
  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[ix, ];
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1]
  block2 <- NN$nn.idx[(n1+1):n,-1]
  i1 <- sum(block1 < n1 + .5)
  i2 <- sum(block2 > n1+.5)
  (i1 + i2) / (k * n)
} ## the function of NN method.
NN_test <- function(x,y,NumR=999,k=3){
N <- c(length(x), length(y)) 
z <- c(x,y)
boot.obj <- boot(data = z, statistic = Tn, R = NumR, sim = "permutation", sizes = N,k=k)
ts <- c(boot.obj$t0,boot.obj$t) 
p.value <- mean(ts>=ts[1])
return(p.value)
}
## Energy test
Energy_test <- function(x,y,NumR=999){
N <- c(length(x), length(y)) 
z <- c(x,y)
boot.obs <- eqdist.etest(z, sizes=N, R=NumR) 
p.value <- boot.obs$p.value
return(p.value)
}
alpha <- 0.05
## Ball test: ball statistic test for equal distributions can be directly achieved by  Ball package.

## -----------------------------------------------------------------------------
set.seed(1212)
m <- 1000; ## number of permutation tests
k<-3; ## boot parameter
p<-2; ## dimension of data
mu <- 0; 
n1 <- n2 <- 50; ## the sample size of x and y
R<-999; ## boot parameter
sigma1 <- 1 
sigma2 <- 2
## x~N(0,1) 50 y~N(0,4) 50
p.values1 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rnorm(n1,mu,sigma1)
y <- rnorm(n2,mu,sigma2)
z <- c(x,y) 
p.values1[i,1] <- NN_test(x,y,R,k) 
p.values1[i,2] <- Energy_test(x,y,R)
p.values1[i,3] <- bd.test(x,y,R=R,seed=i*1212)$p.value 
}
pow <- colMeans(p.values1<alpha)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
power
ggplot(power,aes(methods,pow))+geom_col(fill = 'lightblue')+coord_flip()

## -----------------------------------------------------------------------------
set.seed(1212)
m <- 1000; ## number of permutation tests
k<-3; ## boot parameter
p<-2; ## dimension of data
n1 <- n2 <- 50; ## the sample size of x and y
R<-999; ## boot parameter
mu <- 0.5
sd <- 1.5
#x~N(0,1) 50  y~N(0.5,2.25) 50
p.values2 <- matrix(NA,m,3) 
for(i in 1:m){
  x <- matrix(rnorm(n1*p),ncol=p)
  y <- matrix(rnorm(n2*p,mean=mu,sd=sd),ncol=p)
  z <- rbind(x,y)
p.values2[i,1] <- NN_test(x,y,R,k) 
p.values2[i,2] <- Energy_test(x,y,R)
p.values2[i,3] <- bd.test(x,y,R=R,seed=i*1212)$p.value 
}
pow <- colMeans(p.values2<alpha)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
power
ggplot(power,aes(methods,pow))+geom_col(fill = 'lightblue')+coord_flip()

## -----------------------------------------------------------------------------
set.seed(1212)
m <- 1000; ## number of permutation tests
k<-3; ## boot parameter
n1 <- n2 <- 50; ## the sample size of x and y
df=1 ## freedom of t distribution(heavy-tailed distribtion)
R<-999;## boot parameter
p=0.4
mu1 =-1
mu2 = 1
sigma1 =1
sigma2 =2
#x~t1 50  y~0.3N(-1,1)+0.7N(1,2) 50
p.values3 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rt(n1,df=df)
y <- 0.3 * rnorm(n2,mu1,sigma1) + 0.7*rnorm(n2,mu2,sigma2)
z <- c(x,y) 
p.values3[i,1] <- NN_test(x,y,R,k) 
p.values3[i,2] <- Energy_test(x,y,R)
p.values3[i,3] <- bd.test(x,y,R=R,seed=i*1212)$p.value 
}
pow <- colMeans(p.values3<alpha)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
power
ggplot(power,aes(methods,pow))+geom_col(fill ='lightblue')+coord_flip()

## -----------------------------------------------------------------------------
set.seed(122)
m <- 1000; ## number of permutation tests
k<-3; ## boot parameter
n1 <- 100
n2 <- 10
R<-999; ## boot parameter
mu1 =-1
mu2 = 0
sigma1 =1
sigma2 =2
#x~N(-1,1) 100  y~N(1,2) 10
p.values4 <- matrix(NA,m,3) 
for(i in 1:m){ 
x <- rnorm(n1,mu1,sigma1)
y <- rnorm(n2,mu2,sigma2)
z <- c(x,y)
p.values4[i,1] <- NN_test(x,y,R,k) 
p.values4[i,2] <- Energy_test(x,y,R)
p.values4[i,3] <- bd.test(x,y,R=R,seed=i*1212)$p.value 
}
pow <- colMeans(p.values4<alpha)
power <- data.frame(methods = c('NN','energy','Ball'),pow)
power
ggplot(power,aes(methods,pow))+geom_col(fill ='lightblue')+coord_flip()

## -----------------------------------------------------------------------------
set.seed(2029)
m <- 3000 ## Sample size
x <- numeric(m)
u <- runif(m)
theta=1
eta=0
x[1] <- rnorm(1)
k <- 0

## The following function evaluates the Cauchy(0,1) density.
f <- function(x, theta=1, eta=0){
  out <- theta /(pi * (theta^2+(x-eta)^2))
  return(out)
}
for(i in 2:m){
  xt <- x[i-1]
  y <- rnorm(1,mean=xt) ## At each transition, the candidate point Y is generated from N(mean=xt)
  R <- f(y)*dnorm(xt,mean=y)/(f(xt)*dnorm(y,mean=xt)) ## accept probability
  if(u[i] <= R){
    x[i] <- y
  }else{
    x[i] <- xt
    k <- k+1 ## y is rejected
  }
}
print(k/m)
I <- 1001:m ##  The following code will display a partial plot starting at time index 1000.
plot(I,x[I],type="l")
hist(x[I], probability=TRUE,breaks=100)
plot.x <- seq(min(x[I]),max(x[I]),0.01)
lines(plot.x,f(plot.x))
#compare the deciles
observations <- quantile(x[I],seq(0,1,0.1))
expectations <- qcauchy(seq(0,1,0.1))
decile <- data.frame(observations,expectations)
decile

## -----------------------------------------------------------------------------
qqnorm(x[I])
qqline(x[I])
sequence <- seq(0,1,0.01)
standard.Cauchy <- qcauchy(sequence)
standard.Cauchy <- standard.Cauchy[(standard.Cauchy> -Inf) 
                                 & (standard.Cauchy< Inf)]
hist(x[I], freq = FALSE)
lines(standard.Cauchy, dcauchy(standard.Cauchy), lty = 2)

## -----------------------------------------------------------------------------
set.seed(1212)
f = function(x,y){
choose(n,x) * y^(x + a - 1) * (1 - y)^(n - x + b -1)
} ## x = 0,1,2,...,n and 0 < y < 1
N <- 10000   ## sample size
## fixed a,b and n
a <- 2
b <- 4
n <- 16
X=NULL ## Vector of X values
X[1]=5 ## initial value of X given to be 5
Y=NULL ## Vector of Y values
Y[1]=0.33 ## initial value of Y given to be 0.33
Z=matrix(0,N,2) ## the chain,a bivariate sample
## generate the chain
for( i in 2:N){
   X[i]=rbinom(1,n,Y[i-1])
   Y[i]=rbeta(1,(X[i]+a),(n-X[i]+b))
   Z[i,]=c(X[i],Y[i])
}
sample_X <- Z[(1:N),1]## The first column of Z gives X values of the sample
sample_Y <- Z[(1:N),2]## The second column of Z gives Y values of the sample

sample<-round(c(mean(sample_X),mean(sample_Y)),4) #Estimated mean of marginal random variable X from sample
true<-round(c(n*Y[1],(X[1]+a)/(X[1]+a+n-X[1]+b)),4)
compare <- data.frame(sample,true)
compare

## -----------------------------------------------------------------------------
# Gelman-Rubin statistic
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)

  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

## generates a Metropolis chain for Cauchy(0,1)
## with Normal(mean=xt) proposal distribution
## and starting value X1
Standard_Cauchy_Chain <- function(N, X1){
  X <- numeric(N)
  X[1] <- X1    #初始值
  for(i in 2:N){
    Xt <- X[i-1]
    Y <- rnorm(1,0,abs(Xt))
    r <- dt(Y,1)*dnorm(Xt,0,abs(Y))/dt(Xt,1)/dnorm(Y,0,abs(Xt))
    U <- runif(1)
    if(r > 1) r <- 1
    if(U <= r) X[i] <- Y
    else X[i] <- Xt
  }
  return(X)
}
k <- 4      
N <- 8000
b <- 1000     #burn-in length
X1 <- c(0.1,0.2,0.1,0.2)    #initial value

set.seed(12345)
X <- matrix(0, nrow = k, ncol = N)
for(i in 1:k){
  X[i,] <- Standard_Cauchy_Chain(N, X1[i])
}

# compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))
print(Gelman.Rubin(psi))
for (i in 1:k)
  if(i==1){
    plot((b+1):N,psi[i, (b+1):N],ylim=c(-1,1), type="l",
         xlab='Index', ylab=bquote(phi))
  }else{
      lines(psi[i, (b+1):N], col=i)
  }
#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
# Gelman-Rubin statistic
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)

  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}
Bivariate.Gibbs <- function(N, X1){
  a <- b <- 1
  X <- matrix(0, N, 2)
  X[1,] <- X1    #初始值
  for(i in 2:N){
    X2 <-  X[i-1, 2]
    X[i,1] <- rbinom(1,25,X2)
    X1 <- X[i,1]
    X[i,2] <- rbeta(1,X1+a,25-X1+b)
  }
  return(X)
}
k <- 4          
N <- 8000 
b <- 1000    #burn-in length
X1 <- cbind(c(2,7,10,15),runif(4)) # initial value

set.seed(12345)
X <- matrix(0, nrow=k, ncol=N)
Y <- matrix(0, nrow=k, ncol=N)
for (i in 1:k){
  BG <- Bivariate.Gibbs(N, X1[i,])
  X[i, ] <- BG[,1]
  Y[i, ] <- BG[,2]
}
# Consider X

#compute diagnostic statistics
psi <- t(apply(X, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

# Consider Y

#compute diagnostic statistics
psi <- t(apply(Y, 1, cumsum))
for (i in 1:nrow(psi))
  psi[i,] <- psi[i,] / (1:ncol(psi))

#plot the sequence of R-hat statistics
rhat <- rep(0, N)
for (j in (b+1):N)
  rhat[j] <- Gelman.Rubin(psi[,1:j])
plot(rhat[(b+1):N], type="l", xlab="", ylab="R")
abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
## (a)Write a function to compute the k^{th} term
f=function(a,k){
d=length(a)
((-1)^k/(factorial(k)*2^k))*(sqrt(sum(a*a))^(2*k+2)/((2*k+1)*(2*k+2)))*((gamma((d+1)/2)*gamma(k+1.5))/gamma(k+d/2+1))
}
## (b)Modify the function so that it computes and returns the sum
a <- c(1,2) ## value of a
p <- 0
q <- 0
k <- 1
while(is.nan(p)==FALSE){ 
## we stop calculaton when r starts providing NaN for k^{th} term calculation
q[k]=f(a,k-1) ## the k^{th} term of series
p=q[k]
k=k+1 ## increase the k by 1
}
## (c)Evaluate the sum when $a = (1, 2)^T$
sum(q[-length(q)])

## -----------------------------------------------------------------------------
# find the intervals in which the roots fall  
f <- function(a, k) 1-pt(sqrt(a^2*k/(k+1-a^2)), k)
k <- c(4, 25, 100, 1000)
par(mar=c(1,1,1,1))
for (i in k) {
  g <- function(x) f(x, i-1)-f(x, i)
    a <- seq(0, sqrt(i), .1)
  plot(a, g(a), type = 'l', main = paste('k=',i))
  abline(h = 0)
}

## -----------------------------------------------------------------------------
f <- function(a, k) 1-pt(sqrt(a^2*k/(k+1-a^2)), k)
k <- c(4:25, 100, 500, 1000)
Ak <- numeric(length(k))
i <- 1
for (j in k) {
  g <- function(x) f(x, j-1)-f(x, j)
  Ak[i] <- uniroot(g, lower = 1, upper = 2)$root
  i <- i+1
}
knitr::kable(cbind(k,Ak))

## -----------------------------------------------------------------------------
f <- function(k) 2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
ck <- function(a, k) sqrt(a^2*k/(k+1-a^2))
g <- function(u, k) (1+u^2/k)^(-(k+1)/2)
k <- c(4:25, 100, 500, 1000)
root <- numeric(length(k))
i <- 1
for (j in k) {
  ff <- function(a) f(j)*integrate(function(u) {g(u, j)}, 0, ck(a, j))$value-f(j-1)*integrate(function(u) {g(u, j-1)}, 0, ck(a, j-1))$value 
  root[i] <- uniroot(ff, lower = 1, upper = 2)$root
  i <- i+1
}
knitr::kable(cbind(k, Ak, root))

## -----------------------------------------------------------------------------
tau=1
n=10
n_0=7
y=c(0.54,0.48,0.33,0.43,0.91,0.21,0.85)
lambda.hat <- (tau*(n-n_0)+sum(y))/n_0
lambda.hat

## -----------------------------------------------------------------------------
y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
lambda0 <- 0.5

EMlambda <- function(y, lambda){
  N <- 10000
  tol <- .Machine$double.eps^0.5
  n <- length(y)
  nt <- sum(y==1)
  
  for(i in 1:N){
    total <- sum(y)
    lambda_old <- lambda
    lambda <- (total + nt*lambda_old) / n
    
    if(abs(lambda - lambda_old)/lambda_old <tol) break
  }
  
  return(lambda)
}

EMlambda(y, lambda0)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)
lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared

## -----------------------------------------------------------------------------
# 3
data(mtcars)
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1 / disp),
  mpg ~ disp + wt,
  mpg ~ I(1 / disp) + wt
)
models <- lapply(formulas, function(formula) lm(formula, mtcars))

## -----------------------------------------------------------------------------
lapply(models, function(model) rsq(model))

## -----------------------------------------------------------------------------
# 4
data("mtcars")
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})

## -----------------------------------------------------------------------------
# for loop
models_loop <- list()
for(i in 1:10){
  model <- lm(mpg ~ disp, data = bootstraps[[i]])
  models_loop <- c(models_loop, list(model))
}

## -----------------------------------------------------------------------------
# R^2
lapply(models_loop, function(model) rsq(model))

## -----------------------------------------------------------------------------
# lapply
models_lapply <- lapply(bootstraps, function(data) lm(mpg ~ disp, data))

## -----------------------------------------------------------------------------
# R^2
lapply(models_lapply, function(model) rsq(model))

## -----------------------------------------------------------------------------
cars ##  a numeric data.frame.
vapply(cars, sd, numeric(1))

## -----------------------------------------------------------------------------
## a mixed data.frame"mtcars".
vapply(mtcars[vapply(mtcars, is.numeric, logical(1))],sd, numeric(1)) ## use vapply() twice.

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('
NumericMatrix gibbsc(int N,int burn,int a,int b,int n){
  NumericMatrix X(N,2);
  NumericMatrix XX(N-burn,2);
  float x,y;
  X(0,0)=1;
  X(0,1)=0.5;
  for(int i=1;i<N;i++){
    y=X(i-1,1);
    X(i,0)=rbinom(1,n,y)[0];
    x=X(i,0);
    X(i,1)=rbeta(1,x+a,n-x+b)[0];
  }
  for(int k=0;k<N-burn;k++){
    XX(k,0)=X(k+burn,0);
    XX(k,1)=X(k+burn,1);
  }
  return XX;
}
')
res1<-gibbsc(5000,1000,2,2,10)

## -----------------------------------------------------------------------------
gibbsr<-function(N,burn,a=2,b=2,n=10){
  X <- matrix(0, N, 2) #the chain, a bivariate sample
  ###### generate the chain #####
  X[1, ] <- c(1, .5) #initialize
  for (i in 2:N) {
    y <- X[i-1, 2]
    X[i, 1] <- rbinom(1, n, y)
    x <- X[i, 1]
    X[i, 2] <- rbeta(1,x+a,n-x+b)
    }
  b <- burn + 1
  return(X[b:N, ])
}
res2<-gibbsr(5000,1000,2,2,10)

## -----------------------------------------------------------------------------
a <- ppoints(50)
Q1x <- quantile(res1[,1], a)   #quantiles of Cauchy
Q2x <- quantile(res2[,1], a)
Q1y <- quantile(res1[,2], a)   #quantiles of Cauchy
Q2y <- quantile(res2[,2], a)
qqplot(Q1x, Q2x, main="",
        xlab="Rcpp Quantiles of x", ylab="R Quantiles of x")
    abline(c(0,0),c(1,1),col='blue',lwd=2)
qqplot(Q1y, Q2y, main="",
        xlab="Rcpp Quantiles of y", ylab="R Quantiles of y")
    abline(c(0,0),c(1,1),col='blue',lwd=2)

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(Rcpp=gibbsc(5000,1000,2,2,10),R=gibbsr(5000,1000,2,2,10))
summary(ts)[,c(1,3,5,6)]

