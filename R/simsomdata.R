#' Simulate some datasets
#'
#' This function simulates some survival data with a cure proportion.
#' @param i Seed value, integer (used for reproducible simulations and multi-dataset simulation).
#' @param n Sample size n, integer.
#' @param ng Number of higher-level units (groups), integer.
#' @param coeff A list of size 2, containing a vector of parameter values each, respectively of the indicence and the latency model.
#' @param REvcov A list of size 3, containing the covariance-matrix elements (varcure, frailty and correlation parameter).
#' @export
#' @keywords simulate
#

simsomdata<-function(i=1,n,ng,coeffic,REvcov){
  set.seed(i)
  k = 1
  pcure = 0.20
  pcens = 0.25

  coeffic[[1]]->gamma0
  coeffic[[2]]->beta0
  REvcov[[1]]->varcure
  REvcov[[2]]->varfrail
  REvcov[[3]]->corr

  randeff<-MASS::mvrnorm(ng,c(0,0),matrix(c(varcure,rep(corr*sqrt(varcure*varfrail),2),varfrail),2),empirical=T)
  grp<-sample(1:ng,n,replace=TRUE)

  Xs = rnorm(n)
  Xd = rnorm(n)
  Ds = sample(c(0,1),n,replace=TRUE)
  Dd = sample(c(0,1),n,replace=TRUE)
  Xmat<-cbind(1,Xs,Xd,Ds,Dd)

  phi = exp(Xmat%*%gamma0+randeff[grp,1])
  phi = phi/(1+phi)
  B = (runif(n) <= phi)

  lambdaX = exp(Xmat%*%beta0+randeff[grp,2])
  uu = runif(n)
  T = 100000*(B==0) + (-log(uu)/lambdaX)^(1/k)*(B==1)

  bC = lambdaX/(1-pcens)*(pcens-pcure)
  uu2 = runif(n)
  C =(-log(uu2)/bC)^(1/k)

  Y = apply(cbind(T,C),1,min)
  Delta = (T <= C)
  data = data.frame(cbind(Y=Y,Delta=Delta,Xmat[,-1],grp=grp))
  names(data)<-c('Y','Delta','Xs','Xd','Ds','Dd','grp')

  return(data)
}
