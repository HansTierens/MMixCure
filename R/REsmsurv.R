#' Estimate the baseline survival function S_0(t)
#'
#' This function estimates the baseline survival function given the survival data, covariates and frailty predictions.
#' @param Time A vector of length n, containing the observed event/censoring times.
#' @param Status A vector of length n, containing the event/censoring indicators (event:1 and censoring:0).
#' @param X A nxk matrix of covariates..
#' @param beta A vector of length k, containing the covariate effects.
#' @param frailty A vector of length j, containing all j frailty predictions.
#' @param grp A vector of length n, containing the grouping indicators (to identify the frailty term for each subject).
#' @param w A vector of length n, containing the expected susceptibility indicator.
#' @export
#' @keywords baseline survival
#

REsmsurv <- function(Time,Status,X,beta,frailty,grp,w){
  death_point <- sort(unique(subset(Time, Status==1)))

  coxexp <- exp((beta)%*%t(X)+frailty[grp])
  lambda <- numeric()
  event <- numeric()
  for(i in 1: length(death_point)){
    event[i] <- sum(Status*as.numeric(Time==death_point[i]))
    temp <- sum(as.numeric(Time>=death_point[i])*w*drop(coxexp))
    temp1 <- event[i]
    lambda[i] <- temp1/temp
  }
  HHazard <- numeric()
  for(i in 1:length(Time)){
    HHazard[i] <- sum(as.numeric(Time[i]>=death_point)*lambda)
    if(Time[i]>max(death_point)) HHazard[i] <- Inf
    if(Time[i]<min(death_point)) HHazard[i] <- 0
  }
  survival <- exp(-HHazard)
  list(survival=survival)
}
