#' Simulate (Robust) Standard Errors using a bootstrap procedure or jackknife approach.
#'
#' Calculate (Robust) Standard Errors using either a bootstrap procedure or a Jackknife estimate.
#' @param fit A fitted model of class smcure.
#' @param data The dataset on which the fitted model was fit.
#' @param method Either 'boot' for a bootstrap procedure, or 'jack' for a jackknife procedure.
#' @param k the number of observations to omit from the dataset for the leave-k-out Jackknife estimate on the standard errors. Default = 1, for al leave-one-out jackknife procedure.
#' @param nboot Number of iterations in the sampling-resampling procedure of the bootstrap estimator.
#' @param nmax Maximum number of resampling iterations for the bootstrap and jackknife estimators.
#' @export
#' @keywords robust standard error
#

#########################################
#  RANDOM EFFECTS IN MIXTURE CURE MODEL #
#########################################


sebootjack<-function(fit,data,method=c('boot','jack'),k=1,nboot=100,nmax=1000){
  n=dim(data)[1]
  method=match.arg(method)
  if (method=='boot') {
    nsample <- max(min(nboot,nmax),100)
    k=NULL
    for (i in 1:nsample){
      set.seed(i)
      di=data[sample(1:n,n,rep=T),]
      fiti<-update(fit,data=di)
      if(i==1){
        gs = matrix(fiti$gamma,nrow=1)
        bs = matrix(fiti$beta,nrow=1)
      }
      else{
        gs = rbind(gs,matrix(fiti$gamma,nrow=1))
        bs = rbind(bs,matrix(fiti$beta,nrow=1))
      }
    }
    se.g <- sqrt(nsample/(nsample-1) * colMeans((bs-matrix(rep(colMeans(bs),nsample),nrow=nsample,byrow=T))^2))
    se.b <- sqrt(nsample/(nsample-1) * colMeans((betas-matrix(rep(colMeans(betas),nsample),nrow=nsample,byrow=T))^2))
    info<-noquote(paste('Standard Errors were computed using ',nsample,' Bootstrap Samples',sep=''))
    simse<-list(se.g,se.b,info)
  }

  if(method=='jack'){
    nboot=NULL
    nsample=floor(n/k) ; k=k
    if(floor(n/k)<10){ nsample=100 ; k=floor(n/nsample)}
    if(floor(n/k)>nmax){ nsample=nmax ; k=floor(n/nmax)}

    samplemat = suppressWarnings(matrix(sample(1:n,n,rep=F),ncol=k,nrow=nsample))

    for(i in 1:nsample){
      di=data[-c(samplemat[i,]),]
      fiti<-update(fit,data=di)
      if(i==1){
        gs = matrix(fiti$gamma,nrow=1)
        bs = matrix(fiti$beta,nrow=1)
      }
      else{
        gs = rbind(gs,matrix(fiti$gamma,nrow=1))
        bs = rbind(bs,matrix(fiti$beta,nrow=1))
      }
    }
    se.g <- sqrt((n-k)/k * colMeans((bs-matrix(rep(colMeans(bs),nsample),nrow=nsample,byrow=T))^2))
    se.b <- sqrt((n-k)/k * colMeans((betas-matrix(rep(colMeans(betas),nsample),nrow=nsample,byrow=T))^2))
    info<-noquote(paste('Standard Errors were computed using ',nsample,' ','Jackknife (',k,'-deleted) Samples',sep=''))
    simse<-list(se.g,se.b,info)
  }

  fit=append(fit,simse)
  fit<<-fit
}

