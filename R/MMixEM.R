#' Expectation-Maximization Algorithm (internal function)
#'
#' (Internal function) Implementation and Parameter optimization using the EM-algorithm.
#' @param Time A vector of length n, containing the observed event/censoring time to event of interest.
#' @param Status A vector of length n, containing the event/censoring indicators (1=event, and 0=censoring).
#' @param X A nxk matrix of covariates for the latency model.
#' @param Z A nxk matrix of covariates for the incidence model.
#' @param coefs A stacked vector of length 2k+2j+1, corresponding to gamma, beta, random and frailty, as provided by initialisation stage.
#' @param varmat A 2x2 variance-covariance matrix of the bivariate gaussian random effects distribution.
#' @param grp A vector of length n, containing the grouping indicators (to identify the frailty term for each subject).
#' @param emmax specifieks the maximum number of iterations for the Expectation-Maximization algorithm. If the convergence criterion is not met, the EM iterations will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default emmas = 50.
#' @param eps sets the convergence criterion (tolerance). If the summed squared changes in the parameters between iterations is lower than the specified value, the algorithm is considered to be converged.
#' @param idx Auxiliary parameter containing the identification indices of the sets of parameters.
#' @export
#' @keywords baseline survival
#

REem <- function(Time,Status,X,Z,coefs,varmat,grp,emmax,eps,idx){

  w<-Status

  idxgamma=idx$idxgamma ;idxbeta=idx$idxbeta; idxrandom=idx$idxrandom;idxfrail=idx$idxfrail ; ncl=idx$ncl

  gamma<-coefs[idxgamma]
  beta<-coefs[idxbeta]
  random<-coefs[idxrandom]
  frailty<-coefs[idxfrail]

  s <- REsmsurv(Time,Status,X,beta,frailty,grp,w)$survival  # initial baseline survival estimator (Breslow-type)

  convergence<- 1000 ; i <-1
  Tstartest<-Sys.time()

  while (convergence > eps & i < emmax){

    ## E step

      uncureprob <- matrix(exp(Z%*%gamma+random[grp])/(1+exp(Z%*%gamma+random[grp])),ncol=1) # uncured probability
      survival<-drop(s^(exp(X%*%beta+frailty[grp]))) # conditional survival probability

    w <- Status+(1-Status)*(uncureprob*survival)/((1-uncureprob)+uncureprob*survival) # CONDITIONAL EXPECTATION OF (UN)CURE STATUS

    ## M step

      # BLUP-type log-likelihood
      minusloglikelihood<-function(coeff){
        gamma<-coeff[idxgamma] ; beta<-coeff[idxbeta]
        random<-coeff[idxrandom] ; frailty<-coeff[idxfrail]

        lpc<-Z%*%gamma+random[grp]
        l1<-as.numeric(logLik(glm(w~0+offset(lpc),family = 'binomial')))

        lps<-X%*%beta+frailty[grp]
        l2<-as.numeric(survival::coxph(Surv(Time,Status)~1+offset(log(w))+offset(lps),subset=w!=0,ties="breslow")$loglik)

        l3<--1/(2*(1-varmat[3]^2))*(sum(random^2)/varmat[1]+sum(frailty^2)/varmat[2]-2*varmat[3]*sum(random*frailty)/sqrt(varmat[1]*varmat[2])) +
          -ncl/2*log(varmat[1]*varmat[2]*(1-varmat[3]^2))

        ll<--(l1+l2+l3)
        return(ll)
      }

      # internal looping for VarComp-estimation
      # try<-cbind(seq(-0.99,0.99,0.01),0)
      # for (i in 1:dim(try)[1]){  varmat[3]=try[i,1] ; try[i,2]=minusloglikelihood(coefs)  }
      # varmat[3]<-try[which.min(try[,2]),1]
        #optpar<-nlm(minusloglikelihood,coefs,hessian=F)$estimate
      optpar<-nloptr::bobyqa(coefs,minusloglikelihood,lower=rep(-4*max(1,varmat),length(coefs)),upper=rep(4*max(1,varmat),length(coefs)),
                     control=list(maxeval=4000))$par


      converge.2<- 1000 ; i.2 <-1
      while(converge.2>eps & i.2<10){

        # BLUP-type log-likelihood optimization
          blup<-nlm(minusloglikelihood,optpar,hessian=T)
          newpar<-blup$estimate

        # REML estimation for variance components

          M1<-diag(ncl) ; M0<-matrix(0,ncol=ncl,nrow=ncl)
          J1<-rbind(cbind(M1,M0),cbind(M0,M0)) ; J2<-rbind(cbind(M0,M1),cbind(M1,M0)) ; J3<-rbind(cbind(M0,M0),cbind(M0,M1))

          T33<-solve(blup$hessian)[c(idxrandom,idxfrail),c(idxrandom,idxfrail)]
          a<-newpar[c(idxrandom,idxfrail)] ; A<-a%*%t(a)

          L1<-sum(diag(J1%*%(T33+A))) ; L2<-sum(diag(J2%*%(T33+A)))/2 ; L3<-sum(diag(J3%*%(T33+A)))

          update_varmat<-c(L1/ncl,L3/ncl,L2/sqrt(L1*L3))

          converge.2<-sum((newpar-optpar)^2)+sum((varmat-update_varmat)^2)
          varmat<-update_varmat
          optpar<-newpar
          i.2<-i.2+1
      }

    update_coefs <- optpar
    update_s <-REsmsurv(Time,Status,X,coefs[idxbeta],coefs[idxfrail],grp,w)$survival

    convergence<-sum(c(coefs-update_coefs)^2)+sum((s-update_s)^2)

    coefs <- update_coefs
    asymvar<-diag(solve(blup$hessian)[c(idxgamma,idxbeta),c(idxgamma,idxbeta)])
    s<-update_s
    i <- i+1
  }
  Tstopest<-Sys.time() ; run.time<-Tstopest-Tstartest

  gamma<-coefs[idxgamma]
  segamma<-sqrt(asymvar[idxgamma])
  beta<-coefs[idxbeta]
  sebeta<-sqrt(asymvar[idxbeta])
  random<-coefs[idxrandom]
  frailty<-coefs[idxfrail]
  uncureprob <- matrix(exp(Z%*%gamma+random[grp])/(1+exp(Z%*%gamma+random[grp])),ncol=1)

  em <- list(gamma=gamma,segamma=segamma, beta= beta,sebeta=sebeta,varmat=varmat,random=random,frailty=frailty,Survival=s,uncureprob=uncureprob,runtime=run.time)
  return(em)
}

