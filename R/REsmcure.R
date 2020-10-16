#' Mixture Cure Model with Correlated Random Effects
#'
#' Implementation and Parameter optimization using the EM-algorithm.
#' @param formula Double-sided formula ... .
#' @param cureform Right-hand-side formula for the latent incidence model.
#' @param RE String input identifying the variable name corresponding to the grouping variable.
#' @param rho Logical value. TRUE indicates estimation of the correlation parameter. FALSE sets the correlation to zero.
#' @param data data.frame for the dataset which contains all input data.
#' @param emmax Maximum number of iterations for the Expectation-Maximization Algorithm (default = 50).
#' @param eps Hyperparameter: convergence criterion (tolerance) on squared changes between iterations (default =1e-07).
#' @export
#' @keywords baseline survival
#

REsmcure<-function (formula, cureform,RE='grp',rho=F, data, emmax = 50, eps = 1e-07) {
  call<-match.call()

  n <- dim(data)[1]
  mf <- model.frame(formula, data)
  cvars <- all.vars(cureform)
  Z <- as.matrix(cbind(rep(1, n), data[, cvars]))
  colnames(Z) <- c("(Intercept)", cvars)

  Y <- model.extract(mf, "response")
  Xm <- model.matrix(attr(mf, "terms"), mf)
  X <- as.matrix(Xm[,-1])
  if (!inherits(Y, "Surv"))
    stop("Response must be a survival object")
  Time <- Y[, 1]
  Status <- Y[, 2]
  grp <- as.numeric(factor(data[,RE]))

  gammaname <- colnames(Z)
  idxgamma <- c(1:ncol(Z))
  betaname <- colnames(Xm)[-1]
  idxbeta <- c((1+max(idxgamma)):(max(idxgamma)+ncol(X)))
  ncl<-length(unique(grp))
  idxrandom<-c((1+max(idxbeta)):(max(idxbeta)+ncl))
  idxfrail<-c((1+max(idxrandom)):(max(idxrandom)+ncl))
  idx<-list(idxgamma=idxgamma,idxbeta=idxbeta,idxrandom=idxrandom,idxfrail=idxfrail,ncl=ncl)

  # Generate INITIAL parameter values ####
    w <- Status

    cumo<-lme4::glmer(w~0+(1|grp)+Z,family='binomial')
    gamma<-c(lme4::fixef(cumo))
    random<-c(lme4::ranef(cumo)$grp[,1])
    var_random<-as.numeric(lme4::VarCorr(cumo)$grp[1])

    sumo<-coxme::coxme(Surv(Time, Status) ~ (1|grp)+ X + offset(log(w)), subset = w != 0, ties = "breslow")
    beta <- sumo$coef
    frailty <- c(lme4::ranef(sumo)$grp)
    var_frailty<-as.numeric(sumo$vcoef$grp)

  init<-c(gamma,beta,random,frailty)

  cor=0 ; if(rho==T){cor=ifelse(var(random)<=eps|var(frailty)<=eps,eps,
                                ifelse(length(random)!=length(frailty),eps,cor(random,frailty)))}
  varmat<-c(var_random,var_frailty,cor)
  #####

  emfit <- suppressWarnings(REem(Time, Status, X, Z, init, varmat, grp, emmax, eps,idx))


  fit <- list()
  class(fit) <- c("smcure")
  fit$call <- call

  fit$gname <- gammaname
  fit$gamma <- emfit$gamma
  fit$segamma <- emfit$segamma

  fit$bname <- betaname
  fit$beta <- emfit$beta
  fit$sebeta <- emfit$sebeta

    RE<-cbind(emfit$random,emfit$frailty); colnames(RE)<-c('incidence','cure')
  fit$RE <-RE

    REvar<-matrix(c(emfit$varmat[1],emfit$varmat[3]*sqrt(emfit$varmat[1]*emfit$varmat[2]),emfit$varmat[3]*sqrt(emfit$varmat[1]*emfit$varmat[2]),emfit$varmat[2]),ncol=2,byrow=TRUE)
  fit$REVar<-REvar

  fit$Survival<- emfit$Survival
  fit$CureProb <- 1-emfit$uncureprob
  fit$runtime <- emfit$runtime
  fit$varmat<-emfit$varmat

  return(fit)
}
