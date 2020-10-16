#' Print REsmcure Object
#'
#' Output of REsmcure object.
#' @param fit an object of REsmcure.
#' @export
#' @keywords output
#'

print.REsmcure<-function(fit){
  x<-fit

  cat('\n Call:\n')
  dput(x$call)

  cat('\nCure probability model:\n')
    cure <- array(round(x$gamma,5), c(length(x$gamma), 4)) ; rownames(cure) <- x$gname ; colnames(cure) <- c('Estimate','Std.Error','Z value','Pr(>|Z|)')
    cure[,2]<-round(x$segamma,5)
    cure[,3]<-round(x$gamma / x$segamma,5)
    cure[,4]<-round((1-pnorm(abs(x$gamma/x$segamma)))*2,5)
    print(cure)

  cat('\n\nFailure time distribution model:\n')
    surv <- array(round(x$beta,5), c(length(x$beta), 4)) ; rownames(surv) <- x$bname ; colnames(surv) <- c('Estimate','Std.Error','Z value','Pr(>|Z|)')
    surv[,2]<-round(x$sebeta,5)
    surv[,3]<-round(x$beta / x$sebeta,5)
    surv[,4]<-round((1-pnorm(abs(x$beta/x$sebeta)))*2,5)
    print(surv)

  cat('\nRandom Effects:\n')
    rn<-max(nchar(round(1/min(x$REVar))),3)
    REmat<-cbind(format(diag(x$REVar),digits=rn),c('',round(x$REVar[1,2]/sqrt(prod(diag(x$REVar))),rn))) ; rownames(REmat)<-c('Cure/Incidence','Survival/Latency') ; colnames(REmat)<-c('Variance','Correlation')
    print(as.data.frame(REmat))
  invisible(x)
}
