## Code written by Peter Craigmile
###################################

########################################################
plot.poly <- function (xx.polylist, aux.var, intervals,
                       legend.x, legend.y, ...) {
  ## ======================================================================
  ## Plotting spatial polygons, with an auxiliary variable 'aux.var',
  ## broken down by the 'intervals'.
  ## ======================================================================
  
  if (missing(aux.var)) {
    plot(xx.polylist, ...)
    ##forcefill=FALSE, ...)
  }
  else {
    cols <- grey(seq(0.2, 0.8, length=length(intervals)))
    
    the.cols <- cols[findInterval(aux.var, intervals)]
    plot(xx.polylist, col=the.cols, ...)
    
    ys <- sort(seq(legend.y[1], legend.y[2], len=length(intervals)))
    
    image(range(legend.x), ys,
          rbind(intervals, intervals), col=cols, add=T)
    
    text(min(legend.x), ys, intervals, pos=2, cex=0.9)
  }
  
  invisible()
}
########################################################




CAR.precision <- function (rho, kappa, W, wts = rep(1, nrow(W))) {
  ## ======================================================================
  ## Calculate the precision matrix for a CAR model
  ## with parameters 'rho', 'kappa',
  ## and proximity matrix 'W', and weights 'wts'.
  ## ======================================================================
  
  diag(wts) %*% (diag(nrow(W)) - rho * W) * kappa
}



CAR.m2l.rho <- function (rho, y, W, wts = rep(1, nrow(W)),
                           X=cbind(rep(1, length(y))), penalty=0,
                           debug=FALSE) {
  ## ======================================================================
  ## Calculate minus two times the profile log likelihood plus a penalty
  ## for alpha
  ## ======================================================================
  
  if (debug)
    cat("rho: ", rho, "\n")
  
  n      <- length(y)
  prec   <- CAR.precision(rho, 1, W, wts)
  A      <- crossprod(X, prec)
  z      <- y - X %*% solve(A %*% X) %*% A %*% y
  
  drop(n * log(2*pi) + n * log((crossprod(z, prec) %*% z)/n) -
         log(det(prec)) + penalty)
}



CAR.mle.rho <- function (y, W, wts = rep(1, nrow(W)),
                           X=cbind(rep(1, length(y))),
                           debug=FALSE) {
  ## ======================================================================
  ## Calculate the MLE of alpha.
  ## ======================================================================
  
  optimize(CAR.m2l.rho, CAR.rho.range(W, wts),
           y=y, W=W, X=X, wts=wts, debug=debug)$min
}


CAR.rho.range <- function (W, wts = rep(1, nrow(W))) {
  ## ======================================================================
  ## Calculate the range of a possible alpha values for a CAR model
  ## with proximity matrix 'W' and weights 'wts'.
  ## ======================================================================
  
  1/range(eigen(W)$values)
}




CAR.mle.beta <- function (rho, y, W, wts = rep(1, nrow(W)),
                          X=cbind(rep(1, length(y)))) {
  ## ======================================================================
  ## Calculate the MLE of beta
  ## ======================================================================
  
  A <- crossprod(X, CAR.precision(rho, 1, W, wts))
  
  solve(A %*% X) %*% A %*% y
}



CAR.mle.kappa <- function (rho, y, W, wts = rep(1, nrow(W)),
                          X=cbind(rep(1, length(y)))) {
  ## ======================================================================
  ## Calculate the MLE of tau
  ## ======================================================================
  
  prec   <- CAR.precision(rho, 1, W, wts)
  A      <- crossprod(X, prec)
  mu.hat <- X %*% solve(A %*% X) %*% A %*% y
  z      <- y - as.numeric(mu.hat)
  
  1/drop( (crossprod(z, prec) %*% z)/length(y) )
}


CAR.cov.mle.beta <- function (rho, kappa, W, wts = rep(1, nrow(W)),
                              X=cbind(rep(1, nrow(W)))) {
  ## ======================================================================
  ## Calculate the estimated covariance matrix for the MLE of beta.
  ## ======================================================================
  
  solve(t(X) %*% CAR.precision(rho, kappa, W, wts) %*% X)
}



CAR.mle <- function (y, W, wts = rep(1, nrow(W)),
                     X=cbind(rep(1, nrow(W)))) {
  ## ======================================================================
  ## Estimate the parameters of the CAR model based on the data 'y',
  ## proximity matrix 'W', weights 'wts', and design matrix 'X'.
  ## ======================================================================
  
  
  rho.hat <- CAR.mle.rho(y, W, wts, X=X)
  kappa.hat  <- CAR.mle.kappa(rho.hat, y, W, wts, X=X)
  beta.hat  <- CAR.mle.beta(rho.hat, y, W, wts, X=X)
  
  cov.beta.hat <- CAR.cov.mle.beta(rho.hat, kappa.hat, W, wts, X=X)
  
  npars <- length(beta.hat) + 2
  
  the.loglik <- -0.5*CAR.m2l.rho(rho.hat, y, W, wts, X=X)
  
  structure(list(y=y, W=W, wts=wts, X=X,
                 rho.hat = rho.hat,
                 kappa.hat  = kappa.hat,
                 beta.hat  = beta.hat,
                 cov.beta.hat = cov.beta.hat,
                 loglik = the.loglik,
                 AIC = -2*the.loglik + 2*npars),
            class="CAR")
}




summary.CAR <- function (x) {
  ## ======================================================================
  ## Summarize this CAR model
  ## ======================================================================
  
  est  <- x$beta.hat
  se   <- sqrt(diag(x$cov.beta.hat))
  tval <- est/se
  error.df <- length(x$y) - length(est)
  
  ans <- list()
  ans$coef <- cbind(est, se, tval, 2 * pt(abs(tval),
                                          error.df, lower.tail = FALSE))
  dimnames(ans$coef)[[2]] <-
    c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  
  ans$rho.hat <- x$rho.hat
  ans$kappa.hat  <- x$kappa.hat
  ans$error.df <- error.df
  ans$loglik    <- x$loglik
  ans$AIC       <- x$AIC
  ans
}



