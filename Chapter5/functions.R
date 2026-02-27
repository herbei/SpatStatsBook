
######################################################
## Exp covariance and variogram
######################################################
exp.cov <- function(h, sigma2, tau2, phi){
  if(h==0){
    out=sigma2+tau2
  }
  else{
    out=sigma2 * exp(-h/phi)
  }
  out
}



exp.variog <- function(h, sigma2, tau2, phi){
  if(h==0){
    out=0
  }
  else{
    out=tau2 + sigma2 * (1- exp(-h/phi))
  }
  out
}
###################################################

######################################################
## PowExp covariance and variogram
######################################################
pow.exp.cov <- function(h, sigma2, tau2, phi, nu){
  if(h==0){
    out=sigma2+tau2
  }
  else{
    out=sigma2 * exp(-(h^nu)/phi)
  }
  out
}



pow.exp.variog <- function(h, sigma2, tau2, phi, nu){
  if(h==0){
    out=0
  }
  else{
    out=tau2 + sigma2 * (1- exp(-(h^nu)/phi))
  }
  out
}
###################################################
