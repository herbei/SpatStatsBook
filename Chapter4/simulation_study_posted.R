rm(list=ls())
library(mvtnorm) 
library(ggplot2) 
library(gridExtra)


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




n = 30
x = seq(0, 1, length.out=n)
y = seq(0, 1, length.out=n)
xy = expand.grid(x=x, y=y)


##################################
# Plot various draws from a GP process with different parameter values
##################################
sigma2=0.5
tau2=0.1
phi=1


MUvec <- 0*(xy[,1] + xy[,2])
COVmat <- matrix(NA, nrow = n^2, ncol = n^2) 
for(i in 1:(n^2)){
  for(j in 1:(n^2)){
    h=sqrt( (xy[i,1]-xy[j,1])^2 + (xy[i,2]-xy[j,2])^2 )
    COVmat[i, j] = exp.cov(h,sigma2,tau2,phi) 
  } 
}

Z = MASS::mvrnorm(1, MUvec, COVmat)
pp <- data.frame(Z=Z,x1=xy[,1],x2=xy[,2])

p1 = ggplot(pp,aes(x=x1,y=x2)) +
  geom_raster(aes(fill=Z), interpolate = TRUE) +
  geom_contour(aes(z=Z), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=0, low="blue", mid="green", high="red") 


pdf("figure_name_here.pdf", width=5, height=4)
print(p1)
dev.off()
#######################################


#################################################################################
### A simulation study
##################################################################################
sigma2=0.5
tau2=0.5
phi=1

#################
## Simulate data
################
MUvec <- 0*(xy[,1] + xy[,2])
COVmat <- matrix(NA, nrow = n^2, ncol = n^2) 
for(i in 1:(n^2)){
  for(j in 1:(n^2)){
    h=sqrt( (xy[i,1]-xy[j,1])^2 + (xy[i,2]-xy[j,2])^2 )
    COVmat[i, j] = exp.cov(h,sigma2,tau2,phi) 
  } 
}

Z = MASS::mvrnorm(1, MUvec, COVmat)

#########################
## Subsample the field Z
#########################
NN  = n^2
N=200
idx_d = sample(seq(1:NN),size=N)

Z_d=Z[idx_d]


################################
## Plot the full field and the subsampled one
###############################
pp <- data.frame(Z=Z,x1=xy[,1],x2=xy[,2])

p1 = ggplot(pp,aes(x=x1,y=x2)) +
  geom_raster(aes(fill=Z), interpolate = TRUE) +
  geom_contour(aes(z=Z), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=0, low="blue", mid="green", high="red") 


pdf("figure_name_here.pdf", width=5, height=4)
print(p1)
dev.off()


###############
## Variog function
###############
tt=seq(0,1.4,by=0.01)
VV=rep(0, len=length(tt))
for(i in 1:length(tt)){
  VV[i] = exp.variog(tt[i], sigma2, tau2, phi)
}


#Z = rmvnorm(1, mean = MUvec, sigma = COVmat)
Z=t(Z) # make it a column vector




library(geoR)
sites <- cbind(xy[idx_d,1], xy[idx_d,2])


## Estimate and draw the variogram
vg <- variog(coords=sites, data=Z_d, estimator.type = "modulus")
reml <- likfit(data=Z_d, coords = sites, cov.model="exponential", lik.method="REML", fix.nugget = FALSE, nugget = 0.3, ini.cov.pars = c(0.8,0.3))


#################################
## Plot estimated variogram
#################################
pdf(file="figure_name_here.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,5), xlim=c(-0.0001,1.5), xlab="distance", ylab="variogram", cex=1)
lines(reml, col="blue")
lines(tt[2:length(tt)],VV[2:length(tt)], col="red")
points(tt[1], VV[1], pch=16, col="red")
dev.off()
#######################################################################






