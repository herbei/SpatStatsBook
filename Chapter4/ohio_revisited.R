rm(list=ls())
setwd("~/Dropbox/0Work/0Teach/6530/OLDNotes/Figures/Ch4/")
ohio_temps = read.csv('ohio_temps.csv', header = TRUE)

#remove the outlier
ohio_temps = ohio_temps[-13,]

myLON = ohio_temps$LONGITUDE
myLAT = ohio_temps$LATITUDE
myTemp = ohio_temps$TMAX
myLON2 = myLON^2
myLAT2 = myLAT^2

model1 = lm(myTemp ~ myLAT + myLON + myLON2 + myLAT2 + myLON * myLAT)

## Calculate the residuals from the linear model
res1 <- resid(model1)
## Calculate the fitted values from the linear model
fits1 <- fitted(model1)




library(geoR)
sites <- cbind(myLON, myLAT)
## Estimate and draw the semivariogram cloud
vg <- variog(coords=sites, data=res1, estimator.type = "modulus",  max.dist = 3)


########## Weighted least squares ##############################
wls <- variofit(vg, cov.model="powered.exponential", fix.nugget = FALSE, nugget = 3, ini.cov.pars = c(0.5,0.5), fix.kappa = FALSE, kappa = 1.9, weights="cressie")
pdf(file="fig_04_200.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,6.1), xlab="distance", ylab="variogram", cex=1)
lines(wls)
abline(h=wls$nugget, lty = 2, col="blue")
abline(h=wls$nugget+wls$cov.pars[1], lty=2, col="blue")
dev.off()
##############################################################

########## Weighted least squares - Matern covariance #######
wls1 <- variofit(vg, cov.model="matern", fix.nugget = FALSE, nugget = 3, ini.cov.pars = c(0.5,0.5), fix.kappa = FALSE, kappa = 1.9, weights="cressie")
pdf(file="fig_04_201.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,6.1), xlab="distance", ylab="variogram", cex=1)
lines(wls1, col="red")
abline(h=wls1$nugget, lty = 2, col="blue")
abline(h=wls1$nugget+wls1$cov.pars[1], lty=2, col="blue")
dev.off()
##############################################################

######## Plot both WLS fits
pdf(file="fig_04_202.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,6.1), xlab="distance", ylab="variogram", cex=1)
lines(wls, col="black")
lines(wls1, col="red")
dev.off()
############################################



########## MLE estimation #############################
mle <- likfit(data=res1, coords = sites, cov.model="powered.exponential", fix.nugget = FALSE, nugget = 3.6, ini.cov.pars = c(0.7,0.3), fix.kappa = FALSE, kappa = 1.9, lik.method="ML")
summary(mle)
pdf(file="fig_04_210.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,6.1), xlab="distance", ylab="variogram", cex=1)
lines(wls)
lines(mle, col="blue")
dev.off()
#######################################################################

########## MLE estimation #############################
mle1 <- likfit(data=res1, coords = sites, cov.model="matern", fix.nugget = FALSE, nugget = 2.5, ini.cov.pars = c(0.7,0.3), fix.kappa = FALSE, kappa = 1, lik.method="ML")
pdf(file="fig_04_211.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,6.1), xlab="distance", ylab="variogram", cex=1)
lines(wls, col="black")
lines(mle, col="blue")
lines(mle1, col="red")
dev.off()
#######################################################################






############# MLEs and REML ##########################################
reml1 <- likfit(data=res1, coords = sites, cov.model="powered.exponential", lik.method="REML", fix.nugget = FALSE, nugget = 3.6, ini.cov.pars = c(0.8,0.3), fix.kappa = FALSE, kappa = 1.9)
reml2 <- likfit(data=res1, coords = sites, cov.model="spherical", lik.method="REML", fix.nugget = FALSE, nugget = 1, ini.cov.pars = c(0.8,0.3))
reml3 <- likfit(data=res1, coords = sites, cov.model="exp", lik.method="REML", fix.nugget = FALSE, nugget = 1, ini.cov.pars = c(0.8,0.3))

pdf(file="fig_04_210.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,6.1), xlab="distance", ylab="variogram", cex=1)
#lines(wls)
lines(mle, col="blue")
lines(reml1, col="red")
lines(reml2, col="purple")
lines(reml3, col="green")
dev.off()
#######################################################################








#mydata = cbind(res1, sites)
#write.csv(mydata, file="myresiduals.csv")
