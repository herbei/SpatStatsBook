rm(list=ls())
setwd("~/Dropbox/0Work/0Teach/6530/Notes/Chapter3/")
library("ggplot2")
library("plgp")
library(latex2exp)
library(gridExtra)
library(maps)

ohio_temps = read.csv('ohio_temps.csv', header = TRUE)

# what variables do we work with 
names(ohio_temps)

myLON = ohio_temps$LONGITUDE
myLAT = ohio_temps$LATITUDE
myTemp = ohio_temps$TMAX
myELEV = ohio_temps$ELEVATION



######################################################
pdf(file="fig_03_010.pdf", width=6, height=6)
par(mfrow=c(1,1), cex=0.75, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")

map("state", "ohio", resolution=0, xlim=c(-85.5,-80), ylim=c(38,42.5))
map.axes()
points(myLON, myLAT, type = "n", xlab = "Lon", ylab = "Lat")
text(myLON, myLAT, round(myTemp, 1))
mtext(c("Longitude", "Latitude", "Max Temp, 01/01/22"), side=c(1,2,3), line = 2.5, cex=0.75)
dev.off()
######################################################




######################################################
pdf(file="fig_03_15.pdf", width=4, height=4)
par(mfrow=c(1,1), cex=0.75, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(myLON, myTemp, xlab="Longitude", ylab="Max Temperature", pch = 16, cex=1.2)
dev.off()


pdf(file="fig_03_20.pdf", width=4, height=4)
par(mfrow=c(1,1), cex=0.75, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(myLAT, myTemp, xlab="Latitude", ylab="Max Temperature", pch = 16, cex = 1.2)
dev.off()

pdf(file="fig_03_25.pdf", width=4, height=4)
par(mfrow=c(1,1), cex=0.75, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(myELEV, myTemp, xlab="Elevation", ylab="Max Temperature", pch = 16, cex = 1.2)
#text(myELEV, myTemp, ohio_temps$NAME, 1)
dev.off()
#######################################################




myLON2 = myLON^2

model1 = lm(myTemp ~ myLAT + myLON + myLON2)
summary(model1)

## Calculate the residuals from the linear model
res1 <- resid(model1)

## Calculate the fitted values from the linear model
fits1 <- fitted(model1)


#################################
pdf(file="fig_03_30.pdf", width=8, height=4)
par(mfrow=c(1,2), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")

map("state", "ohio", myborder=0, xlim=c(-85.5,-80), ylim=c(38,42.5))
map.axes()
text(myLON, myLAT, round(res1, 1))
mtext(c("Longitude", "Latitude"), side=c(1,2), line = 2.5, cex=0.5)


qqnorm(res1, main="", xlab="Normal quantiles")
qqline(res1)
dev.off()
###############################################



#################################
pdf(file="fig_03_35.pdf", width=8, height=4)
par(mfrow=c(1,2), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")

plot(myLON, res1, xlab="Latitude", ylab="Residual")
abline(h=0, lty=2)

plot(myLAT, res1, xlab="Longitude", ylab="Residual")
abline(h=0, lty=2)

dev.off()
################################




model2 = lm(myTemp ~ myLAT + myLON + myLON2 + myELEV + myLON * myLAT)
summary(model2)

## Calculate the residuals from the linear model
res2 <- resid(model2)

## Calculate the fitted values from the linear model
fits2 <- fitted(model2)


#################################
pdf(file="fig_03_40.pdf", width=8, height=4)
par(mfrow=c(1,2), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")

map("state", "ohio", myborder=0, xlim=c(-85.5,-80), ylim=c(38,42.5))
map.axes()
text(myLON, myLAT, round(res2, 1))
mtext(c("Longitude", "Latitude"), side=c(1,2), line = 2.5, cex=0.5)


qqnorm(res2, main="", xlab="Normal quantiles")
qqline(res2)
dev.off()
###############################################


#################################
pdf(file="fig_03_41.pdf", width=8, height=4)
par(mfrow=c(1,2), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")

plot(myLON, res2, xlab="Latitude", ylab="Residual")
abline(h=0, lty=2)

plot(myLAT, res2, xlab="Longitude", ylab="Residual")
abline(h=0, lty=2)

dev.off()
################################


library(geoR)


sites <- cbind(myLON, myLAT)
## Estimate and draw the semivariogram cloud
vg.cloud <- variog(coords=sites, data=res2, option="cloud")

################################################################
pdf(file="fig_03_42.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg.cloud)
dev.off()
################################################################


## Estimate and draw the binned empirical semivariogram
vg <- variog(coords=sites, data=res2)


################################################################
pdf(file="fig_03_43.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, pch=16, ylim=c(0,7), xlab="distance", ylab="semivariogram", cex=2)
dev.off()
################################################################

## Estimate the semivariogram Cressie and Hawkins
vg.CH <- variog(coords=sites, data=res2, estimator.type = "modulus", max.dist = 3)


################################################################
pdf(file="fig_03_55.pdf", width=7, height=3)
par(mfrow=c(1,1), cex=0.5, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg.CH$u, vg.CH$v, pch=16, ylim=c(0,7), xlab="distance", ylab="semivariogram", cex=2)
dev.off()
################################################################


vg_01 = variog(coords = sites, data=res2, direction = 0)
vg_02 = variog(coords = sites, data=res2, direction = pi/4)
vg_03 = variog(coords = sites, data=res2, direction = pi/2)
vg_04 = variog(coords = sites, data=res2, direction = 3*pi/4)


##########################################################
pdf(file="fig_03_50.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg_01$u, vg_01$v, type="l", ylim=c(0,20), col="black", xlab = "distance", ylab="semivariogram")
lines(vg_02$u, vg_02$v, type="l", ylim=c(0,20), col="red")
lines(vg_03$u, vg_03$v, type="l", ylim=c(0,2), col="green")
lines(vg_04$u, vg_04$v, type="l", ylim=c(0,20), col="blue")
dev.off()
##########################################################



