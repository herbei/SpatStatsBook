rm(list=ls())
ohio_temps = read.csv('ohio_temps.csv', header = TRUE)

#remove the outlier
ohio_temps = ohio_temps[-13,]


myLON = ohio_temps$LONGITUDE
myLAT = ohio_temps$LATITUDE
myTemp = ohio_temps$TMAX
myLON2 = myLON^2
myLAT2 = myLAT^2


mydata = cbind(myLON, myLAT, myTemp, myLON2, myLAT2)
mydata = as.data.frame(mydata)



model1 = lm(myTemp ~ myLAT + myLON + myLON2 + myLAT2 + myLON * myLAT)

## Calculate the residuals from the linear model
res1 <- resid(model1)
## Calculate the fitted values from the linear model
fits1 <- fitted(model1)


library(geoR)
library(maps)
library(ggplot2)
source('functions.R')

sites <- cbind(myLON, myLAT)
pred.sites = read.csv('krig_sites.csv', header = TRUE)
pred.sites=pred.sites[,-1]

pred_data=cbind(res1,myLON, myLAT)
pred_data = as.data.frame(pred_data)




#########################################################################
plot_res = ggplot(pred_data, aes(x=myLON, y=myLAT, color=res1)) + 
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(size=1.2)+
  geom_point(data=pred.sites,aes(x=p_lon,y=p_lat),size=0.2, colour="grey")+
  #geom_point(data=pred.sites) +
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint = 0, low="blue", mid="green", high="red") +
  #scale_color_gradientn(colours = c("#0000FF","#00FF00", "#FF0000" ))+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")

pdf("fig_05_100.pdf", width=4, height=4)
print(plot_res)
dev.off()
#######################################################################




#############################
## Fit two cov models via REML
#############################
vg <- variog(coords=sites, data=res1, estimator.type = "modulus",  max.dist = 3)

reml1 <- likfit(data=res1, coords = sites, cov.model="powered.exponential", lik.method="REML", fix.nugget = FALSE, nugget = 3.6, ini.cov.pars = c(0.8,0.3), fix.kappa = FALSE, kappa = 1.9)
reml2 <- likfit(data=res1, coords = sites, cov.model="exponential", lik.method="REML", fix.nugget = FALSE, nugget = 1, ini.cov.pars = c(0.8,0.3))



#############################
## Plot the two fitted models
#############################
pdf(file="fig_05_150.pdf", width=8, height=4)
par(mfrow=c(1,1), cex=1, mar=c(5,5,5,1), mgp=c(1.8,0.5,0), bty="L")
plot(vg$u, vg$v, xlim=c(-0.001, 3), ylim=c(0,5), xlab="distance", ylab="variogram", pch=16, cex=1)
lines(reml1, col="blue")
lines(reml2, col="red")
points(0,0, cex=1, pch=16, col="red")
dev.off()
############################



######################################################
### Simple Kriging - mean assumed constant and known
#####################################################

mykrige1_sk = krige.conv(data=res1, coords=sites, loc=pred.sites, krige = krige.control(type.krige="SK", beta=0, cov.model=reml1$cov.model, cov.pars=reml1$cov.pars, nugget=reml1$nugget))

pred_res1_sk = cbind(mykrige1_sk$predict, pred.sites$p_lon, pred.sites$p_lat)
colnames(pred_res1_sk) = c("pred_res", "lon", "lat")
pred_res1_sk=as.data.frame(pred_res1_sk)

plot_pred_res1 = ggplot()+
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(data=pred_res1_sk, aes(x=lon, y=lat, color=pred_res), size=1.2)+
  geom_text(data=pred_data, aes(x=myLON, y=myLAT,label=round(res1,1)), size=2)+
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint = 0, low="blue", mid="green", high="red") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")

pdf("fig_05_180.pdf", width=4, height=4)
print(plot_pred_res1)
dev.off()

####################################################



########################################################
### Ordinary kriging - mean is unknown, variance parameters are known

mykrige1_ok = krige.conv(data=res1, coords=sites, loc=pred.sites, krige = krige.control(type.krige="OK", cov.model=reml1$cov.model, cov.pars=reml1$cov.pars, nugget=reml1$nugget))
mykrige1_ok$beta.est

mykrige2_ok = krige.conv(data=res1, coords=sites, loc=pred.sites, krige = krige.control(cov.model=reml2$cov.model, cov.pars=reml2$cov.pars, nugget=reml2$nugget))




pred_res1_ok = cbind(mykrige1_ok$predict, pred.sites$p_lon, pred.sites$p_lat)
colnames(pred_res1_ok) = c("pred_res", "lon", "lat")
pred_res1_ok=as.data.frame(pred_res1_ok)

pred_res2_ok = cbind(mykrige2_ok$predict, pred.sites$p_lon, pred.sites$p_lat)
colnames(pred_res2_ok) = c("pred_res", "lon", "lat")
pred_res2_ok=as.data.frame(pred_res2_ok)


plot_pred_res1 = ggplot()+
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(data=pred_res1_ok, aes(x=lon, y=lat, color=pred_res), size=1.2)+
  geom_text(data=pred_data, aes(x=myLON, y=myLAT,label=round(res1,1)), size=2)+
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint = 0, low="blue", mid="green", high="red") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")

plot_pred_res2 = ggplot()+
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(data=pred_res2_ok, aes(x=lon, y=lat, color=pred_res), size=1.2)+
  geom_text(data=pred_data, aes(x=myLON, y=myLAT,label=round(res1,1)), size=2)+
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint = 0, low="blue", mid="green", high="red") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")

pdf("fig_05_201.pdf", width=4, height=4)
print(plot_pred_res1)
dev.off()
pdf("fig_05_202.pdf", width=4, height=4)
print(plot_pred_res2)
dev.off()


mydiff = pred_res1_ok$pred_res - pred_res2_ok$pred_res
pred_diff = cbind(mydiff, pred.sites$p_lon, pred.sites$p_lat)
colnames(pred_diff) = c("diff", "lon", "lat")
pred_diff = as.data.frame(pred_diff)

plot_pred_diff = ggplot()+
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(data=pred_diff, aes(x=lon, y=lat, color=diff), size=1.2) +
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint = 0, low="blue", mid="white", high="red") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")

pdf("fig_05_203.pdf", width=4, height=4)
print(plot_pred_diff)
dev.off()

############################################################################################



#################################################################################
###### Kriging Standard errors
################################################################################

pred_vars1 = cbind(mykrige1_ok$krige.var, pred.sites$p_lon, pred.sites$p_lat)
colnames(pred_vars1) = c("Var", "lon", "lat")
pred_vars1=as.data.frame(pred_vars1)

plot_pred_vars1 = ggplot()+
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(data=pred_vars1, aes(x=lon, y=lat, color=Var), size=1.2)+
  geom_text(data=pred_data, aes(x=myLON, y=myLAT,label=round(res1,1)), size=2)+
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint=3.2,low="blue", mid="pink", high="red") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")


pdf("fig_05_300.pdf", width=4, height=4)
print(plot_pred_vars1)
dev.off()
###############################################################################





##############################################################################
## Predicted temperature = trend mean + predicted residuals
#############################################################################
beta_hat = as.matrix(model1$coefficients, ncol=1)
NN = dim(pred.sites)[1]
trend_temps = numeric(NN)
for(i in 1:NN){
      xx = c(1,pred.sites$p_lat[i], pred.sites$p_lon[i], pred.sites$p_lon[i]^2, pred.sites$p_lat[i]^2, pred.sites$p_lat[i]*pred.sites$p_lon[i])
      xx = as.matrix(xx, ncol=1)
      trend_temps[i] = t(xx) %*% beta_hat
}

trend_temps = cbind(pred.sites, trend_temps)
tred_temps=as.data.frame(trend_temps)
names(trend_temps) = c("lon", "lat", "Trend_temps")


#### Plot the data
#########################################################################
plot_data_temps = ggplot(mydata, aes(x=myLON, y=myLAT, color=myTemp)) + 
      annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
      geom_point(size=1.2)+
      geom_text(data=mydata, aes(x=myLON, y=myLAT+1,label=round(myTemp,1)), size=2, color="black")+
      geom_point(data=pred.sites,aes(x=p_lon,y=p_lat),size=0.2, colour="grey")+
      coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
      scale_color_gradient2(midpoint = 50, low="blue", mid="green", high="red") +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
      xlab("Longitude")+
      ylab("Latitude")

pdf("fig_05_301.pdf", width=4, height=4)
print(plot_data_temps)
dev.off()
#######################################################################



#### Plot the trend field
#######################################################################
plot_trend_temps = ggplot()+
  annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
  geom_point(data=trend_temps, aes(x=lon, y=lat, color=Trend_temps), size=1.2)+
  coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
  scale_color_gradient2(midpoint=50,low="blue", mid="green", high="red") +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
  xlab("Longitude")+
  ylab("Latitude")

pdf("fig_05_302.pdf", width=4, height=4)
print(plot_trend_temps)
dev.off()
###############################################################################



### Predicted temperatures = trend + predicted residuals
pred_temps = numeric(NN)
pred_temps = trend_temps$Trend_temps + mykrige1_sk$predict
pred_temps = cbind(pred.sites$p_lon, pred.sites$p_lat, pred_temps)
pred_temps = as.data.frame(pred_temps)
names(pred_temps) = names(trend_temps) = c("lon", "lat", "Pred_temps")


#### Now plot the entire predicted temperature field.
################################################################################
plot_pred_temps = ggplot()+
      annotation_map(map_data("state"),fill = "white", colour = "darkgrey") +
      geom_point(data=pred_temps, aes(x=lon, y=lat, color=Pred_temps), size=1.2)+
      #geom_text(data=mydata, aes(x=myLON, y=myLAT,label=round(myTemp,1)), size=2)+
      coord_quickmap(xlim=c(-85,-80), ylim=c(38.5,42)) +
      scale_color_gradient2(midpoint=50,low="blue", mid="green", high="red") +
      theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
      xlab("Longitude")+
      ylab("Latitude")


pdf("fig_05_303.pdf", width=4, height=4)
print(plot_pred_temps)
dev.off()
 ###############################################################################




