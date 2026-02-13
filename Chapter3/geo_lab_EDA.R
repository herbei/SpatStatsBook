

rm(list=ls())

library(geoR)
library(maps)
library(ggplot2)
library(dplyr)
library(latex2exp)
library(gridExtra)
library(sp)

My_Theme = theme(
     #axis.title.x = element_text(size = 16),
     axis.text.x = element_text(size = 4),
     #axis.text.x=element_blank(),
     #axis.title.x=element_blank(),
     #axis.text.y=element_blank(),
     #axis.title.y=element_blank(),
     axis.text.y = element_text(size = 4),
     #axis.title.y = element_text(size = 16),
     #axis.ticks.x=element_blank(),
     #axis.ticks.y=element_blank()
)


# Plot a map of Chesapeake Bay
plot_pred = ggplot()+
     annotation_map(map_data("state"),fill = "white", colour = "black") +
     coord_quickmap(xlim=c(-77.5,-75.5), ylim=c(36.6,40)) +
     theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
     xlab("Longitude")+
     ylab("Latitude")
#  scale_x_continuous(expand=c(0,0)) +
#  scale_y_continuous(expand=c(0,0))


# Chesapeake Bay  
pdf("fig_03_060.pdf")
print(plot_pred)
dev.off()



ches <- read.table("chesbay.txt", header=TRUE)
names(ches)
boundary <- read.table("boundary.txt", header=TRUE)

pdf(file="fig_03_061.pdf", width=4, height=6)
plot(boundary, xlab="Easting", ylab="Northing", col="gray40", type="l")
points(ches$easting, ches$northing, pch=16, cex=1)
dev.off()




## Show the boundaries in gray.
pdf(file="fig_03_0611.pdf", width=4, height=6)
plot(boundary, xlab="Easting", ylab="Northing", col="gray40", type="l")
## plot the nitrogen values
text(ches$easting, ches$northing, round(ches$nitrogen,1), cex=0.7)
dev.off()



plot1=ggplot()+
     geom_point(aes(x=ches$easting, y=ches$nitrogen), size=0.5)+
     xlab(("Easting"))+
     ylab("Nitrogen")+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() )+
     My_Theme


plot2=ggplot()+
     geom_point(aes(x=ches$northing, y=ches$nitrogen), size=0.5)+
     xlab("Northing")+
     ylab("Nitrogen")+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
     My_Theme


pdf(file="fig_03_062.pdf", width=6, height=3)
print(grid.arrange(plot1, plot2, nrow=1))
dev.off()




log.nitro <- log(ches$nitrogen)


plot1=ggplot()+
     geom_point(aes(x=ches$easting, y=log.nitro), size=0.5)+
     xlab(("Easting"))+
     ylab("log(Nitrogen)")+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() )+
     My_Theme


plot2=ggplot()+
     geom_point(aes(x=ches$northing, y=log.nitro), size=0.5)+
     xlab("Northing")+
     ylab("log(Nitrogen)")+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
     My_Theme


pdf(file="fig_03_063.pdf", width=6, height=3)
print(grid.arrange(plot1, plot2, nrow=1))
dev.off()





############################################
trend.model <- lm(log.nitro ~ ches$northing)
summary(trend.model)

trend.resids <- resid(trend.model)




pdf(file="fig_03_064.pdf", width=8, height=8)
p1 = ggplot()+
     geom_text(
          data = transform(ches, lab = round(trend.resids, 2)),
          aes(x = longitude, y = latitude, label = lab),
          size = 2.5
     ) +
     theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() )+
     xlab("Longitude")+
     ylab("Latitude")

qq_df <- data.frame(resid = trend.resids)
p2 <- ggplot(qq_df, aes(sample = resid)) +
     stat_qq() +
     stat_qq_line() +
     labs(x = "Normal quantiles", y = "Sample quantiles of residuals") +
     theme_bw()

p3 <- ggplot(data.frame(easting = ches$easting, resid = trend.resids),
             aes(x = easting, y = resid)) +
     geom_point() +
     geom_hline(yintercept = 0, linetype = 2) +
     coord_cartesian(ylim = c(-0.4, 0.4)) +
     labs(x = "easting", y = "residuals") +
     theme_bw()
p4 <- ggplot(data.frame(northing = ches$northing, resid = trend.resids),
             aes(x = northing, y = resid)) +
     geom_point() +
     geom_hline(yintercept = 0, linetype = 2) +
     coord_cartesian(ylim = c(-0.4, 0.4)) +
     labs(x = "northing", y = "residuals") +
     theme_bw()

print(grid.arrange(p1, p2, p3, p4, nrow=2))
dev.off()








################################################################################
my.resids <- cbind(ches$easting, ches$northing, trend.resids)
Gnitro <- as.geodata(my.resids, coords.col=1:2, data.col=3)

Euclidean.dist.matrix <- function (sites) {
     ## ======================================================================
     ## Purpose: The function calculate the Euclidean distances among sites.
     ## Assumes: 'sites' are matrices or data frames.
     ## ======================================================================
     
     dd <- as.matrix(dist(sites, upper=TRUE, diag=TRUE))
     dimnames(dd) <- NULL
     dd
}

## calculate the Euclidean distances between the sites.
dists <- Euclidean.dist.matrix(Gnitro$coords)

## Summarize the distances
hist(dists)
summary(as.numeric(dists))

## Estimate the semivariogram
emp.var.cloud <- variog(Gnitro, option="cloud", max.dist=150)
emp.var.rob <- variog(Gnitro, estimator="modulus", max.dist=150)
plot1=ggplot()+
     geom_point(aes(x=emp.var.cloud$u, y=emp.var.cloud$v), size=0.5)+
     xlab(("Distance"))+
     ylab("Variogram cloud")+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank() )+
     My_Theme


plot2=ggplot()+
     geom_point(aes(x=emp.var.rob$u, y=emp.var.rob$v), size=0.5)+
     xlab("Distance")+
     ylab("Robust binned semi-variogram")+
     theme_bw()+
     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+
     My_Theme


pdf("fig_03_070.pdf", width=8, height=3)
print(grid.arrange(plot1, plot2, nrow=1))
dev.off()

pdf("fig_03_071.pdf", width=8, height=4)
emp.var4 <- variog4(Gnitro, estimator.type="modulus", max.dist=150)
plot(emp.var4, lwd=2.5)
dev.off()


Rmy.resids <- cbind(2.5 * ches$easting, ches$northing, trend.resids)
Rnitro <- as.geodata(Rmy.resids, coords.col=1:2, data.col=3)

## Using directional variogram to assess isotropy of the data
emp.var4 <- variog4(Rnitro, estimator.type="modulus", max.dist=150)

pdf("fig_03_072.pdf", width=8, height=4)
plot(emp.var4, lwd=2.5)
dev.off()


pdf("fig_03_073.pdf", width=8, height=4)
emp.var.final <- variog(Rnitro, estimator="modulus", max.dist=150)
plot(emp.var.final, pch=16, cex=1)
dev.off()









