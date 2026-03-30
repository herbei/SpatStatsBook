# Acknowledgement: most of this code was written by Peter Craigmile
###################################################################

rm(list = ls())
source("CAR_RH.R")



## Read the SAT dataset
sat <- read.table("sat.txt", header=TRUE)
colnames(sat)[5] = "region"

## What are the variable names?
names(sat)

## Load in the shape file for the map.
load("sat_US_states.RData")

## Match the state names in the dataset to the shape file
## where.is.state <- pmatch(tolower(US$STATE_NAME), tolower(sat$name))

## Read in the proximity matrix that we will use
W <- as.matrix(read.table("sat_proximity.txt", header=F))



####################################################################
## Produce a Chloropleth map of the SAT verbal score.
us_states <- map_data("state") %>% 
  select(lon = long, lat, group, id = subregion, region)

merged_data=inner_join(us_states, sat, by = "region")

sat.plot = ggplot() +
  geom_polygon(data=merged_data, aes(x=lon, y=lat, group=group, fill = vscore), colour = "white", linewidth=0.2) + 
  coord_quickmap()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))



pdf("fig_06_700.pdf")
print(sat.plot)
dev.off()
####################################################################



#####################################################################
sat.plot1 = ggplot()+
  geom_point(data=sat, aes(x=pc, y=vscore)) +
  xlab("Percentage taking the exam")+
  ylab("Verbal score")+

pdf("fig_06_701.pdf", width=6, height=4)
print(sat.plot1)
dev.off()
####################################################################





#######################################################################
## Fit the ordinary least squares model predicting the verbal
## score from a quadratic of the percent taking the exam.
ols.model <- lm(vscore ~ pc + I(pc^2), data=sat)
summary(ols.model)

## Add the fitted line from the ordinary least squares fit to the scatterplot
## -- [order(sat$pc)] reorders the fitted values to be in increasing order of
## the precentage.
xx=sort(sat$pc)
yy=fitted(ols.model)[order(sat$pc)]

sat.plot2 = ggplot()+
  geom_point(data=sat, aes(x=pc, y=vscore)) +
  geom_line(aes(x=xx, y=yy))+
  xlab("Percentage taking the exam")+
  ylab("Verbal score")

pdf("fig_06_702.pdf", width=6, height=4)
print(sat.plot2)
dev.off()
#######################################################################


## Calculate and plot the residuals
res <- resid(ols.model)

sat1=cbind(sat, res)
merged_data1=inner_join(us_states, sat1, by = "region")
sat.plot3 = ggplot() +
  geom_polygon(data=merged_data1, aes(x=lon, y=lat, group=group, fill = res), colour = "white", linewidth=0.2) + 
  coord_quickmap()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))



pdf("fig_06_703.pdf")
print(sat.plot3)
dev.off()
####################################################################




## Check the normality assumption using a normal Q-Q plot
pdf("fig_06_704.pdf", width=4, height=4)
par(mfrow=c(1,1), cex=0.75, mar=c(4,4,1,1), mgp=c(2,0.5,0), bty="L")

qqnorm(res)
qqline(res)
dev.off()
#####################################################################




## Create the row stochastic proximity matrix
W.tilde <- W/rowSums(W)
## Create the vector of weights
the.weights <- rowSums(W)
## Create the design matrix
X <- model.matrix(~ pc + I(pc^2), data=sat)


## Fit the CAR model
car.model <- CAR.mle(sat$vscore, W.tilde, the.weights, X=X)


## First, what are the plausible values of rho for the CAR model?
## [1] -1.392427  1.000000
CAR.rho.range(car.model$W, car.model$wts)


## Now summarize the model
summary(car.model)





