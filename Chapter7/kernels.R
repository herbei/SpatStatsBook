rm(list=ls())
setwd('~/Dropbox/0Work/0Teach/6530/Notes/Figures/Ch7/')


KK = function(z){
  
  # Uniform
  # out = dunif(z, min=-1, max=1)
  
  # Gaussian
  # out = dnorm(z,0,1)
  
  # Epanechnikov
   out = (3/4) * max(1-z^2, 0)
  
  
}

xx = seq(-5,5, by=0.01)
yy=numeric(length(xx))
for(i in 1:length(xx)){
  yy[i] = KK(xx[i])
}

plot = ggplot() + 
  geom_line(aes(x=xx,y=yy), lwd =0.5) + 
  xlab("z") + 
  ylab("K(z)")

pdf("Fig_07_102.pdf", width=3, height=2)
print(plot)
dev.off()


N=50
X = rnorm(N, mean=1, sd=1)

bw=0.5
f_est1 = numeric(length(xx))
for(i in 1:length(xx)){
  
  S=0
  for(j in 1:N){
    S = S + KK((xx[i] - X[j])/bw)
  }
  S = S/(N*bw)
  
  f_est1[i]=S
}

bw=0.3
f_est2 = numeric(length(xx))
for(i in 1:length(xx)){
  
  S=0
  for(j in 1:N){
    S = S + KK((xx[i] - X[j])/bw)
  }
  S = S/(N*bw)
  
  f_est2[i]=S
}

bw=0.05
f_est3 = numeric(length(xx))
for(i in 1:length(xx)){
  
  S=0
  for(j in 1:N){
    S = S + KK((xx[i] - X[j])/bw)
  }
  S = S/(N*bw)
  
  f_est3[i]=S
}


pdf(file="fig_07_104.pdf", width=4, height=3.5)
plot(X, rep(0.01, length=N), pch=16, cex=1, xlim=c(-2,4), ylim=c(-0.1,1), ylab="", main="")
dev.off()


f_est1 = density(X, bw=0.5, kernel="epanechnikov")
f_est2 = density(X, bw=0.3, kernel="epanechnikov")
f_est3 = density(X, bw=0.05, kernel="epanechnikov")

pdf(file="fig_07_107.pdf", width=7, height=3.5)
par(mfrow=c(1,3), cex=0.8, mar=c(3.5,3.5,2,0), mgp=c(2,0.5,0), bty="L")
plot(f_est1, main="")
plot(f_est2, main="")
plot(f_est3, main="")
dev.off()


plot(X, rep(0.01, length=N), pch=16, cex=0.6, xlim=c(-2,4), ylim=c(-0.1,0.6),xlab="z, bw=0.5", ylab = "f_hat(z)" , main="")
lines(xx, f_est1)
plot(X, rep(0.01, length=N), pch=16, cex=0.6, xlim=c(-2,4), ylim=c(-0.1,0.6),xlab="z, bw=0.3", ylab = "f_hat(z)" , main="")
lines(xx, f_est2)
plot(X, rep(0.01, length=N), pch=16, cex=0.6, xlim=c(-2,4), ylim=c(-0.1,0.6),xlab="z, bw=0.05", ylab = "f_hat(z)" , main="")
lines(xx, f_est3)
dev.off()


###############################################################################################


N=50
XX = numeric(N)
for(i in 1:N){
  if( runif(1,0,1)<0.3){
    XX[i] = rnorm(1,0,1)
  }  
  else{
    XX[i] = rnorm(1,3,0.5)
  }
}
f_est1 = density(XX, bw=1)
f_est2 = density(XX, bw=0.1)
f_est3 = density(XX, bw=0.5)
f_est4 = density(XX, bw=0.3)


pdf(file="fig_07_108.pdf", width=6, height=4)
par(mfrow=c(2,2), cex=0.7, mar=c(3,4,4,1), mgp=c(2,0.5,0), bty="L")
plot(f_est1, type = "l", main="")
plot(f_est3, type = "l", main="")
plot(f_est4, type = "l", main="")
plot(f_est2, type = "l", main="")
dev.off()
 #######################################################################


N=100
X = rexp(N, rate=1)
f_est = density(X, bw=0.5)
pdf(file="fig_07_109.pdf", width=4, height=3.5)
par(mfrow=c(1,1), cex=0.7, mar=c(3,4,4,1), mgp=c(2,0.5,0), bty="L")
plot(f_est, type = "l", main="")
dev.off()


Y = log(X)
fy_est = density(Y, bw=0.8)

xx = exp(fy_est$x)
yy = (fy_est$y)/xx

pdf(file="fig_07_110.pdf", width=8, height=3.5)
par(mfrow=c(1,2), cex=0.7, mar=c(3,4,4,1), mgp=c(2,0.5,0), bty="L")
plot(fy_est, type = "l", main="Density of log(X)")
plot(xx[100:512],yy[100:512], type = "l", main = "density of X", xlim=c(-1,10))
lines(xx,dexp(xx), col="red")
dev.off()












