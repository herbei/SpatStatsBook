rm(list = ls())
setwd('~/Dropbox/0Work/0Teach/6530/Notes/Figures/Ch7/')






pow.exp.cov <- function(h, sigma2, tau2, phi, nu){
  if(h==0){
    out=sigma2+tau2
  }
  else{
    out=sigma2 * exp(-(h^nu)/phi)
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
tau2=0.01
phi=0.5
nu=2

MUvec <- 10+0*(xy[,1] + xy[,2])
COVmat <- matrix(NA, nrow = n^2, ncol = n^2) 
for(i in 1:(n^2)){
  for(j in 1:(n^2)){
    h=sqrt( (xy[i,1]-xy[j,1])^2 + (xy[i,2]-xy[j,2])^2 )
    COVmat[i, j] = pow.exp.cov(h,sigma2,tau2,phi,nu) 
  } 
}

Z = MASS::mvrnorm(1, MUvec, COVmat)
pp <- data.frame(Z=Z,x1=xy[,1],x2=xy[,2])

Z = Z-5;
pp <- data.frame(Z=Z,x1=xy[,1],x2=xy[,2])

p1 = ggplot(pp,aes(x=x1,y=x2)) +
  geom_raster(aes(fill=Z), interpolate = TRUE) +
  geom_contour(aes(z=Z), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=5, low="blue", mid="green", high="red") 



p2 = ggplot(pp,aes(x=x1,y=x2)) +
  geom_raster(aes(fill=exp(Z)), interpolate = TRUE) +
  geom_contour(aes(z=exp(Z)), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=exp(5), low="blue", mid="green", high="red") 


pdf("fig_07_400.pdf", width=5, height=4)
print(p1)
dev.off()

pdf("fig_07_401.pdf", width=5, height=4)
print(p2)
dev.off()

#######################################


LL = exp(Z)
LL
pp1 = data.frame(L = LL, x= xy[,1], y=xy[,2])
C = max(LL)+10


N = rpois(1,C)
XX = matrix(runif(2*N, 0,1), ncol=2)


p3 = ggplot() +
  geom_raster(aes(x = pp1$x, y=pp1$y, fill=pp1$L), interpolate = TRUE) +
  guides(fill=guide_legend(title='Intensity'))+
  #geom_contour(aes(z=exp(Z)), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  geom_point(aes(x=XX[,1], y=XX[,2]), size=1, color="black")+
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=exp(5), low="blue", mid="green", high="red") +
  xlab("X1")+
  ylab("X2")


pdf("fig_07_402.pdf", width=5, height=4)
print(p3)
dev.off()

XX=XX[,c(1,2)]
XX = cbind(XX, rep(0, len=dim(XX)[1]))
for(i in 1:dim(XX)[1]){
  xxx=XX[i,1]
  yyy=XX[i,2]
  AA = pp1[,c(2,3)] - cbind(rep(xxx, dim(pp1)[1]), rep(yyy, dim(pp1)[1]))
  ii=which.min(rowSums(AA^2))
  
  my.prob = pp1[ii,1]/C
  
  if(runif(1,0,1) < my.prob){
    XX[i,3] = 1
  }
  
  
}

idx=XX[,3]==1



p4 = ggplot() +
  geom_raster(aes(x = pp1$x, y=pp1$y, fill=pp1$L), interpolate = TRUE) +
  guides(fill=guide_legend(title='Intensity'))+
  #geom_contour(aes(z=exp(Z)), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  geom_point(aes(x=XX[idx,1], y=XX[idx,2]), size=1, color="black")+
  geom_point(aes(x=XX[!idx,1], y=XX[!idx,2]), size=1, color="red")+
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=exp(5), low="blue", mid="green", high="red") +
  xlab("X1")+
  ylab("X2")


pdf("fig_07_403.pdf", width=5, height=4)
print(p4)
dev.off()

p5 = ggplot() +
  geom_raster(aes(x = pp1$x, y=pp1$y, fill=pp1$L), interpolate = TRUE) +
  guides(fill=guide_legend(title='Intensity'))+
  #geom_contour(aes(z=exp(Z)), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  geom_point(aes(x=XX[idx,1], y=XX[idx,2]), size=1, color="black")+
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=exp(5), low="blue", mid="green", high="red") +
  xlab("X1")+
  ylab("X2")

p6 = ggplot() +
#  geom_raster(aes(x = pp1$x, y=pp1$y, fill=pp1$L), interpolate = TRUE) +
  guides(fill=guide_legend(title='Intensity'))+
  #geom_contour(aes(z=exp(Z)), bins = 12, color = "gray30", linewidth = 0.1, alpha = 0.5) +
  geom_point(aes(x=XX[idx,1], y=XX[idx,2]), size=1, color="black")+
  coord_equal() +
  scale_x_continuous(expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_gradient2(midpoint=exp(5), low="blue", mid="green", high="red") +
  xlab("X1")+
  ylab("X2")+
  theme_bw()



pdf("fig_07_404.pdf", width=5, height=4)
print(p5)
dev.off()

pdf("fig_07_405.pdf", width=5, height=4)
print(p6)
dev.off()



