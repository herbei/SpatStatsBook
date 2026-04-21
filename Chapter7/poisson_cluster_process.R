rm(list = ls())
setwd('~/Dropbox/0Work/0Teach/6530/Notes/Figures/Ch7/')




N=20
XX=matrix(runif(2*N, 0,1), ncol=2)

plot1 = ggplot()+
  geom_point(aes(x= XX[,1], y=XX[,2]), pch=16)+
  theme_bw()+
  xlab("X1")+
  ylab("X2")

pdf("fig_07_300.pdf", width=4.5, height=4.5)
print(plot1)
dev.off()


#Simulate number of offspring
N1=10
p1=0.5
KK=rbinom(dim(XX)[1], N1, p1)

DD=data.frame(XX[,1], XX[,2], KK)
colnames(DD)=c("x", "y", 'Offspring')
plot2 = ggplot(DD, aes(x=x, y=y,label=Offspring))+
  geom_point()+
  geom_text(hjust=-0.5, vjust=0.1, color="red")+
  theme_bw()+
  xlab("X1")+
  ylab("X2")

pdf("fig_07_301.pdf", width=4.5, height=4.5)
print(plot2)
dev.off()


XXX=numeric(0)
R1 = 0.1

for(i in 1:dim(XX)[1]){
  
  if(KK[i]>0){
    for(j in 1:KK[i]){
      myR = R1 * runif(1,0,1); myA = 2*pi*runif(1,0,1)
      pp1 = XX[i,1] + myR*cos(myA)
      pp2 = XX[i,2] + myR*sin(myA)
      XXX = rbind(XXX, t(c(pp1,pp2)))
    }
  }
  
}


plot3 = ggplot()+
  geom_point(aes(x= XX[,1], y=XX[,2]), pch=16, color="black", size=2)+
  geom_point(aes(x= XX[,1], y=XX[,2]), shape=21, color="black", size=25, stroke=0.1)+
  geom_point(aes(x=XXX[,1], y=XXX[,2]), pch=16, size=1, color="red")+
#  geom_point(aes(x=0.3, y=0.5), pch=16)+
#  geom_point(aes(x=0.3, y=0.7), pch=16)+
#  geom_point(aes(x=0.3, y=0.6), shape=21, size=25)+
  theme_bw()+
  xlab("X1")+
  ylab("X2")
  
  pdf("fig_07_302.pdf", width=4.5, height=4.5)
  print(plot3)
  dev.off()
  

  
plot4 = ggplot()+
    #geom_point(aes(x= XX[,1], y=XX[,2]), pch=16, color="black", size=2)+
    #geom_point(aes(x= XX[,1], y=XX[,2]), shape=21, color="black", size=25, stroke=0.1)+
    geom_point(aes(x=XXX[,1], y=XXX[,2]), pch=16, size=1, color="red")+
    #  geom_point(aes(x=0.3, y=0.5), pch=16)+
    #  geom_point(aes(x=0.3, y=0.7), pch=16)+
    #  geom_point(aes(x=0.3, y=0.6), shape=21, size=25)+
    theme_bw()+
    xlab("X1")+
    ylab("X2")
  
  pdf("fig_07_303.pdf", width=4.5, height=4.5)
  print(plot4)
  dev.off()
  
  