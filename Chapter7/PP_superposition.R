


lam1 = 3
N1=rpois(1,lam1)
Z1=matrix(runif(2*N1,0,1), nrow=N1, ncol=2)

lam2 = 10
N2=rpois(1,lam2)
Z2=matrix(runif(2*N2,0,1), nrow=N2, ncol=2)


Z = rbind(Z1,Z2)

pdf(file="fig_07_025.pdf", width=7, height=3.5)
par(mfrow=c(1,3), cex=0.8, mar=c(3,3,2,1), mgp=c(2,0.5,0), bty="L", pty="s")
plot(Z1[,1], Z1[,2], xlab="", ylab="", xlim=c(0,1), pch=16, cex=1, main="lam=3")
plot(Z2[,1], Z2[,2], xlab="", ylab="", xlim=c(0,1), pch=17, cex=1, main="lam=10")
plot(Z[,1], Z[,2], xlab="", ylab="", xlim=c(0,1), pch=18, cex=1, main="lam=3+10=13")
dev.off()


