


lam = 20
N=rpois(1,lam)
Z=matrix(runif(2*N,0,1), nrow=N, ncol=2)


p = 0.2
Z1 = numeric(0)
Z2 = numeric(0)
for(i in 1:N){
      if (runif(1,0,1)<p){
            Z1 = rbind(Z1,Z[i,])
      }
      else{
            Z2 = rbind(Z2, Z[i,])
      }
}

pdf(file="fig_07_026.pdf", width=7, height=3.5)
par(mfrow=c(1,3), cex=0.8, mar=c(3,3,2,1), mgp=c(2,0.5,0), bty="L", pty="s")
plot(Z[,1], Z[,2], xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pch=16, cex=1, main="lam = 20")
plot(Z1[,1], Z1[,2], xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pch=17, cex=1, main="lam = 0.2 * 20")
plot(Z2[,1], Z2[,2], xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), pch=18, cex=1, main="lam = 0.8 * 20")
dev.off()


