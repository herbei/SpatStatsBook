rm(list = ls())
library(png)



find.neighbors<-function(i,N){
  ## Assume an N x N matrix; i = 1,2, ... , N^2
  ## function returns the N-S-E-W neighbors of i
  
  if(i%%N >0){
    my.row = (i %/% N) + 1
    my.col = i %% N
  }else{
    my.row = i%/%N
    my.col = N
  }  
  out=numeric(0)
  if(my.row>1){out=c(out, (my.row-2)*N + my.col)}
  if(my.row<N){out=c(out, my.row*N + my.col)}
  if(my.col>1){out=c(out, i-1)}
  if(my.col<N){out=c(out,i+1)}
  return(out)
}



YY=readPNG('wolf_noisy2.png')
Y = as.vector(t(YY)) # convert XX to a vector (by row)
  
  
  
sig2=0.01
rho=0.8
N=dim(YY)[1]

# initialize X
XX=YY
X = as.vector(t(XX)) # convert XX to a vector (by row)
kappa=3 # starting value



# compute the number of neighbors for each pixel
NN = matrix(4, nrow=N, ncol=N)
NN[1,] = NN[1,]-1
NN[N,]=NN[N,]-1
NN[,1]=NN[,1]-1
NN[,N]=NN[,N]-1
D = as.vector(NN)

MCS=30000 # number of MCMC iterations
MCO=50 # thin the output by this much

# some variables for storing MCMC output
sample_count=1
kappa_keep=numeric(0)
XXmean = XX

#begin MCMC loop
for(k in 1:MCS){
  
  # sample kappa given X
  a = (N^2)/2+1
  b=0.1 + (rho/2)* (sum(sum(diff(XX)^2)) + sum(sum(diff(t(XX))^2))) + ((1-rho)/2)*sum(D*(X^2))
  kappa = rgamma(1, shape=a, rate=b)
  
  
  # sample X[i] given X[-i], kappa
  for(i in 1:(N^2)){
    Ni = D[i]
    ci = kappa*Ni + 1/sig2
    idx = find.neighbors(i,N)
    di = Y[i]/sig2 + kappa*rho*sum(X[idx])
    X[i] = rnorm(1,mean=di/ci, sd=sqrt(1/ci))
  }

  
  XX = matrix(X, nrow=N, ncol=N, byrow=TRUE)

  
  if(k %% MCO == 0){
    print(k)
    sample_count=sample_count+1
    kappa_keep = c(kappa_keep, kappa)
    XXmean = XXmean + XX
  }
  
}





post.mean  = XXmean/sample_count
writePNG(post.mean, 'postmean.png')


pdf(file="figure.pdf", width=5, height=3.5)
par(mfrow=c(1,1), cex=0.8, mar=c(3.5,3.5,0,0), mgp=c(2,0.5,0), bty="L")
plot(kappa_keep, ylim=c(5.8,6.3), main="", xlab="iteration", ylab="kappa", type="l" )
abline(h=mean(kappa_keep), col="red", lwd=2)
dev.off()

pdf(file="another_figure.pdf", width=4, height=3.5)
par(mfrow=c(1,1), cex=0.8, mar=c(3.5,0,0,0), mgp=c(2,0.5,0), bty="L")
hist(kappa_keep, xlim=c(5.8, 6.2), xlab="kappa", main="")
dev.off()

