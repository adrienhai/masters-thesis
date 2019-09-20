library(mvtnorm)
library(ggplot2)
library(coda)
library(ggpubr)

batch_MCMC<-function(z,s,Time,theta_ini) {
  z<-as.matrix(z)
  N=dim(z)[1]
  
  shuffle=sample((1:N),N,replace=FALSE)
  d=length(theta_ini)
  theta=list()
  S=list()
  
  ptm <- proc.time()
  for (i in (1:s)) {
    S[[i]]=shuffle[(((i-1)*(N-N%%s)/s+1):(i*(N-N%%s)/s))]#set of indices
    theta[[i]]=matrix(NA,nrow=d,ncol=Time+1)
    theta[[i]][,1]=theta_ini
    for (t in (1:Time)) {
      tp=rnorm(d,theta[[i]][,t],0.13)
      ratio=(sum(llik(S[[i]],z,tp))+(1/s)*log_prior(tp))-
        (sum(llik(S[[i]],z,theta[[i]][,t]))+(1/s)*log_prior(theta[[i]][,t]))
      if(log(runif(1))<=ratio) 
        theta[[i]][,t+1]=tp
      else
        theta[[i]][,t+1]=theta[[i]][,t]
    }
    
  }
  print(paste("Each parallel computing was completed in ", (proc.time()[3]-ptm[3])/s," sec", sep=""))
  return(theta)
}

gaussian_kernel<-function(theta) {
  ptm<-proc.time()
  Time=dim(as.matrix(theta[[1]]))[2]
  d=dim(as.matrix(theta[[1]]))[1]
  s=length(theta)
  t=sample((1:Time),s,replace = TRUE)
  DGP3_neiswanger=matrix(NA,nrow=d,ncol=Time)
  for(i in (1:Time)) {
    h=1.7e-1
    for(k in (1:s)) {
      c<-t
      c[k]=sample((1:Time),1,replace = TRUE)
      theta_bar_t=0
      theta_bar_c=0
      w_t=0
      w_c=0
      for (j in (1:s)) {
        theta_bar_t=theta_bar_t+theta[[j]][,t[j]]
        theta_bar_c=theta_bar_c+theta[[j]][,c[j]]
      }
      theta_bar_t=(1/s)*theta_bar_t
      theta_bar_c=(1/s)*theta_bar_c
      
      for (j in (1:s)) {
        w_t=w_t + dmvnorm(theta[[j]][,t[j]],mean=theta_bar_t,sigma=h^2*diag(d), log=TRUE)
        w_c=w_c + dmvnorm(theta[[j]][,c[j]],mean=theta_bar_c,sigma=h^2*diag(d), log=TRUE)
      }
      
      if(log(runif(1))<w_c-w_t) {
        t<-c
        theta_bar_t<-theta_bar_c
      }
      
    }
    DGP3_neiswanger[,i]=rmvnorm(1,mean=theta_bar_t,sigma=(h^2/s)*diag(d))
    
  }
  print(paste("Combining part completed in ", proc.time()[3]-ptm[3]," sec", sep=""))
  
  return(DGP3_neiswanger)
}

