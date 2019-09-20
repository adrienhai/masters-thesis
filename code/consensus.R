library(coda)
library(ggplot2)
library(ggpubr)
library(forecast)

consensus_batch<-function(z,s,Time,theta_ini){
  z<-as.matrix(z)
  N=dim(z)[1]
  shuffle=sample((1:N),N,replace=FALSE)
  d=length(theta_ini)
  theta=list()
  S=list()
  W=list()
  
  ptm <- proc.time()
  
  for (i in (1:s)) {
    S[[i]]=shuffle[(((i-1)*(N-N%%s)/s+1):(i*(N-N%%s)/s))]#set of indices
    theta[[i]]=matrix(NA,nrow=d,ncol=Time+1)
    theta[[i]][,1]=theta_ini
    llk_current=sum(S[[i]],z,theta_ini)
    for (t in (1:Time)) {
      tp=rnorm(d,theta[[i]][,t],0.1)
      llk_proposal=sum(llik(S[[i]],z,tp))
      ratio=(llk_proposal+(1/s)*log_prior(tp))-
        (llk_current+(1/s)*log_prior(theta[[i]][,t]))
      if(log(runif(1))<=ratio)  {
        theta[[i]][,t+1]=tp
        llk_current=llk_proposal
      }
      else
        theta[[i]][,t+1]=theta[[i]][,t]
    }
    
    W[[i]]=solve(var(t(theta[[i]])))
  }
  print(paste("Each parallel computing was completed in ", (proc.time()[3]-ptm[3])/s," sec", sep=""))
  return(list(theta=theta,W=W))
}

consensus<-function(theta,W){
  ptm<-proc.time()
  Time=dim(theta[[1]])[2]
  d=dim(theta[[1]])[1]
  s=length(theta)
  
  norm_const=Reduce('+',W)
  norm_const=solve(norm_const)
  out=0
  for (i in (1:s)) {
    out=out+W[[i]]%*%theta[[i]]
  }
  out=norm_const%*%out
  print(paste("Combining part completed in ", proc.time()[3]-ptm[3]," sec", sep=""))
  return(out)
}

