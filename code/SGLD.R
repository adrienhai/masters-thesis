library(rootSolve)
library(ggplot2)
library(ggpubr)


SGLD<-function(z, eps=0.0001, Time=1000, m=1000, theta_ini=0) {
  
  z<-as.matrix(z)
  N=dim(z)[1]
  d=length(theta_ini)
  theta=matrix(NA,nrow=d, ncol=Time+1)
  theta[,1]=theta_ini
  
  for(i in (2:(Time+1))) {
    batch=sample((1:N),m)
    theta[,i]=theta[,i-1] + (eps/2)*(gradient(log_prior,theta[,i-1])
                                     +(N/m)*apply(gradient(llik,x=theta[,i-1],i=batch,z=z),2,sum)
                                     +rnorm(d,mean=0,eps))
  }
  
  return(theta)
}

