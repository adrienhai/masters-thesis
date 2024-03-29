library(plyr)
library(pracma)
library(ggplot2)
library(ggpubr)
library(coda)
dist_eucl<-function(x1,x2) {
  out=0
  if (is.vector(x1)) {
    out=sum((x1-x2)^2)
  }
  
  else {
    out=sweep(x1,2,x2)
    out=apply(out^2,1,sum)
  }

  return(sqrt(out))
}

clustering <- function(z,eps) {
  # y<-scale(y) SCALE DATA SET
  # x<-scale(x)
  r<-as.matrix(z)
  
  n=dim(r)[1]
  I=rep(0,n) #Indicator variable, 1 if observation i has been given a cluster
  z=0 #centroids
  C=0 #indices inside each clusters
  k=0
  index=seq(1,n,1)
  
  for (j in (1:n)) {
    if(I[j]==0) {
      
      
      p<-which(dist_eucl(r[as.logical(1-I),],r[j,])<=eps, arr.ind = TRUE)
      C_current=index[p]
      index=index[-p]
      
      if (length(C_current)==1)
        z=cbind(r[C_current,],z)
      else
        z=cbind((1/length(C_current))*apply(as.matrix(r[C_current,]),2,sum),z)
      I[C_current]=1
      C=rbind.fill(as.data.frame(t(C_current)),as.data.frame(C))
      k=k+1
    }
  }
  
  return(list(K=k,z=z,C=t(C)))
}

data=cbind(x,z)
data[which(data[,d+1]==1),d+1]=100
data[which(data[,d+1]==0),d+1]=-100
out<-clustering(data,5.3)
for(k in (1:out$K)) {
  out$z[21,k]=ifelse (out$z[21,k]==100,1,0)
}
data[which(data[,d+1]==100),d+1]=1
data[which(data[,d+1]==-100),d+1]=0

#Computing second constant B see appendix of paper
B=list()
for (k in (1:(out$K))) {
  B[[k]]=matrix(0,dim(data)[2],dim(data)[2])
  for (i in na.omit(out$C[,k])) {
    b=data[i,]-out$z[,k]
    B[[k]]=B[[k]]+b%*%t(b)
  }
}
quiroz<-function(data,z,K,C,m, Time, theta_ini) {
  
  N=dim(data)[1]
  d=length(theta_ini)
  
  d_b=rep(0,m)
  d_bp=rep(0,m)
  theta=matrix(NA,nrow=length(theta_ini), ncol=Time+1)
  theta[,1]=theta_ini
  
  for(t in (1:Time)) {
    
    bp=rnorm(d,theta[,t],0.02)
    u=sample((1:N),m,replace=TRUE)
    for (i in (1:m)) {
      coord=which(out$C==u[i],arr.ind = TRUE)[2]

      d_b[i]=llik(data[u[i],],theta[,t])-llik(z[,coord],theta[,t])-
        gradient_llik(z=z[,coord],theta=theta[,t]) %*% (data[u[i],]-z[,coord])-
        0.5*t(data[u[i],]-z[,coord])%*%hessian_llik(z=z[,coord],theta=theta[,t])%*%(data[u[i],]-z[,coord])
      d_bp[i]=llik(data[u[i],],bp)-llik(z[,coord],bp)-
        gradient_llik(z=z[,coord],theta=bp)%*%(data[u[i],]-z[,coord])-
        0.5*t(data[u[i],]-z[,coord])%*%hessian_llik(z=z[,coord],theta=bp)%*%(data[u[i],]-z[,coord])
      
    }
    mu_b=mean(d_b)
    sigma_b=var(d_b)
    
    mu_bp=mean(d_bp)
    sigma_bp=var(d_bp)
    
    first_term_b=0
    first_term_bp=0
    third_term_b=0
    third_term_bp=0
    for (k in (1:K)) {
      Nk=length(na.omit(out$C[,k]))
      first_term_b=first_term_b+Nk*llik(z[,k],theta[,t])
      first_term_bp=first_term_bp+Nk*llik(z[,k],bp)
      
      H_bp=hessian_llik(z=z[,k],theta=bp)
      H_b=hessian_llik(z=z[,k],theta=theta[,t])
      third_term_b=third_term_b+H_b*B[[k]]
      third_term_bp=third_term_bp+H_bp*B[[k]]
    }
    q_b=first_term_b+0.5*sum(rowSums(third_term_b))
    q_bp=first_term_bp+0.5*sum(rowSums(third_term_bp))
    
    l_hat_b=q_b+N*mu_b-(N^2)/(2*m)*sigma_b
    l_hat_bp=q_bp+N*mu_bp-(N^2)/(2*m)*sigma_bp
    
    prob=(l_hat_bp+log_prior(bp))-(l_hat_b+log_prior(theta[,t]))
    if(log(runif(1))<=prob) {
      theta[,t+1]=bp
    }
    else {
      theta[,t+1]=theta[,t]
    }
    print(t)
  }
  return(theta)
}

