
library(mvtnorm)
library(coda)
library(ggpubr)
library(statip)


korattikara<-function(z, Time, e, m,theta_ini){
  
  z<-as.matrix(z)
  N=dim(z)[1]
  d=length(theta_ini)
  theta=matrix(NA,d,Time+1)
  theta[,1]=theta_ini
  number_llik_eval=0
  ptm <- proc.time()
  
  for (t in (1:Time)) {
    tp=rnorm(d,theta[,t],0.06)
    u=runif(1)
    lbar=0
    lsqbar=0
    n=0
    done=FALSE
    mu0=(1/N)*(log(u)+log_prior(theta[,t])-log_prior(tp))
    batch=0
    index=(1:N)
    
    while(done==FALSE) {
      draw=sample(length(index),min(m,N-n))
      batch=base::c(index[draw],batch)
      n=n+min(m,N-n)
      index=index[-draw]
      
      
      l=llik(batch,z,tp)-llik(batch,z,theta[,t])
      lbar=mean(l)
      lsqbar=mean(l^2)
      sd_batch=sqrt((n/(n-1))*(lsqbar-lbar^2))
      sd_hat=sd_batch/sqrt(n)*sqrt(1-(n-1)/(N-1))
      
      delta=1-pt(abs((lbar-mu0)/sd_hat), n-1)
      
      if(delta<e) {
        if(lbar>mu0) { 
          theta[,t+1]=tp
        }
        else {
          theta[,t+1]=theta[,t]
        }
        done=TRUE 
        number_llik_eval=number_llik_eval+n
      }
      
    }
  }
  print(paste("Completed in ", proc.time()[3]-ptm[3]," sec", sep=""))
  return(list(theta=theta,number_llik_eval=number_llik_eval))
}

