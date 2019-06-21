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

Time=20000
s=3
DGP3_consensus=consensus_batch(z,s,Time,theta)
plot(DGP3_consensus$theta[[1]][2,],type="l")
1 - rejectionRate(as.mcmc(DGP3_consensus$theta[[1]][2,]))
ggplot(data=as.data.frame(DGP3_consensus$theta[[2]][2,]), aes(DGP3_consensus$theta[[2]][2,]))+
  geom_density()

#Remove burn-in
DGP3_consensus$theta[[1]][2,]=DGP3_consensus$theta[[1]][2,-(1:1000)]
DGP3_consensus$theta[[1]][4,]=DGP3_consensus$theta[[1]][4,-(1:1000)]

DGP3_consensus=consensus(DGP3_consensus$theta,DGP3_consensus$W)

#Cost for each machine : 2T((n-n%%s)/s)/n

ESS_consensus=effectiveSize(as.mcmc(DGP3_consensus[2,]))/(Time*((n-n%%s)/s)/n)
1 - rejectionRate(as.mcmc(consensus_batch(z,s,Time,theta)$theta[[1]][1,]))

p1<-ggplot(data=as.data.frame(cbind(DGP3_MH_thi[,2],DGP3_MH_thi[,1])))+
  stat_density2d(geom="contour", bins=4, size= 1, aes(DGP3_MH_thi[,1],DGP3_MH_thi[,2],color="MH"))+
  stat_density2d(data=as.data.frame(cbind(DGP3_consensus[4,],DGP3_consensus[2,])), aes(DGP3_consensus[4,],DGP3_consensus[2,],color="consensus"),geom="contour", bins=4, size= 1)+
  ylab(expression(theta[2]))+
  xlab(expression(theta[4]))+
  scale_colour_manual(name="",
                      values=c("consensus"="dodgerblue4",
                               "MH"="red3"),
                      labels=c("MH","Consensus"))



p2<-ggplot(data=as.data.frame(DGP3_consensus[1,]))+
  geom_line(aes(y=DGP3_consensus[2,], x=(1:20001),color="Consensus"))+
  xlab("Time")+
  ylab(expression(theta[2]))+
  scale_colour_manual(name="",
                      values=c("Consensus"="dodgerblue4"))
 ggarrange(p1, p2, nrow=2, common.legend = TRUE, legend="bottom")

statip::hellinger(DGP3_consensus[2,],DGP3_MH[2,])
ggarrange(p1, p2, nrow=2, common.legend = TRUE, legend="bottom")

##Getting it right
hist(rnorm(1000,0,0.3))
DGP3_GITconsensus=numeric(100)

for (i in (1:100)) {
  phi=rnorm(d,mean=0,sd=1)
  pr=1/(1+exp(-x%*%phi))
  y_sim=rbinom(n,1,pr)
  temp=consensus_batch(y_sim,3,3000,theta)
  DGP3_GITconsensus[i]=consensus(temp$theta,temp$W)[2,sample((500:3000),1)]
  print(i)
}
hist(DGP3_GITconsensus)
ks.test(DGP3_GITconsensus,pnorm,0,1)
prior_sim=rnorm(1e5)
ggplot(data=as.data.frame(prior_sim))+
  geom_density(aes(prior_sim, color="Prior"))+
  geom_density(data=as.data.frame(DGP3_GITconsensus),
               aes(DGP3_GITconsensus,color="Marginal"))+
  scale_colour_manual(name="",
                      values=c("Prior"="red3",
                               "Marginal"="dodgerblue4"),
                      labels=c(expression(paste("Marginal of", " ", theta)),"Prior"))+
  theme(legend.position="bottom")+
  xlab(expression(theta))+
  ylab("Density")

N=seq(2,6,1)
N=10^N
consensus_std_llk_eval=(2*Time*((N-N%%s)/s))

plot(y=consensus_std_llk_eval/T,N)
