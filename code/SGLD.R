library(rootSolve)
library(ggplot2)
library(ggpubr)
#TRAVAILLER SUR FONCTION SAMPLE POUR QU'ELLE SOIT PLUS EFFICACE SUR GRAND echantillon


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
  #return(list(b0=theta[1,], b1=theta[2,]))
}

global <- new.env(parent=emptyenv()) 
global$sizeOfBatch<-0


DGP3_SGLD<-SGLD(z, eps=0.022,Time=1,m=500,rep(0,d))
sum(global$sizeOfBatch)/n
DGP3_SGLD=DGP3_SGLD[(1:d),-(1:1000)]

ESS_SGLD<-effectiveSize(as.mcmc(DGP3_SGLD[2,]))/(2*20000*500/n)
plot(DGP3_SGLD[1,],type="l")
plot(DGP3_SGLD[2,],type="l")

ggplot(data=as.data.frame(cbind(DGP3_SGLD[1,])))+
  geom_density(aes(DGP3_SGLD[1,],color="SGLD"))
  
p1<-ggplot(data=as.data.frame(cbind(DGP3_MH[(1000:1.5e6),2],DGP3_MH[(1000:1.5e6),1])))+
  stat_density2d(geom="contour", bins=4, size= 1, aes(DGP3_MH[(1000:1.5e6),1],DGP3_MH[(1000:1.5e6),2],color="MH"))+
  stat_density2d(data=as.data.frame(cbind(DGP3_SGLD[4,],DGP3_SGLD[2,])), aes(DGP3_SGLD[4,],DGP3_SGLD[2,],color="SGLD"),geom="contour", bins=4, size= 1)+
  ylab(expression(theta[2]))+
  xlab(expression(theta[4]))+
  scale_colour_manual(name="",
                      values=c("SGLD"="dodgerblue4",
                               "MH"="red3"),
                      labels=c("MH","SGLD"))
p1<-ggplot(data=as.data.frame(cbind(DGP3_MH_thi[,2],DGP3_MH_thi[,1])))+
  stat_density2d(geom="contour", bins=4, size= 1, aes(DGP3_MH_thi[,1],DGP3_MH_thi[,2],color="MH"))+
  stat_density2d(data=as.data.frame(cbind(DGP3_SGLD[4,],DGP3_SGLD[2,])), aes(DGP3_SGLD[4,],DGP3_SGLD[2,],color="SGLD"),geom="contour", bins=4, size= 1)+
  ylab(expression(theta[2]))+
  xlab(expression(theta[4]))+
  scale_colour_manual(name="",
                      values=c("SGLD"="dodgerblue4",
                               "MH"="red3"),
                      labels=c("MH","SGLD"))

p2<-ggplot(data=as.data.frame(DGP3_SGLD[1,]))+
  geom_line(aes(y=DGP3_SGLD[2,], x=(1:19001),color="SGLD"))+
  xlab("Time")+
  ylab(expression(theta[2]))+
  scale_colour_manual(name="",
                      values=c("SGLD"="dodgerblue4"))
ggarrange(p1, p2, nrow=2, common.legend = TRUE, legend="bottom")

hellinger(DGP3_SGLD[2,],DGP3_MH[2,])


##Getting it right
hist(rnorm(1000,0,0.3))
DGP3_GITSGLD=numeric(100)

for (i in (20:100)) {
  phi=rnorm(d,mean=0,sd=1)
  pr=1/(1+exp(-x%*%phi))
  y_sim=rbinom(n,1,pr)
  DGP3_GITSGLD[i]=SGLD(y_sim, eps=0.022,Time=2000,m=500,rep(0,d))[2,sample((1000:2000),1)]
  print(i)
}
hist(DGP3_GITSGLD)
ks.test(DGP3_GITSGLD,pnorm,0,1) #0.8906

prior_sim=rnorm(1e5)
ggplot(data=as.data.frame(prior_sim))+
  geom_density(aes(prior_sim, color="Prior"))+
  geom_density(data=as.data.frame(DGP3_GITSGLD),
               aes(DGP3_GITSGLD,color="Marginal"))+
  scale_colour_manual(name="",
                      values=c("Prior"="red3",
                               "Marginal"="dodgerblue4"),
                      labels=c(expression(paste("Marginal of", " ", theta)),"Prior"))+
  theme(legend.position="bottom")+
  xlab(expression(theta))+
  ylab("Density")

statip::hellinger(DGP3_SGLD[2,],DGP3_MH[2,])

seq(1,2e6,50)
DGP3_MH_thi=DGP3_MH[seq(1,1.9e6,125),]
