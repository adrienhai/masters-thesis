#The posterior must be very peaked : many data points. Comparer les valeurs 
# de beta1 posterieur avec y[(1:100)] et y[(1:100000)]
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

global <- new.env(parent=emptyenv()) 
global$sizeOfBatch<-0

DGP3_korattikara<-korattikara(z,2e4,0.05,700,rep(0,d))
DGP3_korattikara$theta=DGP3_korattikara$theta[(1:d),-(1:1000)]

sum(global$sizeOfBatch)/n

ESS_korattikara<-effectiveSize(as.mcmc(DGP3_korattikara$theta[2,]))/(DGP3_korattikara$number_llik_eval/n)

1 - rejectionRate(as.mcmc(DGP3_korattikara$theta[1,]))
plot(DGP3_korattikara$theta[2,],type="l")

p1<-ggplot(data=as.data.frame(cbind(DGP3_MH_thi[,2],DGP3_MH_thi[,1])))+
  stat_density2d(geom="contour", bins=4, size= 1, aes(DGP3_MH_thi[,1],DGP3_MH_thi[,2],color="MH"))+
  stat_density2d(data=as.data.frame(cbind(DGP3_korattikara$theta[4,],DGP3_korattikara$theta[2,])), 
                 aes(DGP3_korattikara$theta[4,],DGP3_korattikara$theta[2,],color="korattikara"),geom="contour", bins=4, size= 1)+
  ylab(expression(theta[2]))+
  xlab(expression(theta[4]))+
  scale_colour_manual(name="",
                      values=c("korattikara"="dodgerblue4",
                               "MH"="red3"),
                      labels=c("MH","ApMHT"))



p2<-ggplot(data=as.data.frame(DGP3_korattikara$theta[2,]))+
  geom_line(aes(y=DGP3_korattikara$theta[2,], x=(1:18001),color="ApMHT"))+
  xlab("Time")+
  ylab(expression(theta))+
  scale_colour_manual(name="",
                      values=c("ApMHT"="dodgerblue4"))

ggarrange(p1, p2, nrow=2, common.legend = TRUE, legend="bottom")

hellinger(DGP3_korattikara$theta[2,],DGP3_MH[2,])
##########Efficiency korattikara#############
efficiency_korattikara=rep(0,length(seq(1000,100000,2000)))
i=1
for(n in seq(1000,100000,2000)) {
  efficiency_korattikara[i]=korattikara(z[(1:n),],1000,0.05,1000,theta)$n/n
  i=i+1
}
i=51
for(n in seq(1e5,1e6,1e5)) {
  efficiency_korattikara[i]=korattikara(z[(1:n),],1000,0.05,1000,theta)$n/n
  i=i+1
}

ggplot(as.data.frame(cbind(efficiency_korattikara[-c(1,2,3,4,5)],c(seq(11e3,1e5,2000),seq(1e5,1e6,1e5)))),
       aes(x=c(seq(11e3,1e5,2000),seq(1e5,1e6,1e5)), y=efficiency_korattikara[-c(1,2,3,4,5)]))+
  geom_smooth(se=FALSE, color="dodgerblue4", method=loess)+
  xlab("Data set size N")+
  ylab(expression(mean[n]))+
  geom_hline(yintercept = 0.02, color="coral")+
  scale_y_continuous(breaks = c(0.02, 0.1, 0.2, 0.3))

########Getting it right###########
hist(rnorm(10000,0,0.3))
DGP3_GITkorattikara=numeric(1000)

for (i in (0:300)) {
  phi=rnorm(1,mean=0,sd=0.3)
  y_sim=rnorm(n,mean = abs(phi),sigma)
  DGP3_GITkorattikara[i]=korattikara(y_sim,5000,0.05,1000,0)$theta[1,sample((500:5000),1)]
  print(i)
}
hist(DGP3_GITkorattikara)
ks.test(DGP3_GITkorattikara,pnorm,0,1) 

prior_sim=rnorm(1e5)
ggplot(data=as.data.frame(prior_sim))+
  geom_density(aes(prior_sim, color="Prior"))+
  geom_density(data=as.data.frame(DGP3_GITkorattikara),
               aes(DGP3_GITkorattikara,color="Marginal"))+
  scale_colour_manual(name="",
                      values=c("Prior"="red3",
                               "Marginal"="dodgerblue4"),
                      labels=c(expression(paste("Marginal of", " ", theta)),"Prior"))+
  theme(legend.position="bottom")+
  xlab(expression(theta))+
  ylab("Density")


#######Hellinger#########

statip::hellinger(DGP3_korattikara$theta[2,],DGP3_MH[2,])
 