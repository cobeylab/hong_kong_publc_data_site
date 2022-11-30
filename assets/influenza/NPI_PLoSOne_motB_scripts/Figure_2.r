#
# R syntax to reproduce information for Figure 2 (HK part) from:
#
# Cowling BJ, Ip DKM, Fang VJ, et al.
# Modes of transmission of influenza B virus in households
# PLoS ONE, 2014 (in press).
#
# Last updated by ang VJ and Cowling BJ.
# August 20, 2014
#


#### mcmc functionlogprior = function(current)
logprior = function(current)
{
   current$logprior=dunif(current$theta.a,0,1,log=TRUE)+dunif(current$theta.c,0,1,log=TRUE)+
                    dunif(current$pi1,0,1,log=TRUE)+dunif(current$pi2,0,1,log=TRUE)+dunif(current$pi3,0,1,log=TRUE)+
                    dunif(current$lambda1,0,1,log=TRUE)+dunif(current$lambda2,0,1,log=TRUE)+dunif(current$lambda3,0,1,log=TRUE)+
                    dunif(current$phi,0,10,log=TRUE)
   current
}

# the likelihood function
loglikelihood = function(current,data)
{
     hu.fun <- function(T){
       hazard <- list(NA)
       hazard[[1]] <- current$phi*current$lambda1^current$phi*T^(current$phi-1)*(1-current$beta1)^data$X_hh    # contact
       hazard[[2]] <- current$phi*current$lambda2^current$phi*T^(current$phi-1)*(1-current$beta2)^data$X_fm    # droplet
       hazard[[3]] <- current$phi*current$lambda3^current$phi*T^(current$phi-1)    # aerosol
       hazard
     }
     Su.fun <- function(T){                                                                                                         
       Spart <- (current$lambda1*T)^current$phi*(1-current$beta1)^data$X_hh+
                (current$lambda2*T)^current$phi*(1-current$beta2)^data$X_fm+
                (current$lambda3*T)^current$phi
     }
     huT0 <- hu.fun(data$T0); huT1 <- hu.fun(data$T1); huT2 <- hu.fun(data$T2)
     hu <- (huT0[[1]]+data$mark1*huT1[[1]]+data$mark2*huT2[[1]])*current$pi1^data$ILI*(1-current$pi1)^(1-data$ILI)+
           (huT0[[2]]+data$mark1*huT1[[2]]+data$mark2*huT2[[2]])*current$pi2^data$ILI*(1-current$pi2)^(1-data$ILI)+
           (huT0[[3]]+data$mark1*huT1[[3]]+data$mark2*huT2[[3]])*current$pi3^data$ILI*(1-current$pi3)^(1-data$ILI)
     Su <- exp(-Su.fun(data$T0)-data$mark1*Su.fun(data$T1)-data$mark2*Su.fun(data$T2))
     St <- exp(-Su.fun(data$trunc_day))
     theta <- (1-data$child)*current$theta.a+data$child*current$theta.c

     likelihood <- ((1-theta)*Su*hu)^data$delta*(theta+(1-theta)*Su)^(1-data$delta)/(theta+(1-theta)*St)
     current$loglikelihood = log(likelihood)
     current
}

bad_config=function(current)#this checks no negative parameters, etc
{
  ok_config=TRUE

  if(current$pi1<0|current$pi1>1)ok_config=FALSE
  if(current$pi2<0|current$pi2>1)ok_config=FALSE
  if(current$pi3<0|current$pi3>1)ok_config=FALSE
  if(current$pi1>current$pi3)ok_config=FALSE
  if(current$pi2>current$pi3)ok_config=FALSE
  if(current$lambda1<0)ok_config=FALSE
  if(current$lambda2<0)ok_config=FALSE
  if(current$lambda3<0)ok_config=FALSE
  if(current$phi<0)ok_config=FALSE
  if(current$theta.a<0|current$theta.a>1)ok_config=FALSE
  if(current$theta.c<0|current$theta.c>1)ok_config=FALSE
  bad=FALSE;if(!ok_config)bad=TRUE
  bad
}

metropolis=function(old,current,data)
{
  REJECT=bad_config(current)
  if(!REJECT)
  {
    current=logprior(current)
    current=loglikelihood(current,data)
    lu=log(runif(1))
    if(lu>(current$logprior+sum(current$loglikelihood)-old$logprior-sum(old$loglikelihood)))REJECT=TRUE
  }
  if(REJECT)current=old
  current
}

mcmc=function(data,MCMC_iterations=1000,BURNIN_iterations=100,THINNING=1,r1,r2)
{
  current=list(beta1=r1, beta2=r2,
               phi=runif(1),lambda1=runif(1),lambda2=runif(1),lambda3=runif(1),pi1=0.1,pi2=0.1,pi3=0.9,theta.a=runif(1),theta.c=runif(1),
               logprior=0,loglikelihood=0)
  dump=list(phi=c(),lambda1=c(),lambda2=c(),lambda3=c(),pi1=c(),pi2=c(),pi3=c(),theta.a=c(),theta.c=c(),
            accept=matrix(NA,ncol=11,nrow=MCMC_iterations))
  sigma = list(phi=0.8,lambda1=0.15,lambda2=0.12,lambda3=0.1,pi1=0.3,pi2=0.25,pi3=0.25,theta.a=0.05,theta.c=0.15)  # hk flu B
  #sigma = list(phi=0.25,lambda1=0.15,lambda2=0.08,lambda3=0.08,pi1=0.3,pi2=0.3,pi3=0.3,theta.a=0.1,theta.c=0.2)  # bk flu B
  current=logprior(current)
  current=loglikelihood(current,data)
  for (iteration in (-BURNIN_iterations+1):MCMC_iterations)
  {
    old=current; current$phi=rnorm(1,current$phi,sigma$phi); current=metropolis(old,current,data)
    old=current; current$lambda1=rnorm(1,current$lambda1,sigma$lambda1); current=metropolis(old,current,data)
    old=current; current$lambda2=rnorm(1,current$lambda2,sigma$lambda2); current=metropolis(old,current,data)
    old=current; current$lambda3=rnorm(1,current$lambda3,sigma$lambda3); current=metropolis(old,current,data)
    old=current; current$pi1=rnorm(1,current$pi1,sigma$pi1); current=metropolis(old,current,data)
    old=current; current$pi2=rnorm(1,current$pi2,sigma$pi2); current=metropolis(old,current,data)
    old=current; current$pi3=rnorm(1,current$pi3,sigma$pi3); current=metropolis(old,current,data)
    old=current; current$theta.a=rnorm(1,current$theta.a,sigma$theta.a); current=metropolis(old,current,data)
    old=current; current$theta.c=rnorm(1,current$theta.c,sigma$theta.c); current=metropolis(old,current,data)

    if(iteration>0) #dump to file
    {
      rr <- iteration
      dump$phi[rr]=current$phi
      dump$lambda1[rr]=current$lambda1
      dump$lambda2[rr]=current$lambda2
      dump$lambda3[rr]=current$lambda3
      dump$pi1[rr]=current$pi1
      dump$pi2[rr]=current$pi2
      dump$pi3[rr]=current$pi3
      dump$theta.a[rr]=current$theta.a
      dump$theta.c[rr]=current$theta.c

      if(rr>1){
       dump$accept[rr,3] <- 1*(dump$phi[rr]!=dump$phi[rr-1])
       dump$accept[rr,4] <- 1*(dump$lambda1[rr]!=dump$lambda1[rr-1])
       dump$accept[rr,5] <- 1*(dump$lambda2[rr]!=dump$lambda2[rr-1])
       dump$accept[rr,6] <- 1*(dump$lambda3[rr]!=dump$lambda3[rr-1])
       dump$accept[rr,7] <- 1*(dump$pi1[rr]!=dump$pi1[rr-1])
       dump$accept[rr,8] <- 1*(dump$pi2[rr]!=dump$pi2[rr-1])
       dump$accept[rr,9] <- 1*(dump$pi3[rr]!=dump$pi3[rr-1])
       dump$accept[rr,10] <- 1*(dump$theta.a[rr]!=dump$theta.a[rr-1])
       dump$accept[rr,11] <- 1*(dump$theta.c[rr]!=dump$theta.c[rr-1])
      }
    }
  }
  dump
}


# MCMC
set.seed(12345)
iter=10000;burnin=iter/10
r1 <- c(0:9/10,0.99); r2 <- c(0:9/10,0.99); n <- length(r1)
fit <- list(NA)
aero.prop <- matrix(NA,ncol=n,nrow=n)
for(i in 1:n){
print(i)
 for(j in 1:n){
  fit[[(i-1)*n+j]] <- mcmc(hkdata,iter,burnin,THINNING=1,r1[i],r2[j])
 }
}

for(i in 1:n){
  for(j in 1:n){
      aero.prop[i,j] <- colMeans(prop.table(cbind(fit[[(i-1)*n+j]]$lambda1^fit[[(i-1)*n+j]]$phi,
                                                  fit[[(i-1)*n+j]]$lambda2^fit[[(i-1)*n+j]]$phi,
                                                  fit[[(i-1)*n+j]]$lambda3^fit[[(i-1)*n+j]]$phi),1))[3]
  }
}

##
##   PLOT
##

hkb <- as.data.frame(aero.prop)

hkb <- as.matrix(hkb);
r1 <- c(0:9/10,0.99); r2 <- c(0:9/10,0.99);

library("fields")

windows(width=4,height=4)
par(mar=c(4,5,2,0))
contour(predict.surface(Tps(as.matrix(expand.grid(x=r1,y=r2)),as.vector(hkb))),
        labcex=1,main="Hong Kong",axes=F,levels=1:10/10,
        labels=c("10%","20%","30%","40%","50%","60%","70%","80%","90%","100%"))
axis(1,at=0:5/5,labels=c("0%","20%","40%","60%","80%","100%"))
axis(2,at=0:5/5,labels=c("0%","20%","40%","60%","80%","100%"),las=1)
mtext("Efficacy of hand hygiene in
reducing contact transmission",side=1,line=3,cex=1)
mtext("Efficacy of surgical mask in
reducing droplet transmission",side=2,line=3,cex=1)

#
# End of script.
#
