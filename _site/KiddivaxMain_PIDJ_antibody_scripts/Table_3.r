#
# R syntax to reproduce information for Table 3 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Humoral antibody response after receipt of inactivated seasonal
# influenza vaccinations one year apart in children
# PIDJ, 2012.
#
# Last updated by Ng S, Fang VJ and Cowling BJ.
# February 15, 2013
#

require(Hmisc)
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_PIDJ_antibody_scripts/dataframe.r")

tab <- matrix(rep(NA,10*13), ncol=13,
              dimnames=list(c("Age 7-8","Age 9-15","Not 08-09 TIV","Received 08-09 TIV","Not infected","Infected","Before Oct 1","After Oct 1","Within 28 days","Beyond 28 days"),
                            c("n","sH1","CI_low","CI_up","sH3","CI_low","CI_up","B","CI_low","CI_up","pH1","CI_low","CI_up")))

##
dat <- sero[sero$intervention=="TIV",]

dat$prevax.sH1.c <-log2(dat$prevax.sH1/2.5)     #rescale titer
dat$prevax.sH3.c <-log2(dat$prevax.sH3/2.5)
dat$prevax.pH1.c <-log2(dat$prevax.pH1/2.5)
dat$prevax.B.Floride.c <-log2(dat$prevax.B.Floride/2.5)
dat$prevax.B.Brisbane.c <-log2(dat$prevax.B.Brisbane/2.5)
dat$postvax.sH1.c <-log2(dat$postvax.sH1/2.5)
dat$postvax.sH3.c <-log2(dat$postvax.sH3/2.5)
dat$postvax.pH1.c <-log2(dat$postvax.pH1/2.5)
dat$postvax.B.Floride.c <-log2(dat$postvax.B.Floride/2.5)
dat$postvax.B.Brisbane.c <-log2(dat$postvax.B.Brisbane/2.5)
dat$pilot.sh1.pre.c <-log2(dat$pilot.sh1.pre/2.5)
dat$pilot.sh3.pre.c <-log2(dat$pilot.sh3.pre/2.5)
dat$pilot.sh1.postv.c <-log2(dat$pilot.sh1.postv/2.5)
dat$pilot.sh3.postv.c <-log2(dat$pilot.sh3.postv/2.5)

dat$prevax_calday <- as.numeric(as.Date(dat$prevax, "%d/%m/%Y")- as.Date("29/8/2009", format= "%d/%m/%Y"))
dat$group <- paste(dat$pilot.intervention,dat$intervention,sep="")

# MI
set.seed(12347)
m <- 10
dat.i <- aregImpute( ~ I(group) + I(age) + I(male) + I(aftervactime) + I(prevax_calday)
                      + I(prevax.sH1.c) + I(prevax.sH3.c) + I(prevax.pH1.c) + I(prevax.B.Floride.c) + I(prevax.B.Brisbane.c)
                      + I(postvax.sH1.c) + I(postvax.sH3.c) + I(postvax.pH1.c) + I(postvax.B.Floride.c) + I(postvax.B.Brisbane.c)
                      + I(pilot.ph1s) + I(pilot.sh1) + I(pilot.sh3) + I(pilot.ili) + I(pilot.ari)
                      + I(pilot.sh1.pre.c) + I(pilot.sh3.pre.c) + I(pilot.sh1.postv.c) +  I(pilot.sh3.postv.c), data=dat, n.impute=m)

dat.nomiss <- list(dat, dat, dat, dat, dat, dat, dat, dat, dat, dat)

for(i in 1:10){
    dat.nomiss[[i]]$postvax.sH1[is.na(dat.nomiss[[i]]$postvax.sH1)] <- 2.5 * (2^dat.i$imputed$postvax.sH1.c[,i])
    dat.nomiss[[i]]$postvax.sH3[is.na(dat.nomiss[[i]]$postvax.sH3)] <- 2.5 * (2^dat.i$imputed$postvax.sH3.c[,i])
    dat.nomiss[[i]]$postvax.pH1[is.na(dat.nomiss[[i]]$postvax.pH1)] <- 2.5 * (2^dat.i$imputed$postvax.pH1.c[,i])
    dat.nomiss[[i]]$postvax.B.Brisbane[is.na(dat.nomiss[[i]]$postvax.B.Brisbane)] <- 2.5 * (2^dat.i$imputed$postvax.B.Brisbane.c[,i])
}

##
pre <- c("prevax.sH1","prevax.sH3","prevax.B.Brisbane","prevax.pH1" )
postv <- c("postvax.sH1","postvax.sH3","postvax.B.Brisbane","postvax.pH1")
infect <- c("pilot.sh1","pilot.sh3","pilot.B","pilot.ph1s")

## function for combining results
combine.mi <- function(titer1,titer2,k,n.impute) {

  model <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
  for (i in 1:n.impute){
    if(k!=3) model[[i]] <- glm(log(dat.nomiss[[i]][,titer2]/dat.nomiss[[i]][,titer1])~
                                cut(age,c(0,8,18))+factor(group)+dat.nomiss[[i]][,infect[k]]+cut(prevax_calday,c(0,32,100))+cut(aftervactime,c(0,28,100)),
                                data=dat.nomiss[[i]], family = gaussian)
    else model[[i]] <- glm(log(dat.nomiss[[i]][,titer2]/dat.nomiss[[i]][,titer1])~
                                cut(age,c(0,8,18))+factor(group)+cut(prevax_calday,c(0,32,100))+cut(aftervactime,c(0,28,100)),
                                data=dat.nomiss[[i]], family = gaussian)
  }

  Q <-matrix(c(model[[1]]$coef,model[[2]]$coef,model[[3]]$coef,model[[4]]$coef,model[[5]]$coef,     # Coefficients
               model[[6]]$coef,model[[7]]$coef,model[[8]]$coef,model[[9]]$coef,model[[10]]$coef), byrow=F, ncol=10)
  Qbar <-rowMeans(Q)  #overall estimate
  U <- matrix(c(summary(model[[1]])$coef[,2],summary(model[[2]])$coef[,2],summary(model[[3]])$coef[,2],   #standard error associated with Qj
                summary(model[[4]])$coef[,2],summary(model[[5]])$coef[,2],summary(model[[6]])$coef[,2],
                summary(model[[7]])$coef[,2],summary(model[[8]])$coef[,2],summary(model[[9]])$coef[,2],
                summary(model[[10]])$coef[,2]), byrow=F, ncol=10)

  Ubar <-rowMeans(U)  #within imputation variance
  B <-(1/(n.impute-1)) * rowSums((Q - Qbar)^2)  #between imputation variance
  T. <-Ubar + (1+ 1/n.impute)*B   #Total variance
  d.f. <- (n.impute-1) * (1+ n.impute*Ubar/ ((n.impute+1) * B))^2   #degree of freedom
  Q.names <- names(model[[1]]$coef)
  output <- data.frame(OR = exp(Qbar),
  		                   lowerCI = exp(Qbar - qt(0.975, df=d.f.)*sqrt(T.)),
  	                     upperCI = exp(Qbar + qt(0.975, df=d.f.)*sqrt(T.)),
  	                     row.names=Q.names)
  round(output,2)
}

# n
tab[,1] <- c(table(cut(dat$age,c(0,8,18))),table(dat$group),NA,NA,table(cut(dat$prevax_calday,c(0,32,100))),table(cut(dat$aftervactime,c(0,28,100))))
# model estiamtes
tab[c(1,3,5,7,9),c(2,5,8,11)] <- 1.00; tab[5,8] <- NA
for(i in 1:4){
  if(i!=3) tab[c(2,4,6,8,10),2:4+3*(i-1)] <- as.matrix(combine.mi(pre[i],postv[i],i,m)[-1,])
  else tab[c(2,4,8,10),2:4+3*(i-1)] <- as.matrix(combine.mi(pre[i],postv[i],i,m)[-1,])
}

tab

#
# End of script.
#


