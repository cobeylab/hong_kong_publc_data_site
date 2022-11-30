#
# R syntax to reproduce Table S3 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

require(Hmisc)
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_scripts/dataframe.r")

mid2 <- index.pre[!is.na(index.pre$mid.season.pH1)&!is.na(index.pre$mid.season.sH3)&!is.na(index.pre$mid.season.B.Brisbane),]; mid2$agegp <- cut(mid2$age,c(0,9,17))

 #------
 mid2$infect12.pH1 <- 1*(mid2$mid.season.pH1/mid2$postvax.pH1>=4)        #(1)
 mid2$infect23.pH1 <- 1*(mid2$post.season.pH1/mid2$mid.season.pH1>=4)    #(2)
 mid2$infect13.pH1 <- 1*(mid2$post.season.pH1/mid2$postvax.pH1>=4)       #(4)
 #------
 mid2$infect12.sH3 <- 1*(mid2$mid.season.sH3/mid2$postvax.sH3>=4)
 mid2$infect23.sH3 <- 1*(mid2$post.season.sH3/mid2$mid.season.sH3>=4)
 mid2$infect13.sH3 <- 1*(mid2$post.season.sH3/mid2$postvax.sH3>=4)
 #------
 mid2$infect12.B.Brisbane <- 1*(mid2$mid.season.B.Brisbane/mid2$postvax.B.Brisbane>=4)
 mid2$infect23.B.Brisbane <- 1*(mid2$post.season.B.Brisbane/mid2$mid.season.B.Brisbane>=4)
 mid2$infect13.B.Brisbane <- 1*(mid2$post.season.B.Brisbane/mid2$postvax.B.Brisbane>=4)
 #------

set.seed(12345)
mid.i <- aregImpute( ~ factor(intervention)+ infect12.pH1+infect23.pH1+infect13.pH1+infect12.sH3+infect23.sH3+infect13.sH3+
                          infect12.B.Brisbane+infect23.B.Brisbane+infect13.B.Brisbane+
                          swab.pH1 + swab.sH3 + swab.B + ARI + ILI +
                          male+ I(agegp)+ vac0809, data=mid2, n.impute=10)

mid.nomiss <- list(mid2, mid2, mid2, mid2, mid2, mid2, mid2, mid2, mid2, mid2)

for(i in 1:10){
    mid.nomiss[[i]]$infect12.pH1[is.na(mid.nomiss[[i]]$infect12.pH1)] <- mid.i$imputed$infect12.pH1[,i]
    mid.nomiss[[i]]$infect23.pH1[is.na(mid.nomiss[[i]]$infect23.pH1)] <- mid.i$imputed$infect23.pH1[,i]
    mid.nomiss[[i]]$infectws.pH1 <- 1*(mid.nomiss[[i]]$infect12.pH1==1|mid.nomiss[[i]]$infect23.pH1==1)
    mid.nomiss[[i]]$infect13.pH1[is.na(mid.nomiss[[i]]$infect13.pH1)] <- mid.i$imputed$infect13.pH1[,i]
    mid.nomiss[[i]]$infect12.sH3[is.na(mid.nomiss[[i]]$infect12.sH3)] <- mid.i$imputed$infect12.sH3[,i]
    mid.nomiss[[i]]$infect23.sH3[is.na(mid.nomiss[[i]]$infect23.sH3)] <- mid.i$imputed$infect23.sH3[,i]
    mid.nomiss[[i]]$infectws.sH3 <- 1*(mid.nomiss[[i]]$infect12.sH3==1|mid.nomiss[[i]]$infect23.sH3==1)
    mid.nomiss[[i]]$infect13.sH3[is.na(mid.nomiss[[i]]$infect13.sH3)] <- mid.i$imputed$infect13.sH3[,i]
    mid.nomiss[[i]]$infect12.B.Brisbane[is.na(mid.nomiss[[i]]$infect12.B.Brisbane)] <- mid.i$imputed$infect12.B.Brisbane[,i]
    mid.nomiss[[i]]$infect23.B.Brisbane[is.na(mid.nomiss[[i]]$infect23.B.Brisbane)] <- mid.i$imputed$infect23.B.Brisbane[,i]
    mid.nomiss[[i]]$infectws.B.Brisbane <- 1*(mid.nomiss[[i]]$infect12.B.Brisbane==1|mid.nomiss[[i]]$infect23.B.Brisbane==1)
    mid.nomiss[[i]]$infect13.B.Brisbane[is.na(mid.nomiss[[i]]$infect13.B.Brisbane)] <- mid.i$imputed$infect13.B.Brisbane[,i]
}

# function to combined the imputed results
combine.mi <- function(data,m){       # data: hhID, intervention, event
  nTIV <- sum(data[[1]]$intervention=="TIV")
  npla <- sum(data[[1]]$intervention=="placebo")

 	# calculate 95% CI
	meanTIV <- varTIV <- meanp <- varp <- rep(NA,m)
  for (i in 1:m){
    meanTIV[i] <- sum(data[[i]]$event[data[[i]]$intervention=="TIV"])/nTIV
    varTIV[i] <- meanTIV[i]*(1-meanTIV[i])/nTIV
    meanp[i] <- sum(data[[i]]$event[data[[i]]$intervention=="placebo"])/npla
    varp[i] <- meanp[i]*(1-meanp[i])/npla
  }
  CI.95 <- function(mu,var){
    m <- length(mu)
    estQ <- mean(mu)
    Ubar <- mean(var)
    B <- sum((mu-estQ)^2)/(m-1)
    T <- (1+1/m)*B+Ubar
    degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
    CI.low <- max(estQ-qt(0.975,df=degf)*sqrt(T),0)
    CI.up <- estQ+qt(0.975,df=degf)*sqrt(T)
    round(c(estQ,CI.low,CI.up),2)
  }

  # calculate p-value
  w <- rep(NA,m)
  for (i in 1:m){
    ctab <- table(data[[i]]$event,data[[i]]$intervention=="placebo")
    w[i] <- chisq.test(ctab)$statistic
  }

  k <- (nrow(ctab)-1)*(ncol(ctab)-1)
  r2 <- (m+1)*sum((sqrt(w)-sqrt(mean(w)))^2)/(m*(m-1))
  W2 <- (mean(w)/k-(m+1)*r2/(m-1))/(1+r2)
  v2 <- k^(-3/m)*(m-1)*(1+1/r2)^2     # degree of freedom (combined chi-square test)
  p.value = round(pf(W2,df1=k,df2=v2,lower.tail=FALSE),2)

  # write the output
  output1 <- c(CI.95(meanTIV,varTIV),CI.95(meanp,varp),p.value)

  # vaccine effect
  logrr <- var <- rep(NA,m)
  for(i in 1:m){
  X <- table(data[[i]]$intervention,data[[i]]$event)
  logrr[i] <- log((X[2,2]/sum(X[2,]))/(X[1,2]/sum(X[1,])))
  var[i] <- 1/X[2,2]-1/sum(X[2,])+1/X[1,2]-1/sum(X[1,])
  }
    m <- length(logrr)
    estQ <- mean(logrr)
    Ubar <- mean(var)
    B <- sum((logrr-estQ)^2)/(m-1)
    T <- (1+1/m)*B+Ubar
    degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
    CI.low <- estQ-qt(0.975,df=degf)*sqrt(T)
    CI.up <- estQ+qt(0.975,df=degf)*sqrt(T)

  # write the output
  output2 <- round(c(1-exp(estQ),1-exp(CI.up),1-exp(CI.low)),2)
  output <- c(output1, output2)
  names(output) <- c("TIV%","CI.low","CI.up","placebo%","CI.low","CI.up","p-value","VE","CI.low","CI.up")
	output
}

tab <- matrix(NA,ncol=10,nrow=12, dimnames=list(
              c("(1)pH1","(1)sH3","(1)B","(2)pH1","(2)sH3","(2)B","(3)pH1","(3)sH3","(3)B","(4)pH1","(4)sH3","(4)B"),
              c("TIV","CI_low","CI_up","Placebo","CI_low","CI_up","p-value","VE","CI_low","CI_up")))


temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect12.pH1}; tab[1,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect12.sH3}; tab[2,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect12.B.Brisbane}; tab[3,] <- combine.mi(temp,10)

temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect23.pH1}; tab[4,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect23.sH3}; tab[5,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect23.B.Brisbane}; tab[6,] <- combine.mi(temp,10)

temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infectws.pH1}; tab[7,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infectws.sH3}; tab[8,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infectws.B.Brisbane}; tab[9,] <- combine.mi(temp,10)

temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect13.pH1}; tab[10,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect13.sH3}; tab[11,] <- combine.mi(temp,10)
temp <- mid.nomiss; for (i in 1:10){temp[[i]]$event <- temp[[i]]$infect13.B.Brisbane}; tab[12,] <- combine.mi(temp,10)

tab

#
# End of script.
#


