#
# R syntax to reproduce information for Table 3 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 11, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe.r")

##################################################
#Multiple Imputation for missing index info      #
##################################################
require(Hmisc)
index.pre$agegp<-cut(index.pre$age,c(0,8,15))

set.seed(2)
index.i <- aregImpute( ~ factor(intervention) 
                         + I(sh1w) + I(sh1s) + I(sh3w) + I(sh3s) + I(ph1s) + I(ariw) + I(iliw) + I(aris) + I(ilis) + I(ilip) + I(day.p)
                         + I(flo) + I(flo.w) + I(flo.s) + I(male) + I(agegp)+ I(vac09)
                          , data=index.pre,n.impute=10)


index.nomiss <- list(index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre)

for(i in 1:10){
    index.nomiss[[i]]$flo[is.na(index.nomiss[[i]]$flo)] <- index.i$imputed$flo[,i]
    index.nomiss[[i]]$flo.w[is.na(index.nomiss[[i]]$flo.w)] <- index.i$imputed$flo.w[,i]
    index.nomiss[[i]]$flo.s[is.na(index.nomiss[[i]]$flo.s)] <- index.i$imputed$flo.s[,i]
    index.nomiss[[i]]$sh1w[is.na(index.nomiss[[i]]$sh1w)] <- index.i$imputed$sh1w[,i]
    index.nomiss[[i]]$sh3w[is.na(index.nomiss[[i]]$sh3w)] <- index.i$imputed$sh3w[,i]
    index.nomiss[[i]]$sh3s[is.na(index.nomiss[[i]]$sh3s)] <- index.i$imputed$sh3s[,i]
    index.nomiss[[i]]$ph1s[is.na(index.nomiss[[i]]$ph1s)] <- index.i$imputed$ph1s[,i]
    index.nomiss[[i]]$sh1s[is.na(index.nomiss[[i]]$sh1s)] <- index.i$imputed$sh1s[,i]
    index.nomiss[[i]]$iliw[is.na(index.nomiss[[i]]$iliw)] <- index.i$imputed$iliw[,i]
    index.nomiss[[i]]$ariw[is.na(index.nomiss[[i]]$ariw)] <- index.i$imputed$ariw[,i]
    index.nomiss[[i]]$ilis[is.na(index.nomiss[[i]]$ilis)] <- index.i$imputed$ilis[,i]
    index.nomiss[[i]]$aris[is.na(index.nomiss[[i]]$aris)] <- index.i$imputed$aris[,i]
    index.nomiss[[i]]$ilip[is.na(index.nomiss[[i]]$ilip)] <- index.i$imputed$ilip[,i]
    index.nomiss[[i]]$day.p[is.na(index.nomiss[[i]]$day.p)] <- index.i$imputed$day.p[,i]
    }


###################################################
#Multiple Imputation for missing member info      #
###################################################
# Re-do agegp for regression
member.pre$agegp<-cut(member.pre$age,c(0,15,45,100))

set.seed(1)
member.i <- aregImpute( ~ factor(intervention)+
                         sh1w + sh1s + sh3w + sh3s + ph1s + ariw + iliw + aris + ilis + ilip + I(day.p) +
                         flo + flo.w + flo.s + male+ I(agegp)+ vac09, data=member.pre,n.impute=10)

member.nomiss <- list(member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre)

for(i in 1:10){
    member.nomiss[[i]]$flo.s[is.na(member.nomiss[[i]]$flo.s)] <- member.i$imputed$flo.s[,i]
    member.nomiss[[i]]$flo.w[is.na(member.nomiss[[i]]$flo.w)] <- member.i$imputed$flo.w[,i]
    member.nomiss[[i]]$flo[is.na(member.nomiss[[i]]$flo)] <- member.i$imputed$flo[,i]
    member.nomiss[[i]]$sh1w[is.na(member.nomiss[[i]]$sh1w)] <- member.i$imputed$sh1w[,i]
    member.nomiss[[i]]$sh3w[is.na(member.nomiss[[i]]$sh3w)] <- member.i$imputed$sh3w[,i]
    member.nomiss[[i]]$sh3s[is.na(member.nomiss[[i]]$sh3s)] <- member.i$imputed$sh3s[,i]
    member.nomiss[[i]]$ph1s[is.na(member.nomiss[[i]]$ph1s)] <- member.i$imputed$ph1s[,i]
    member.nomiss[[i]]$sh1s[is.na(member.nomiss[[i]]$sh1s)] <- member.i$imputed$sh1s[,i]
    member.nomiss[[i]]$iliw[is.na(member.nomiss[[i]]$iliw)] <- member.i$imputed$iliw[,i]
    member.nomiss[[i]]$ariw[is.na(member.nomiss[[i]]$ariw)] <- member.i$imputed$ariw[,i]
    member.nomiss[[i]]$ilis[is.na(member.nomiss[[i]]$ilis)] <- member.i$imputed$ilis[,i]
    member.nomiss[[i]]$aris[is.na(member.nomiss[[i]]$aris)] <- member.i$imputed$aris[,i]
    member.nomiss[[i]]$ilip[is.na(member.nomiss[[i]]$ilip)] <- member.i$imputed$ilip[,i]
    member.nomiss[[i]]$vac09[is.na(member.nomiss[[i]]$vac09)] <- member.i$imputed$vac09[,i]
    member.nomiss[[i]]$day.p[is.na(member.nomiss[[i]]$day.p)] <- member.i$imputed$day.p[,i]
    }


index.pre$gp[index.pre$intervention=="TIV"] <-"1"
index.pre$gp[index.pre$intervention=="placebo"] <-"2"
member.pre$gp[member.pre$intervention=="TIV"] <-"1"
member.pre$gp[member.pre$intervention=="placebo"] <-"2"

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
  output <- c(CI.95(meanTIV,varTIV),CI.95(meanp,varp),p.value)
  names(output) <- c("TIV%","CI.low","CI.up","placebo%","CI.low","CI.up","p.value")
	output
}

tab3 <-matrix(rep(NA,7*22), ncol=7, dimnames=list(
              c("vaccinees", "sero.sh1", "sero.sh3","sero.ph1", "sero.b","pcr.sh1", "pcr.sh3","pcr.ph1","pcr.b", "ili","ari",
               "contacts", "sero.sh1", "sero.sh3", "sero.ph1","sero.b", "pcr.sh1", "pcr.sh3", "pcr.ph1", "pcr.b", "ili", "ari"),
              c("TIV","CI_low","CI_up", "placebo", "CI_low","CI_up", "p.value")))

temp <- index.nomiss

ff <-function (x) {
  for (i in 1:10){
       temp[[i]]$event <- x
       temp[[i]] <- temp[[i]][,c("hhID","intervention","event")]
  }
  combine.mi(temp,10)
}

tab3[1,c(1,4)] <-table(index.pre$gp)
tab3[2,] <-ff(1*(temp[[i]]$sh1w==1 | temp[[i]]$sh1s==1))
tab3[3,] <-ff(1*(temp[[i]]$sh3w==1 | temp[[i]]$sh3s==1))
tab3[4,] <-ff(temp[[i]]$ph1s)
tab3[5,] <-ff(1*(temp[[i]]$flo.w==1 |temp[[i]]$flo.s==1))
tab3[6,] <-ff(temp[[i]]$swab.sh1w==1)
tab3[7,] <-ff(temp[[i]]$swab.sh3w==1)
tab3[8,] <-ff(temp[[i]]$swab.ph1==1)
tab3[8,5:6] <- round(binom.test(sum(index.pre$swab.ph1[index.pre$gp==2]),length(index.pre$swab.ph1[index.pre$gp==2]))$conf[1:2],2)
tab3[9,] <-ff(temp[[i]]$swab.bw==1)
tab3[9,2:3] <- round(binom.test(sum(index.pre$swab.bw[index.pre$gp==1]),length(index.pre$swab.bw[index.pre$gp==1]))$conf[1:2],2)
tab3[10,] <-ff(1*(temp[[i]]$iliw==1 | temp[[i]]$ilis==1))
tab3[11,] <-ff(1*(temp[[i]]$ariw==1 | temp[[i]]$aris==1))

temp <- member.nomiss

tab3[12,c(1,4)] <-table(member.pre$gp)
tab3[13,] <-ff(1*(temp[[i]]$sh1w==1 | temp[[i]]$sh1s==1))
tab3[14,] <-ff(1*(temp[[i]]$sh3w==1 | temp[[i]]$sh3s==1))
tab3[15,] <-ff(temp[[i]]$ph1s)
tab3[16,] <-ff(1*(temp[[i]]$flo.w==1 |temp[[i]]$flo.s==1))
tab3[17,] <-ff(temp[[i]]$swab.sh1w==1)
tab3[18,] <-ff(temp[[i]]$swab.sh3w==1)
tab3[19,] <-ff(temp[[i]]$swab.ph1==1)
tab3[20,] <-ff(temp[[i]]$swab.bw==1)
tab3[20,2:3] <- round(binom.test(sum(member.pre$swab.bw[member.pre$gp==1]),length(member.pre$swab.bw[member.pre$gp==1]))$conf[1:2],2)
tab3[20,5:6] <- round(binom.test(sum(member.pre$swab.bw[member.pre$gp==2]),length(member.pre$swab.bw[member.pre$gp==2]))$conf[1:2],2)
tab3[20,7] <- fisher.test(matrix(c(0,189,0,123),ncol=2))$p.value
tab3[21,] <-ff(1*(temp[[i]]$iliw==1 | temp[[i]]$ilis==1))
tab3[22,] <-ff(1*(temp[[i]]$ariw==1 | temp[[i]]$aris==1))

tab3

#
# End of script.
#
