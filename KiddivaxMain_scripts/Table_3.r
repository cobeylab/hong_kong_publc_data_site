#
# R syntax to reproduce Table 3 from:
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

tab <- matrix(rep(NA,9*10), ncol=10, dimnames=list(
              c("Vaccinees","pcr.ph1","pcr.sh3","pcr.b","sero.ph1","sero.sh3","sero.b", "ARI","ILI"),
              c("TIV","CI_low","CI_up", "placebo", "CI_low","CI_up", "p.value","VE","CI_low","CI_up")))

index.pre$agegp <- cut(index.pre$age,c(0,9,17))
index.pre$followup <- as.numeric(as.Date(as.character(index.pre$end.date),format="%d/%m/%Y")-as.Date(as.character(index.pre$start.date),format="%d/%m/%Y"))

set.seed(12345)
index.i <- aregImpute( ~ factor(intervention) +
                         pH1 + sH3 + B + swab.pH1 + swab.sH3 + swab.B + I(nARI) + I(nILI) +
                         male + I(agegp)+ vac0809 + followup, data=index.pre, n.impute=10)

index.nomiss <- list(index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre)

for(i in 1:10){
    index.nomiss[[i]]$pH1[is.na(index.nomiss[[i]]$pH1)] <- index.i$imputed$pH1[,i]
    index.nomiss[[i]]$sH3[is.na(index.nomiss[[i]]$sH3)] <- index.i$imputed$sH3[,i]
    index.nomiss[[i]]$B[is.na(index.nomiss[[i]]$B)] <- index.i$imputed$B[,i]
}

cumi2 <- function(group,followup,nevent){
  da <- data.frame(gp=group,len=followup,n=nevent)
  n1 <- sum(da$n[da$gp=="TIV"]); N1 <- sum(da$len[da$gp=="TIV"]/365);
  n2 <- sum(da$n[da$gp=="placebo"]); N2 <- sum(da$len[da$gp=="placebo"]/365)
  ci1 <- n1/N1;  ci2 <- n2/N2

  da$len[da$len==0] <- 0.1
  m11 <- glm(n ~ as.factor(gp)+offset(log(len)), family="poisson", data=da)
  m12 <- glm(n ~ offset(log(len)), family="poisson", data=da)
  dev.dff <- anova(m12, m11); p.value <- 1-pchisq(dev.dff[2,4], dev.dff[2,3])
  #
  m1 <- glm(n ~ as.factor(gp)+offset(log(len)), family="poisson", data=da); 
  da$gp2 <- 1*(da$gp=="placebo") ; m2 <- glm(n ~ as.factor(gp2)+offset(log(len)), family="poisson", data=da)
  ci1.CI <- exp(c(m2$coef[1]-1.96*sqrt(diag(vcov(m2)))[1]+log(365),m2$coef[1]+1.96*sqrt(diag(vcov(m2)))[1]+log(365)))
  ci2.CI <- exp(c(m1$coef[1]-1.96*sqrt(diag(vcov(m1)))[1]+log(365),m1$coef[1]+1.96*sqrt(diag(vcov(m1)))[1]+log(365)))
  
  RR <- exp(m1$coef[2]); RR.CI <- exp(c(m1$coef[2]-1.96*sqrt(diag(vcov(m1)))[2],m1$coef[2]+1.96*sqrt(diag(vcov(m1)))[2]))
  p.value <- summary(m11)$coef[2,4]
  round(c(exp(m2$coef[1]+log(365)),ci1.CI,exp(m1$coef[1]+log(365)),ci2.CI,p.value,1-RR,rev(1-RR.CI)),2)  
}

# function to combined the imputed results
cumi2.mi <- function(data,m){
  model1 <- model2 <- list(NA); beta1 <- beta2 <- RR <- var1 <- var2 <- var.RR <- NA
  for(i in 1:m){
    data[[i]]$followup[data[[i]]$followup==0] <- 0.1;     data[[i]]$group <- 1*(data[[i]]$intervention=="placebo")
    model1[[i]] <- glm(event~as.factor(group)+offset(log(followup)), family="poisson", data=data[[i]])            
    model2[[i]] <- glm(event~as.factor(intervention)+offset(log(followup)), family="poisson", data=data[[i]])
    beta1[i] <- summary(model1[[i]])$coef[1,1]; beta2[i] <- summary(model2[[i]])$coef[1,1]; RR[i] <- summary(model2[[i]])$coef[2,1]
    var1[i] <- summary(model1[[i]])$coef[1,2]^2; var2[i] <- summary(model2[[i]])$coef[1,2]^2; var.RR[i] <- summary(model2[[i]])$coef[2,2]^2 
  }
  T1 <- (1+1/m)*(sum((beta1-mean(beta1))^2)/(m-1))+mean(var1); 
  degf1 <- (m-1)*(1 +mean(var1)/((1+1/m)*(sum((beta1-mean(beta1))^2)/(m-1))))*(1+mean(var1)/((1+1/m)*(sum((beta1-mean(beta1))^2)/(m-1))))
  ci1.CI <- c(exp(mean(beta1)+log(365)-qt(0.975,df=degf1)*sqrt(T1)),exp(mean(beta1)+log(365)+qt(0.975,df=degf1)*sqrt(T1)))
  T2 <- (1+1/m)*(sum((beta2-mean(beta2))^2)/(m-1))+mean(var2); 
  degf2 <- (m-1)*(1 +mean(var2)/((1+1/m)*(sum((beta2-mean(beta2))^2)/(m-1))))*(1+mean(var2)/((1+1/m)*(sum((beta2-mean(beta2))^2)/(m-1))))
  ci2.CI <- c(exp(mean(beta2)+log(365)-qt(0.975,df=degf2)*sqrt(T2)),exp(mean(beta2)+log(365)+qt(0.975,df=degf2)*sqrt(T2)))

  T.RR <- (1+1/m)*(sum((RR-mean(RR))^2)/(m-1))+mean(var.RR); 
  degf.RR <- (m-1)*(1 +mean(var.RR)/((1+1/m)*(sum((RR-mean(RR))^2)/(m-1))))*(1+mean(var.RR)/((1+1/m)*(sum((RR-mean(RR))^2)/(m-1))))
  RR.CI <- c(exp(mean(RR)-qt(0.975,df=degf.RR)*sqrt(T.RR)),exp(mean(RR)+qt(0.975,df=degf.RR)*sqrt(T.RR)))
  p.value <- 2*(1-pt(abs(mean(RR)/sqrt(T.RR)),df=degf.RR))
  round(c(exp(mean(beta1)+log(365)),ci1.CI,exp(mean(beta2)+log(365)),ci2.CI,p.value,1-exp(mean(RR)),rev(1-RR.CI)),2)  
}

temp <- index.nomiss

tab[1,c(1,4)] <- rev(table(index.pre$intervention))
tab[2,1:10] <- cumi2(index.pre$intervention,index.pre$followup,index.pre$swab.pH1)
tab[3,1:10] <- cumi2(index.pre$intervention,index.pre$followup,index.pre$swab.sH3)
tab[4,1:10] <- cumi2(index.pre$intervention,index.pre$followup,index.pre$swab.B)

for(i in 1:10){temp[[i]]$event <- temp[[i]]$pH1}; tab[5,1:10] <- cumi2.mi(temp,10)
for(i in 1:10){temp[[i]]$event <- temp[[i]]$sH3}; tab[6,1:10] <- cumi2.mi(temp,10)
for(i in 1:10){temp[[i]]$event <- temp[[i]]$B}; tab[7,1:10] <- cumi2.mi(temp,10)   

tab[8,1:10] <- cumi2(index.pre$intervention,index.pre$followup,index.pre$nARI)
tab[9,1:10] <- cumi2(index.pre$intervention,index.pre$followup,index.pre$nILI)

tab

#
# End of script.
#

