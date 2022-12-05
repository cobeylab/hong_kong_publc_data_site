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

source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_scripts/dataframe.r")

cumi <- function(X){
  p1 <- X[2,2]/sum(X[2,]); p2 <- X[1,2]/sum(X[1,])
  CI1 <- binom.test(X[2,2],sum(X[2,]))$conf[1:2]; CI2 <- binom.test(X[1,2],sum(X[1,]))$conf[1:2]
  p.value <- chisq.test(X)$p.value; p.value2 <- fisher.test(X)$p.value
  round(c(sum(X[2,]),p1,CI1,sum(X[1,]),p2,CI2,p.value,p.value2),2)
}

tab <- matrix(rep(NA,4*10), ncol=10, dimnames=list(
              c("Postvax titer against pH1N1: <=1:20","1:40","1:80",">=1:160"),
              c("TIV","CumInc","CI_low","CI_up", "placebo","CumInc","CI_low","CI_up", "p.value1","p.value2")))

table(index.pre$intervention[!is.na(index.pre$pH1)])
tmp <- index.pre[index.pre$postvax.pH1<=20&!is.na(index.pre$postvax.pH1),]; tab[1,1:10] <- cumi(table(tmp$intervention,tmp$pH1))
tmp <- index.pre[index.pre$postvax.pH1==40&!is.na(index.pre$postvax.pH1),]; tab[2,1:10] <- cumi(table(tmp$intervention,tmp$pH1))
tmp <- index.pre[index.pre$postvax.pH1==80&!is.na(index.pre$postvax.pH1),]; tab[3,1:10] <- cumi(table(tmp$intervention,tmp$pH1))
tmp <- index.pre[index.pre$postvax.pH1>=160&!is.na(index.pre$postvax.pH1),]; tab[4,1:10] <- cumi(cbind(table(tmp$intervention,tmp$pH1),c(0,0)))

tab

#
# End of script.
#
