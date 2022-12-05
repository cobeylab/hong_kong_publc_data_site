#
# R syntax to reproduce information for Table S2 from:
#
# Ng S, Ip DKM, Fang VJ, et al.
# The effect of age and recent influenza vaccination history on the immunogenicity and efficacy
# of 2009-10 seasonal trivalent inactivated influenza vaccination in children
# PLoS ONE 2013 (in press).
#
# Last updated by Ng S, Fang VJ and Cowling BJ.
# March 8, 2013
#

require (Hmisc)
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_PLoSONE_scripts/dataframe.r")

dat <- dat[dat$age%in%9:17,]

# MI
set.seed(1)
m=10
dat.i <- aregImpute(  ~ I(age)+I(chron2)+I(vac0910.pre)+I(vac0809)+I(vac0708)+I(swab.pH1)+I(swab.sH3)+I(swab.B)
                        +I(prevax.sH1.c) +I(prevax.sH3.c)+I(prevax.B.Brisbane.c)+I(prevax.pH1.c)
                        +I(postvax.sH1.c) +I(postvax.sH3.c)+I(postvax.B.Brisbane.c)+I(postvax.pH1.c)
                        +I(post.season.sH1.c) +I(post.season.sH3.c)+I(post.season.B.Brisbane.c)+I(post.season.pH1.c), data=dat, n.impute=m)

dat.nomiss <- list(dat, dat, dat, dat, dat, dat, dat, dat, dat, dat)

for(i in 1:m){
    dat.nomiss[[i]]$vac0809[is.na(dat.nomiss[[i]]$vac0809)] <- dat.i$imputed$vac0809[,i]
    dat.nomiss[[i]]$vac0708[is.na(dat.nomiss[[i]]$vac0708)] <- dat.i$imputed$vac0708[,i]
    dat.nomiss[[i]]$prevax.sH1[is.na(dat.nomiss[[i]]$prevax.sH1.c)] <- (2^dat.i$imputed$prevax.sH1.c[,i])*2.5
    dat.nomiss[[i]]$prevax.sH3[is.na(dat.nomiss[[i]]$prevax.sH3.c)] <- (2^dat.i$imputed$prevax.sH3.c[,i])*2.5
    dat.nomiss[[i]]$prevax.B.Brisbane[is.na(dat.nomiss[[i]]$prevax.B.Brisbane.c)] <- (2^dat.i$imputed$prevax.B.Brisbane.c[,i])*2.5
    dat.nomiss[[i]]$prevax.pH1[is.na(dat.nomiss[[i]]$prevax.pH1.c)] <- (2^dat.i$imputed$prevax.pH1.c[,i])*2.5
    dat.nomiss[[i]]$postvax.sH1[is.na(dat.nomiss[[i]]$postvax.sH1.c)] <- (2^dat.i$imputed$postvax.sH1.c[,i])*2.5
    dat.nomiss[[i]]$postvax.sH3[is.na(dat.nomiss[[i]]$postvax.sH3.c)] <- (2^dat.i$imputed$postvax.sH3.c[,i])*2.5
    dat.nomiss[[i]]$postvax.B.Brisbane[is.na(dat.nomiss[[i]]$postvax.B.Brisbane.c)] <- (2^dat.i$imputed$postvax.B.Brisbane.c[,i])*2.5
    dat.nomiss[[i]]$postvax.pH1[is.na(dat.nomiss[[i]]$postvax.pH1.c)] <- (2^dat.i$imputed$postvax.pH1.c[,i])*2.5
}
for(i in 1:m){dat.nomiss[[i]] <- dat.nomiss[[i]][dat.nomiss[[i]]$intervention=="TIV",]}
dat <- dat[dat$intervention=="TIV",]

# functions
GMT.prop <- function(titer){

      # Geometric mean titer
      e0m <- e1m <- e0sd <- e1sd <- n0 <- n1 <- rep(NA,m); res <- list(NA,NA,NA)
      for(comp in 1:3){

        for(k in 1:m){
           e0m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer]))
           e0sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer]))
           n0[k] <- sum(dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0)
           if(comp==1){
             e1m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer]))
             e1sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer]))
             n1[k] <- sum(dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0)
           }
           else if(comp==2){
             e1m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer]))
             e1sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer]))
             n1[k] <- sum(dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1)
           }
           else if(comp==3){
             e1m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer]))
             e1sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer]))
             n1[k] <- sum(dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1)
           }
        }

        mean1 <- mean(e1m);	mean0 <- mean(e0m)
      	meanQ <- mean(e1m-e0m)
      	Ubar <- mean(e1sd^2/n1+e0sd^2/n0)
      	B <- sum((e1m-e0m-meanQ)^2)/(m-1)
  	    T <- (1 + 1/m) * B + Ubar
      	degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
      	res[[comp]] <- cbind(round(exp(mean1)),round(2*(1 - pt(abs(meanQ)/sqrt(T),df=degf)),2))
      }
      gmt <- cbind(round(exp(mean0)),res[[1]],res[[2]],res[[3]])

      # Proportion titer >=40
      mean1 <- mean0 <- w <- rep(NA,m); res2 <- list(NA,NA,NA)
      for(comp in 1:3){
          for(k in 1:m){
             mean0[k] <- mean(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer]>=40)
             if(comp==1){
               mean1[k] <- mean(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer]>=40)
               ctab <- table(dat.nomiss[[k]][dat.nomiss[[k]]$vac0809==0,titer]>=40,dat.nomiss[[k]]$vac0708[dat.nomiss[[k]]$vac0809==0])
             }
             if(comp==2){
               mean1[k] <- mean(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer]>=40)
               ctab <- table(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0,titer]>=40,dat.nomiss[[k]]$vac0809[dat.nomiss[[k]]$vac0708==0])
             }
             if(comp==3){
               mean1[k] <- mean(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer]>=40)
               ctab <- table(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708+dat.nomiss[[k]]$vac0809!=1,titer]>=40,
                             dat.nomiss[[k]]$vac0708[dat.nomiss[[k]]$vac0708+dat.nomiss[[k]]$vac0809!=1])
             }
             w[k] <- chisq.test(ctab)$statistic
          }
          k2 <- (nrow(ctab)-1)*(ncol(ctab)-1)
          r2 <- (m+1)*sum((sqrt(w)-sqrt(mean(w)))^2)/(m*(m-1))
          W2 <- (mean(w)/k2-(m+1)*r2/(m-1))/(1+r2)
          v2 <- k2^(-3/m)*(m-1)*(1+1/r2)^2     # degree of freedom (combined chi-square test)
          p.value = round(pf(W2,df1=k2,df2=v2,lower.tail=FALSE),2)
          res[[comp]] <- c(mean(mean1),p.value)
      }
      prop <- round(c(mean(mean0),res[[1]],res[[2]],res[[3]]),2)

      gmt.prop <- rbind(gmt,prop)
      gmt.prop
}

GMTR <- function(titer1,titer2){
      # Geometric mean titer ratio
      e0m <- e1m <- e0sd <- e1sd <- n0 <- n1 <- rep(NA,m); res <- list(NA,NA,NA)
      for(comp in 1:3){

        for(k in 1:m){
           e0m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer1]))
           e0sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0,titer1]))
           n0[k] <- sum(dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==0)
           if(comp==1){
             e1m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer1]))
             e1sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0,titer1]))
             n1[k] <- sum(dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==0)
           }
           else if(comp==2){
             e1m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer1]))
             e1sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==0&dat.nomiss[[k]]$vac0809==1,titer1]))
             n1[k] <- sum(dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1)
           }
           else if(comp==3){
             e1m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer1]))
             e1sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1,titer1]))
             n1[k] <- sum(dat.nomiss[[k]]$vac0708==1&dat.nomiss[[k]]$vac0809==1)
           }
        }

        mean1 <- mean(e1m);	mean0 <- mean(e0m)
      	meanQ <- mean(e1m-e0m)
      	Ubar <- mean(e1sd^2/n1+e0sd^2/n0)
      	B <- sum((e1m-e0m-meanQ)^2)/(m-1)
  	    T <- (1 + 1/m) * B + Ubar
      	degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
      	res[[comp]] <- cbind(round(exp(mean1)),round(2*(1 - pt(abs(meanQ)/sqrt(T),df=degf)),2))
      }

      gmtr <- cbind(round(exp(mean0)),res[[1]],res[[2]],res[[3]])
      gmtr
}

## construct the table
tab <- matrix(rep(NA,21*7), ncol=7,
              dimnames=list(c("n","sH1:preGMT","titer>=40","postGMT","titer>=40","GMTR","sH3:preGMT","titer>=40","postGMT","titer>=40","GMTR",
                              "B:preGMT","titer>=40","postGMT","titer>=40","GMTR","pH1:preGMT","titer>=40","postGMT","titer>=40","GMTR"),
                            c("Reference","Compare1","p","Compare2","p","Compare3","p")))
tab[1,c(1,2,4,6)] <- c(table(dat$vac0708,dat$vac0809))

pre <- c("prevax.sH1","prevax.sH3","prevax.B.Brisbane","prevax.pH1")
postv <- c("postvax.sH1","postvax.sH3","postvax.B.Brisbane","postvax.pH1")

for(i in 1:4){
  tab[2:3+(i-1)*5,] <- GMT.prop(pre[i])
  tab[4:5+(i-1)*5,] <- GMT.prop(postv[i])
  tab[6+(i-1)*5,] <- GMTR(pre[i],postv[i])
}

tab

#
# End of script.
#


