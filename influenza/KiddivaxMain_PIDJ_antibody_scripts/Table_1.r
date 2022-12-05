#
# R syntax to reproduce information for Table 1 from:
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

# Multiple Imputation

set.seed(12347)
m <- 10
dat.i <- aregImpute( ~I(group) + I(age) + I(male) + I(aftervactime) + I(prevax_calday) 
                      + I(prevax.sH1.c) + I(prevax.sH3.c) + I(prevax.pH1.c) + I(prevax.B.Floride.c) + I(prevax.B.Brisbane.c)
                      + I(postvax.sH1.c) + I(postvax.sH3.c) + I(postvax.pH1.c) + I(postvax.B.Floride.c) + I(postvax.B.Brisbane.c)
                      + I(pilot.ph1s) + I(pilot.sh1) + I(pilot.sh3) + I(pilot.ili) + I(pilot.ari) 
                      + I(pilot.sh1.pre.c) + I(pilot.sh3.pre.c) + I(pilot.sh1.postv.c) +  I(pilot.sh3.postv.c), data=dat, n.impute=m)

dat.nomiss <- list(dat, dat, dat, dat, dat, dat, dat, dat, dat, dat)

for(i in 1:m){
    dat.nomiss[[i]]$postvax.sH1[is.na(dat.nomiss[[i]]$postvax.sH1)] <- 2.5 * (2^dat.i$imputed$postvax.sH1.c[,i])
    dat.nomiss[[i]]$postvax.sH3[is.na(dat.nomiss[[i]]$postvax.sH3)] <- 2.5 * (2^dat.i$imputed$postvax.sH3.c[,i])
    dat.nomiss[[i]]$postvax.pH1[is.na(dat.nomiss[[i]]$postvax.pH1)] <- 2.5 * (2^dat.i$imputed$postvax.pH1.c[,i])
    dat.nomiss[[i]]$postvax.B.Brisbane[is.na(dat.nomiss[[i]]$postvax.B.Brisbane)] <- 2.5 * (2^dat.i$imputed$postvax.B.Brisbane.c[,i])
}

#
#

tab <- matrix(rep(NA,48*10), ncol=10,
              dimnames=list(c("sH1:GMT","","    titer>=40","","    titer>=160","","    GMTrise","","    ratio>=4","","    ratio>=8","",
                              "sH3:GMT","","    titer>=40","","    titer>=160","","    GMTrise","","    ratio>=4","","    ratio>=8","",
                              "sB: GMT","","    titer>=40","","    titer>=160","","    GMTrise","","    ratio>=4","","    ratio>=8","",
                              "pH1:GMT","","    titer>=40","","    titer>=160","","    GMTrise","","    ratio>=4","","    ratio>=8",""),
                            c("Group","n","Predose 1","p","Postdose 1","p","Predose 2","p","Postdose 2","p")))
tab[,1] <- rep(c("TIV-TIV","PL-TIV"),24)
tab[,2] <- rep(table(dat$group)[2:1],24)

# functions
GMT <- function(titer,miss){
  if(miss==0){
      tmp40 <- table(dat[titer]>=40,dat$group,exclude=NULL)[1:2,1:2]; index40 <- 1*(sum(tmp40[1,])*sum(tmp40[,1])/sum(tmp)<5|sum(tmp40[1,])*sum(tmp40[,2])/sum(tmp)<5|sum(tmp40[2,])*sum(tmp40[,1])/sum(tmp)<5|sum(tmp40[2,])*sum(tmp40[,2])/sum(tmp)<5)
      tmp160 <- table(dat[titer]>=160,dat$group,exclude=NULL)[1:2,1:2]; index160 <- 1*(sum(tmp160[1,])*sum(tmp160[,1])/sum(tmp)<5|sum(tmp160[1,])*sum(tmp160[,2])/sum(tmp)<5|sum(tmp160[2,])*sum(tmp160[,1])/sum(tmp)<5|sum(tmp160[2,])*sum(tmp160[,2])/sum(tmp)<5)
      gmt <- matrix(c(round(c(exp(mean(log(dat[dat$group=="TIVTIV",titer]))),exp(mean(log(dat[dat$group=="placeboTIV",titer]))))),
                      round(c(mean(dat[dat$group=="TIVTIV",titer]>=40),mean(dat[dat$group=="placeboTIV",titer]>=40),
                             mean(dat[dat$group=="TIVTIV",titer]>=160),mean(dat[dat$group=="placeboTIV",titer]>=160)),2),
                      round(c(wilcox.test(dat[dat$group=="TIVTIV",titer],dat[dat$group=="placeboTIV",titer])$p.value,NA,
                              sum(fisher.test(tmp40)$p.value*index40,chisq.test(tmp40)$p.value*(1-index40),na.rm=T),NA,
                              sum(fisher.test(tmp160)$p.value*index160,chisq.test(tmp160)$p.value*(1-index160),na.rm=T),NA),2)),ncol=2)
  }
  if(miss==1){
      # Geometric mean titer
      e11m <- e10m <- e11sd <- e10sd <- rep(NA,m)
      for(k in 1:m){
         e11m[k] <- mean(log(dat.nomiss[[k]][dat$group=="TIVTIV",titer]))
         e10m[k] <- mean(log(dat.nomiss[[k]][dat$group=="placeboTIV",titer]))
         e11sd[k] <- sd(log(dat.nomiss[[k]][dat$group=="TIVTIV",titer]))
         e10sd[k] <-sd(log(dat.nomiss[[k]][dat$group=="placeboTIV",titer]))      
      }
      n11 <- sum(dat$group=="TIVTIV"); n10 <- sum(dat$group=="placeboTIV")
    	mean11 <- mean(e11m);	mean10 <- mean(e10m)
    	meanQ <- mean(e11m-e10m)
    	Ubar <- mean(e11sd^2/n11+e10sd^2/n10)
    	B <- sum((e11m-e10m-meanQ)^2)/(m-1)
	    T <- (1 + 1/m) * B + Ubar
    	degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
    	gmt1 <- cbind(round(c(exp(mean11),exp(mean10))),round(c(2*(1 - pt(abs(meanQ)/sqrt(T),df=degf)),NA),2))

      # Proportion titer >=40 or 160
      prop.titer <- function(cutoff){
        mean11 <- mean10 <- rep(NA,m)
        for(i in 1:m){
          mean11[i] <- mean(dat.nomiss[[i]][dat.nomiss[[i]]$group=="TIVTIV",titer]>=cutoff)
          mean10[i] <- mean(dat.nomiss[[i]][dat.nomiss[[i]]$group=="placeboTIV",titer]>=cutoff)
        }                                                                                                                                    
        w <- rep(NA,m)
        for (i in 1:m){
          ctab <- table(dat.nomiss[[i]][,titer]>=cutoff,dat.nomiss[[i]]$group=="TIVTIV")
          w[i] <- chisq.test(ctab)$statistic
        }     
        k <- (nrow(ctab)-1)*(ncol(ctab)-1)
        r2 <- (m+1)*sum((sqrt(w)-sqrt(mean(w)))^2)/(m*(m-1))
        W2 <- (mean(w)/k-(m+1)*r2/(m-1))/(1+r2)
        v2 <- k^(-3/m)*(m-1)*(1+1/r2)^2     # degree of freedom (combined chi-square test)
        p.value = round(pf(W2,df1=k,df2=v2,lower.tail=FALSE),2)
        result <- round(matrix(c(mean(mean11),mean(mean10),p.value,NA),ncol=2),2)
      }
      gmt <- rbind(gmt1,prop.titer(40),prop.titer(160))  
  }
  gmt
}

GMTR <- function(titer1,titer2,miss){
  if(miss==0){
      tmp4 <- table(dat[titer2]/dat[titer1]>=4,dat$group,exclude=NULL)[1:2,1:2]; index4 <- 1*(sum(tmp4[1,])*sum(tmp4[,1])/sum(tmp)<5|sum(tmp4[1,])*sum(tmp4[,2])/sum(tmp)<5|sum(tmp4[2,])*sum(tmp4[,1])/sum(tmp)<5|sum(tmp4[2,])*sum(tmp4[,2])/sum(tmp)<5)
      tmp8 <- table(dat[titer2]/dat[titer1]>=8,dat$group,exclude=NULL)[1:2,1:2]; index8 <- 1*(sum(tmp8[1,])*sum(tmp8[,1])/sum(tmp)<5|sum(tmp8[1,])*sum(tmp8[,2])/sum(tmp)<5|sum(tmp8[2,])*sum(tmp8[,1])/sum(tmp)<5|sum(tmp8[2,])*sum(tmp8[,2])/sum(tmp)<5)
      gmtr <- cbind(c(round(c(exp(mean(log(dat[dat$group=="TIVTIV",titer2]/dat[dat$group=="TIVTIV",titer1]))),exp(mean(log(dat[dat$group=="placeboTIV",titer2]/dat[dat$group=="placeboTIV",titer1])))),1),
                      round(c(mean(dat[dat$group=="TIVTIV",titer2]/dat[dat$group=="TIVTIV",titer1]>=4),mean(dat[dat$group=="placeboTIV",titer2]/dat[dat$group=="placeboTIV",titer1]>=4),
                              mean(dat[dat$group=="TIVTIV",titer2]/dat[dat$group=="TIVTIV",titer1]>=8),mean(dat[dat$group=="placeboTIV",titer2]/dat[dat$group=="placeboTIV",titer1]>=8)),2)),
                      round(c(wilcox.test(dat[dat$group=="TIVTIV",titer2]/dat[dat$group=="TIVTIV",titer1],dat[dat$group=="placeboTIV",titer2]/dat[dat$group=="placeboTIV",titer1])$p.value,NA,
                              sum(fisher.test(tmp4)$p.value*index4,chisq.test(tmp4)$p.value*(1-index4),na.rm=T),NA,
                              sum(fisher.test(tmp8)$p.value*index8,chisq.test(tmp8)$p.value*(1-index8),na.rm=T),NA),2))
  }
  if(miss==1){
      # Geometric mean titer ratio
      e11m <- e10m <- e11sd <- e10sd <- rep(NA,m)
      for (k in 1:m) {
        e11m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$group=="TIVTIV",titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$group=="TIVTIV",titer1]))
        e10m[k] <- mean(log(dat.nomiss[[k]][dat.nomiss[[k]]$group=="placeboTIV",titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$group=="placeboTIV",titer1]))
        e11sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$group=="TIVTIV",titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$group=="TIVTIV",titer1]))
        e10sd[k] <- sd(log(dat.nomiss[[k]][dat.nomiss[[k]]$group=="placeboTIV",titer2]/dat.nomiss[[k]][dat.nomiss[[k]]$group=="placeboTIV",titer1]))
      }
      n11<- sum(dat$group=="TIVTIV"); n10 <- sum(dat$group=="placeboTIV")
    	mean11 <- mean(e11m);	mean10 <- mean(e10m)
    	meanQ <- mean(e11m-e10m)
    	Ubar <- mean(e11sd^2/n11+e10sd^2/n10)
    	B <- sum((e11m-e10m-meanQ)^2)/(m-1)
    	T <- (1 + 1/m) * B + Ubar
    	degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
    	gmtr1 <- cbind(round(c(exp(mean11),exp(mean10)),1),round(c(2*(1 - pt(abs(meanQ)/sqrt(T),df=degf)),NA),2))

      # Proportion titer ratio >=4 or 8
      prop.titer.ratio <- function(cutoff){
        mean11 <- mean10 <- rep(NA,m)
        for(i in 1:m){
          mean11[i] <- mean(dat.nomiss[[i]][dat.nomiss[[i]]$group=="TIVTIV",titer2]/dat.nomiss[[i]][dat.nomiss[[i]]$group=="TIVTIV",titer1]>=cutoff)
          mean10[i] <- mean(dat.nomiss[[i]][dat.nomiss[[i]]$group=="placeboTIV",titer2]/dat.nomiss[[i]][dat.nomiss[[i]]$group=="placeboTIV",titer1]>=cutoff)
        }                                                                                                                                    
        w <- rep(NA,m)
        for (i in 1:m){
          ctab <- table(dat.nomiss[[i]][,titer2]/dat.nomiss[[i]][,titer1]>=cutoff,dat.nomiss[[i]]$group=="TIVTIV")
          w[i] <- chisq.test(ctab)$statistic
        }     
        k <- (nrow(ctab)-1)*(ncol(ctab)-1)
        r2 <- (m+1)*sum((sqrt(w)-sqrt(mean(w)))^2)/(m*(m-1))
        W2 <- (mean(w)/k-(m+1)*r2/(m-1))/(1+r2)
        v2 <- k^(-3/m)*(m-1)*(1+1/r2)^2     # degree of freedom (combined chi-square test)
        p.value = round(pf(W2,df1=k,df2=v2,lower.tail=FALSE),2)
        result <- round(matrix(c(mean(mean11),mean(mean10),p.value,NA),ncol=2),2)
      }
      gmtr <- rbind(gmtr1,prop.titer.ratio(4),prop.titer.ratio(8))    
  }
  gmtr
}

## construct the table
pre <- c("pilot.sh1.pre","pilot.sh3.pre","pilot.FluB.Florida.pre","pilot.ph1.pre","prevax.sH1","prevax.sH3","prevax.B.Brisbane","prevax.pH1" ) 
postv <- c("pilot.sh1.postv","pilot.sh3.postv","pilot.FluB.Florida.postv","pilot.ph1.postv","postvax.sH1","postvax.sH3","postvax.B.Brisbane","postvax.pH1")

for(i in 1:4){
  tab[1:6+12*(i-1),3:4] <- GMT(pre[i],0)
  tab[1:6+12*(i-1),5:6] <- GMT(postv[i],0)
  tab[1:6+12*(i-1),7:8] <- GMT(pre[i+4],0)
  tab[1:6+12*(i-1),9:10] <- GMT(postv[i+4],1)
  tab[7:12+12*(i-1),5:6] <- GMTR(pre[i],postv[i],0)
  tab[7:12+12*(i-1),9:10] <- GMTR(pre[i+4],postv[i+4],1)
}

tab

#
# End of script.
#
