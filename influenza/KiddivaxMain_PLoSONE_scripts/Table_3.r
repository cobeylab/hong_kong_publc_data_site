#
# R syntax to reproduce information for Table 3 from:
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

swab.Victoria <- unique(swab[swab$FluB.subtype=="Victoria" ,c("hhID","member")]); swab.Victoria$swab.Victoria <- 1
dat <- merge(dat, swab.Victoria, by=c("hhID","member"),all.x=T); dat$swab.Victoria[is.na(dat$swab.Victoria)] <- 0

dat$followup <- as.numeric(as.Date(as.character(dat$end.date),format="%d/%m/%Y")-as.Date(as.character(dat$start.date),format="%d/%m/%Y"))

# count ARI/ILI episodes
cnames <- c("nARI","nILI")
for(k in 1:2){
  if(k==1) epi <- symp[symp$ARI==1,]
  if(k==2) epi <- symp[symp$ILI==1,]
  
  epi$day <- dates(as.character(epi$date),format="d/m/Y")-dates("1/1/2010",format="d/m/Y")
  temp <- unique(epi[1:2])
  temp$subject <- 1:nrow(temp)
  epi <- merge(epi,temp,by=c("hhID","member"),all.x=T)
  
  for (i in 1:nrow(temp)){
     init.row <- nrow(epi[epi$subject<i,])+1
     epi$epi[init.row] <- 1
     j <- 1
     epi.row <- init.row
     while (j <= nrow(epi[epi$subject==i,])-1){
          if(epi$day[init.row+j]-epi$day[init.row+j-1]<7){
             epi$epi[init.row+j] <- epi$epi[init.row+j-1]
          }
          else {
             epi$epi[init.row+j] <- epi$epi[init.row+j-1]+1
             epi.row <- init.row+j
          }
          j <- j+1
     }
     temp$nepi[i] <- epi$epi[init.row+j-1]
  }
  dat <- merge(dat,temp[c(1,4)],all.x=T); dat$nepi[is.na(dat$nepi)] <- 0
  names(dat)[ncol(dat)] <- cnames[k]
}


#### construct the model
tab <- matrix(rep(NA, 7*18), nrow=18, dimnames=list(c("6-8_Model1: Placebo","TIV1","TIV2","Model2: Placebo","TIV1","TIV2","Model3: Placebo","TIV1","TIV2",
                                                     "9-17_Model1: Placebo","TIV1","TIV2","Model2: Placebo","TIV1","TIV2","Model3: Placebo","TIV1","TIV2"),
                                                   c("RR","CI_low","CI_up","p-value","VE","CI_low","CI_up")))

for(k in 1:2){
  if(k==1) dats <- dat[dat$age%in%6:8,]
  else if(k==2) dats <- dat[dat$age%in%9:17,]
  # MI
  set.seed(1)
  m <- 10
  dat.i <- aregImpute( ~  I(age)+I(chron2)+I(male)+I(vac0910.pre)+I(vac0809)+I(vac0708)+I(swab.Victoria)+I(swab.sH3)+I(swab.B)
                          +I(prevax.sH1.c) +I(prevax.sH3.c)+I(prevax.B.Brisbane.c)+I(prevax.pH1.c)
                          +I(postvax.sH1.c) +I(postvax.sH3.c)+I(postvax.B.Brisbane.c)+I(postvax.pH1.c)
                          +I(post.season.sH1.c) +I(post.season.sH3.c)+I(post.season.B.Brisbane.c)+I(post.season.pH1.c)
                          +I(nILI)+I(nARI)+I(followup), data=dats, n.impute=m)
  dat.nomiss <- list(dats, dats, dats, dats, dats, dats, dats, dats, dats, dats) 
  for(i in 1:m){
      dat.nomiss[[i]]$vac0809[is.na(dat.nomiss[[i]]$vac0809)] <- dat.i$imputed$vac0809[,i]
      dat.nomiss[[i]]$vac0708[is.na(dat.nomiss[[i]]$vac0708)] <- dat.i$imputed$vac0708[,i]
  }
    
  # function to combined the imputed results
  cumi2.mi <- function(idata,m,rr,HxDef){
  
    model<- list(NA); RR <- var.RR <- NA
    for(i in 1:m){
      idata[[i]]$followup[idata[[i]]$followup==0] <- 0.1;
      idata[[i]]$va<-NA
      if (HxDef==1) {
                    idata[[i]]$va[idata[[i]]$vac0910.pre==1 & idata[[i]]$vac0708==0] <- 2  
                    idata[[i]]$va[idata[[i]]$vac0910.pre==1 & idata[[i]]$vac0708==1] <- 3
                    }
      
      if (HxDef==2) {
                    idata[[i]]$va[idata[[i]]$vac0910.pre==1 & idata[[i]]$vac0809==0] <- 2
                    idata[[i]]$va[idata[[i]]$vac0910.pre==1 & idata[[i]]$vac0809==1] <- 3
                    }
  
      if (HxDef==3) {
                    idata[[i]]$va[idata[[i]]$vac0910.pre==1 & idata[[i]]$vac0708==0 & idata[[i]]$vac0809==0] <- 2  
                    idata[[i]]$va[idata[[i]]$vac0910.pre==1 & idata[[i]]$vac0708+idata[[i]]$vac0809!=0] <- 3
                    }  
      idata[[i]]$va[idata[[i]]$vac0910.pre==0]<-1
      
      model[[i]] <- glm(swab.Victoria~factor(va)+chron2+offset(log(followup)), family="poisson", data=idata[[i]])
      
      RR[i] <- summary(model[[i]])$coef[rr,1]
      var.RR[i] <- summary(model[[i]])$coef[rr,2]^2
    }
  
    T.RR <- (1+1/m)*(sum((RR-mean(RR))^2)/(m-1))+mean(var.RR);
    degf.RR <- (m-1)*(1 +mean(var.RR)/((1+1/m)*(sum((RR-mean(RR))^2)/(m-1))))*(1+mean(var.RR)/((1+1/m)*(sum((RR-mean(RR))^2)/(m-1))))
    RR.CI <- c(exp(mean(RR)-qt(0.975,df=degf.RR)*sqrt(T.RR)),exp(mean(RR)+qt(0.975,df=degf.RR)*sqrt(T.RR)))
    p.value <- 2*(1-pt(abs(mean(RR)/sqrt(T.RR)),df=degf.RR))
    round(c(exp(mean(RR)),RR.CI,p.value),2)
  }
     
  for (j in 2:3){
     for(u in 1:3){
      tmp <- cumi2.mi(dat.nomiss,10,j,u)
      tab[j+(u-1)*3+(k-1)*9,] <- c(tmp,1-tmp[c(1,3,2)])
     }
  }
}

tab[1:6*3-2,1] <- 1.00; tab[6,] <- NA
tab

#
# End of script.
#


