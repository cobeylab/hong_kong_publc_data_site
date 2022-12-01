#
# Script to reproduce results in main text from:
#
# Ip DKM, Schutten M, Fang VJ, Fung ROP, et al.
# Validation of self-swab for virologic confirmation of influenza virus infections
# in a community setting.
# Journal of Infectious Diseases, 2012; 205(4): 631-634
#
# Last updated by Vicky Fang and Ben Cowling
# Feb 9, 2012

require(nlme)

source("http://www.hku.hk/bcowling/influenza/NPImain_scripts/Analyzed_hh.r")
dir <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""))
clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
av <- read.csv(paste(dir, "antiviral_m.csv", sep=""))

# Reshape qPCR results
mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
hc2 <- data.frame(hhID = baseflu$hhID, member = baseflu$member)
hc.temp <- reshape(hc[1:4], timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="qPCR")
hc2 <- merge(hc2,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hc2) <- c("hhID","member","V1","V2","V3")
hc2$exd_index <- hculture$exd_index

# construct the data frame
index <- hc2[hc2$member==0&hc2$exd_index==0,]
hchar$v2_day[hchar$hhID==12] <- 3 # data error
index <- merge(index,hchar[c(1,7:9)],all.x=TRUE,by="hhID")
index <- merge(index,clinic[c("hhID","QVres")],all.x=T)
av <- av[av$member==0&av$av=="tamiflu",]; av$tamiflu <- 1
index <- merge(index,av[c("hhID","tamiflu")],all.x=T)
index$tamiflu[is.na(index$tamiflu)] <- 0
index <- merge(index,baseflu[baseflu$member==0,c("hhID","age","male")],all.x=T)

for(i in 1:2){

  index.flu <- index[index$QVres==i,]

  # reshape the dataset
  vl1 <- index.flu[c(1,3,7)]
  names(vl1) <- c("hhID","vl","day")
  vl1$visit <- 1
  vl2 <- index.flu[c(1,4,8)]
  names(vl2) <- c("hhID","vl","day")
  vl2$visit <- 2
  vl3 <- index.flu[c(1,5,9)]
  names(vl3) <- c("hhID","vl","day")
  vl3$visit <- 3

  vl <- rbind(vl1,vl2,vl3)
  vl$detect.mark <- 1*(!is.na(vl$vl)&vl$vl<=900&vl$vl>0)
  vl$vl[!is.na(vl$vl)&vl$vl<=900] <- 450
  vl$logvl <- log10(vl$vl)

  vl <- merge(vl,index[c("hhID","age","male","tamiflu")],all.x=T)
  vl$child <- 1*(vl$age<=15)
  vl$daysq <- vl$day^2
  vl <- vl[order(vl$hhID,vl$visit),]
  nvl <- length(unique(vl$hhID))
  vl$ppch <- rep(c(0,2,0),nvl)
  vl.re <- vl[-(1:nvl*3-1),]
  vl.re.fit <- vl.re[!is.na(vl.re$day)&!is.na(vl.re$logvl),]

  # lme
  fit_lme <- lme(logvl~day+child+child*day+tamiflu+tamiflu*day, data=vl.re.fit, random=~1|hhID)
  pd2 <- data.frame(hhID=vl$hhID[1:nvl*3-1],day=vl$day[1:nvl*3-1],daysq=vl$daysq[1:nvl*3-1],child=vl$child[1:nvl*3-1],
                    tamiflu=vl$tamiflu[1:nvl*3-1])
  pd2 <- pd2[!is.na(pd2$day),]
  pd2$predict_logvl <- predict(fit_lme,pd2)
  pd2 <- merge(pd2,vl[c("hhID","day","logvl")],by=c("hhID","day"),all.x=T)

  if(i==1) npi.fluA <- list(pcr=vl,fit=fit_lme,pred=pd2)
  if(i==2) npi.fluB <- list(pcr=vl,fit=fit_lme,pred=pd2)
}

###

# Demographics
table(index$QVres)  # number of subjects with flu A and flu B
mean(index$age<=15); range(index$age)   # age distribution
mean(index$male==0); mean(index$tamiflu)  # female, tamiflu

# Swabs collected in visit 2 with detectable virus
tmp <- index; tmp$mark <- 1*(tmp$V2>0)
table(tmp$QVres,tmp$mark)

# Compare between observed and predicted viral load in visit 2

# Flu A
dataA <- npi.fluA$pred
dataA <- dataA[dataA$logvl>log10(450)&!is.na(dataA$logvl),]
dA <- dataA$logvl-dataA$predict_logvl
round(c(mean(dA),mean(dA)-qt(0.975,length(dA)-1)*sd(dA)/sqrt(length(dA)),mean(dA)+qt(0.975,length(dA)-1)*sd(dA)/sqrt(length(dA))),2)   #log_10 copies/mL

# Flu B
dataB <- npi.fluB$pred
dataB <- dataB[dataB$logvl>log10(450)&!is.na(dataB$logvl),]
dB <- dataB$logvl-dataB$predict_logvl
round(c(mean(dB),mean(dB)-qt(0.975,length(dB)-1)*sd(dB)/sqrt(length(dB)),mean(dB)+qt(0.975,length(dB)-1)*sd(dB)/sqrt(length(dB))),2)   #log_10 copies/mL

# End of script.



