#
# Script to reproduce information in Figure 1 (B,D,F,H) from:
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

  vl <- merge(vl,baseflu[baseflu$member==0,c("hhID","age","male")],all.x=T)
  vl$child <- 1*(vl$age<=15)
  vl <- merge(vl,index[c("hhID","tamiflu")],all.x=T)
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

#
#  plot (2x2) HK NPI (flu A & B), logVL & residuals
#

color <- c(gray(0.2),"black")
cex.pt <- c(0.9,1.0)

plot.pcr <- function(datalist,morder){
  par(mar=c(4,4,1,1))
  plot(NA,xlim=c(0,16.8),ylim=c(0,10),axes=FALSE,xlab="",ylab="",main="")
  points(jitter(datalist$pcr$day),jitter(datalist$pcr$logvl),pch=datalist$pcr$ppch+1,col=color[datalist$pcr$ppch/2+1],cex=cex.pt[datalist$pcr$ppch/2+1])
  axis(1,at=0:4*3)
  axis(2,at=0:5*2,las=1)
  mtext("Day from illness onset",side=1,line=2.5,at=6,cex=0.8)
  mtext(expression(paste(log[10],"VL",sep="")),side=2,line=2.5,cex=0.8)
  mtext(morder,side=3,at=0,font=2)
}

plot.pred <- function(datalist,morder,detect.limit){
  datalist$pred <- datalist$pred[datalist$pred$logvl>log10(detect.limit/2),]
  par(mar=c(4,4,1,1))
  plot(NA,xlim=c(3,10),ylim=c(-4,4),axes=FALSE,xlab="",ylab="")
  curve(0*x,col="gray",lty=2,to=8,add=TRUE)
  points(datalist$pred$predict_logvl,datalist$pred$logvl-datalist$pred$predict_logvl,pch=3)
  axis(1,at=3:8)
  axis(2,at=0:4*2-4,las=1)
  mtext(expression(paste("Predicted ",log[10],"VL",sep="")),side=1,line=2.5,at=5.5,cex=0.8)
  mtext(expression(paste("Observed ",log[10],"VL - Predicted ",log[10],"VL",sep="")),side=2,line=2.5,cex=0.8)
  mtext(morder,side=3,at=3,font=2)
  dr <- cut(datalist$pred$logvl-datalist$pred$predict_logvl,breaks=0:8-4)
  dvd <- max(table(dr))/1.6
  for (i in 1:length(table(dr))){
  polygon(c(0,table(dr)[[i]],table(dr)[[i]],0)/dvd+8.3,c(-4,-4,-3,-3)+1*(i-1),lwd=0.8)
  }
}

windows(width=9,height=8)
layout(matrix(c(1,2,3:6), ncol=2, byrow=T),heights = c(1.5,5,5))

par(mar=c(2,2,2,2))
plot(0, type="n", axes=FALSE, xlab="", ylab="",xlim=c(0,10),ylim=c(0,10))
mtext("Influenza A",cex=1.5,line=-1)
lines(c(0,0,10,10),c(0,3,3,0))

par(mar=c(2,2,2,2))
plot(0, type="n", axes=FALSE, xlab="", ylab="",xlim=c(0,10),ylim=c(0,10))
mtext("Influenza B",cex=1.5,line=-1)
lines(c(0,0,10,10),c(0,3,3,0))

set.seed(12345)
plot.pcr(npi.fluA,"B")
plot.pcr(npi.fluB,"D")

plot.pred(npi.fluA,"F",900)
plot.pred(npi.fluB,"H",900)

# End of script.

