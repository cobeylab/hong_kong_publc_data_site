#
# R syntax to reproduce information for Figure 1 from:
#
# Cowling BJ, Ip DKM, Fang VJ, et al.
# Aerosol transmission is an important mode of influenza A virus spread
# Nature Communications, 2013 (in press).
#
# Last updated by ang VJ and Cowling BJ.
# March 22, 2013
#

source("http://www.hku.hk/bcowling/influenza/NPI_NatureComm_mot_scripts/dataframe.r")

# HK hazard from raw data
hkhazard <- data.frame(t=0:8,h1c=NA,h1h=NA,h1m=NA,h0c=NA,h0h=NA,h0m=NA)
for(m in 1:3){
 if(m==1) temp <- hkdata[hkdata$X_hh==0,]
 else if(m==2) temp <- hkdata[hkdata$X_hh==1&hkdata$X_fm==0,]
 else temp <- hkdata[hkdata$X_fm==1,]
 for(k in 1:2){
   if(k==1) temp$event <- 1*(temp$delta==1&temp$ILI==1)
   else temp$event <- 1*(temp$delta==1&temp$ILI==0)

   for(i in 1:nrow(hkhazard)){
      hkhazard[i,m+(k-1)*3+1] <- sum(temp$T0==i-1&temp$event==1)/sum(temp$T0>=i-1)
   }
 }
}
hkhazard[is.na(hkhazard)] <- 0

# plot (for HK only)
windows(width=5.5*1.5,height=5.5/2)
layout(matrix(1:3,ncol=3,byrow=T))
par(mar=c(4,4,3,1))

plot(NA,xlim=c(0,8),ylim=c(0,0.12),xlab="Day since index case illness onset",ylab="Cumulative hazard",main="",axes=F)
lines(hkhazard$t,cumsum(hkhazard$h1c),lty=1)
lines(hkhazard$t,cumsum(hkhazard$h0c),lty=2)
axis(1,at=0:4*2);axis(2,at=0:4*0.03,las=1)
legend(0,0.11,legend=c("Infection without fever plus cough","Infection with fever plus cough"),lty=2:1,bty="n",border=NA,cex=0.9)
mtext("Control",side=3,line=-1,cex=0.8)

plot(NA,xlim=c(0,8),ylim=c(0,0.12),xlab="Day since index case illness onset",ylab="",main="",axes=F)
lines(hkhazard$t,cumsum(hkhazard$h1h),lty=1)
lines(hkhazard$t,cumsum(hkhazard$h0h),lty=2)
axis(1,at=0:4*2);axis(2,at=0:4*0.03,las=1)
mtext("Hand hygiene",side=3,line=-1,cex=0.8)
mtext("Hong Kong",side=3,line=1.5,font=2)

plot(NA,xlim=c(0,8),ylim=c(0,0.12),xlab="Day since index case illness onset",ylab="",main="",axes=F)
lines(hkhazard$t,cumsum(hkhazard$h1m),lty=1)
lines(hkhazard$t,cumsum(hkhazard$h0m),lty=2)
axis(1,at=0:4*2);axis(2,at=0:4*0.03,las=1)
mtext("Mask plus hand hygiene",side=3,line=-1,cex=0.8)

#
# End of script.
#

