#
# R syntax to reproduce information for Figure 1 (HK part) from:
#
# Cowling BJ, Ip DKM, Fang VJ, et al.
# Modes of transmission of influenza B virus in households
# PLoS ONE, 2014 (in press).
#
# Last updated by ang VJ and Cowling BJ.
# August 20, 2014
#


source("http://www.hku.hk/bcowling/influenza/NPI_PLoSOne_motB_scripts/dataframe.r")

# HK

# hazard from raw data
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

# plot
windows(width=8,height=4)
layout(matrix(1:2,ncol=2,byrow=T))
par(mar=c(3.5,4,2.5,1))

###
plot(NA,xlim=c(0,8),ylim=c(0,0.16),xlab="",ylab="Cumulative hazardb",main="",axes=F)
lines(hkhazardb$t,cumsum(hkhazardb$h0c),lty=1)
lines(hkhazardb$t,cumsum(hkhazardb$h0h),lty=2)
lines(hkhazardb$t,cumsum(hkhazardb$h0m),lty=3)
axis(1,at=0:4*2);axis(2,at=0:4*0.04,las=1)
mtext("Confirmed infection without fever plus cough",side=3,line=-1,cex=0.65)
mtext("Hong Kong",side=3,line=1,font=2,at=0)

plot(NA,xlim=c(0,8),ylim=c(0,0.16),xlab="",ylab="",main="",axes=F)
lines(hkhazardb$t,cumsum(hkhazardb$h1c),lty=1)
lines(hkhazardb$t,cumsum(hkhazardb$h1h),lty=2)
lines(hkhazardb$t,cumsum(hkhazardb$h1m),lty=3)
axis(1,at=0:4*2);axis(2,at=0:4*0.04,las=1)
legend(4,0.15,legend=c("Control","Hand hygiene","Mask + HH"),lty=1:3,bty="n",border=NA,cex=0.9)
mtext("Confirmed infection with fever plus cough",side=3,line=-1,cex=0.65)


