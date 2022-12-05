#
# R syntax to reproduce information for Figure 2 from:
#
# Lau LLH, Cowling BJ, Fang VJ, Chan KH, Lau EHY, et al.
# Viral shedding and clinical illness in naturally acquired influenza virus infections
# JID, 2010 (in press).
#
# Last updated by Vicky Fang and Lincoln Lau
# Jan 30, 2010

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"
source("http://www.hku.hk/bcowling/influenza/Vshed_scripts/JID_dataframe.r")

hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
demog <- read.csv(paste(dir, "adherence_m.csv", sep=""))
clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
geometric <- function(x)exp(mean(log(x)))
demog$adult <- 1*(demog$age>15)
demog$adult[demog$hhID==148&demog$member==3] <- 0

plotdata <- merge(plotdata, clinic[c(2,13)], all.x=T)

# qPCR
pcr <- merge(plotdata,hchar[c(1,7:9)],all.x=T)
pcr <- pcr[!is.na(pcr$v2_day),]
v1_day <- pcr$v1_day
pcr$v1_day <- pcr$v1_day-pcr$any_onset-v1_day
pcr$v2_day <- pcr$v2_day-pcr$any_onset-v1_day
pcr$v3_day <- pcr$v3_day-pcr$any_onset-v1_day
temp1 <- pcr[c(1:2,11,3,7,10)]
temp2 <- pcr[c(1:2,12,4,8,10)]
temp3 <- pcr[c(1:2,13,5,9,10)]
names(temp1) <- names(temp2) <- names(temp3) <- c("hhID","member","day","pcr","TCID","QVres")
pcr2 <- rbind(temp1,temp2,temp3)
pcr2 <- merge(pcr2,demog[c("hhID","member","adult")],by=c("hhID","member"),all.x=T)
pcr2$pcr[!is.na(pcr2$pcr)&pcr2$pcr<450] <- 450
pcr2 <- pcr2[!is.na(pcr2$day),]
pcr2 <- pcr2[!is.na(pcr2$day),]
pcr2$logpcr <- log10(pcr2$pcr)


###                                             ###
###         PLOT vs Quantitative Culture        ###
###                                             ###

# FLU A
windows(height=8,width=14)
layout(matrix(1:8, ncol=4, byrow=TRUE))
par(mar=c(1.7,2,2,1), oma=c(2,2,1,0))


# Day -1 : Day 0
for (i in 1:7){
    plot(pcr2$TCID[pcr2$day == (i-2) & pcr2$QVres == 1 & pcr2$adult == 1], pcr2$logpcr[pcr2$day == (i-2) & pcr2$QVres == 1 & pcr2$adult == 1]
         ,type ='p', axes=FALSE, xlim=c(0,6), ylim=c(2,10), bty='n', cex=1, xlab="", ylab="")
    points(pcr2$TCID[pcr2$day == (i-2) & pcr2$QVres == 1 & pcr2$adult == 0], pcr2$logpcr[pcr2$day == (i-2) & pcr2$QVres == 1 & pcr2$adult == 0]
           , cex=0.8, pch=3)
    axis(1, pos=2, lwd=1, font=1, at=c(0:6), label=c(0:6), mgp=c(1.5,0.6,0))
    axis(2, pos=0, lwd=1, font=1, at=c(2:10),
        label=c(expression(10^2),expression(10^3),expression(10^4),
        expression(10^5),expression(10^6),expression(10^7),expression(10^8),expression(10^9),expression(10^10)),cex.axis=1,las=1, mgp=c(1.5,0.7,0))
    text(0.75, 10, paste("Day",i-2), cex=1.2, font=2)
    
    if(length(na.exclude(pcr2$TCID[pcr2$day == (i-2) & pcr2$QVres == 1]))>3){
      lm <- lm(pcr2$logpcr[pcr2$day == (i-2) & pcr2$QVres == 1] ~ pcr2$TCID[pcr2$day == (i-2) & pcr2$QVres == 1])
      x1 <- range(pcr2$TCID[pcr2$day == (i-2) & pcr2$QVres == 1],na.rm=T)[1] - 0.1
      x2 <- range(pcr2$TCID[pcr2$day == (i-2) & pcr2$QVres == 1],na.rm=T)[2] + 0.1
      y1 <- summary(lm)[[4]][1] + x1 * summary(lm)[[4]][2]
      y2 <- summary(lm)[[4]][1] + x2 * summary(lm)[[4]][2]
      lines(c(x1,x2),c(y1,y2), lwd=1.15)
      mtext(paste("slope =",round(summary(lm)[[4]][2],2)), side=1, line=-3, cex=0.6, at=5)
      mtext(paste("p-value =",round(summary(lm)[[4]][8],2)), side=1, line=-2.2, cex=0.6, at=5)
    }
}

plot(0,0, type ='p', axes=FALSE, xlim=c(34.0, 41), ylim=c(2,10), bty='n', cex=1, xlab="", ylab="")
legend (36,6.5, c("child","adult"), pch=c(3,1), col=c("black","black"), cex=1.0, bty="n")

mtext("Viral load in NTS (copies/ml)", WEST <- 2, cex=0.7, line=0.2, at=0.75, font=2, outer=TRUE)
mtext("Viral load in NTS (copies/ml)", WEST <- 2, cex=0.7, line=0.2, at=0.25, font=2, outer=TRUE)
mtext(expression(bold(paste(log[10]," ",TCID[50], sep=""))), SOUTH <- 1, cex=0.7, line=-0.2, at=0.13, font=2, outer=TRUE)
mtext(expression(bold(paste(log[10]," ",TCID[50], sep=""))), SOUTH <- 1, cex=0.7, line=-0.2, at=0.38, font=2, outer=TRUE)
mtext(expression(bold(paste(log[10]," ",TCID[50], sep=""))), SOUTH <- 1, cex=0.7, line=-0.2, at=0.63, font=2, outer=TRUE)
mtext(expression(bold(paste(log[10]," ",TCID[50], sep=""))), SOUTH <- 1, cex=0.7, line=-23.6, at=0.88, font=2, outer=TRUE)

#
# End of script
#

