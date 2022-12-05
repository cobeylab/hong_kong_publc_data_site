
#
# R syntax to reproduce information for Figure 1 from:
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

plotdata$T1 <- 0
plotdata <- merge(plotdata, clinic[c(2,13)], all.x=T)
data <- merge(plotdata[c(1,2,6,10)],symp, by=c("hhID", "member"))


# Symptom Plot
for (i in 1:dim(data)[1]){
  data$up[i] <- sum( data$sthroat[i], data$rnose[i], na.rm=TRUE)/2
  data$lo[i] <- sum( data$cough[i], data$phlegm[i], na.rm=TRUE)/2
  data$sy[i] <- sum( data$fever[i], data$headache[i], data$pmuscle[i],na.rm=TRUE)/3
}

data2 <- data[order(data$hhID, data$member),c(1,2,5,3:4,6,15:17)]   # check the definition of symptom score
data2$day <- data2$day-data2$any_onset                                   
data2 <- merge(data2,demog[c("hhID","member","adult")],by=c("hhID","member"),all.x=T)

# qPCR plot
pcr <- merge(plotdata,hchar[c(1,7:9)],all.x=T)
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
pcr2$TCID[!is.na(pcr2$TCID)&pcr2$TCID==0] <- 0.3
pcr2 <- pcr2[!is.na(pcr2$day),]
pcr2 <- pcr2[!is.na(pcr2$day),]
pcr2$logpcr <- log10(pcr2$pcr)

# calculate the mean viral load, TCID
meandata <- data.frame(day=sort(unique(pcr2$day)))
for (i in 1:nrow(meandata)){
   meandata$pcrA[i] <- log10(geometric(pcr2$pcr[pcr2$day==meandata$day[i]&pcr2$QVres==1]))
   meandata$pcrB[i] <- log10(geometric(pcr2$pcr[pcr2$day==meandata$day[i]&pcr2$QVres==2]))
   meandata$TCIDA[i] <- mean(na.exclude(pcr2$TCID[pcr2$day==meandata$day[i]&pcr2$QVres==1]))
   meandata$TCIDB[i] <- mean(na.exclude(pcr2$TCID[pcr2$day==meandata$day[i]&pcr2$QVres==2]))
}

# calculate the mean symptom score, fever
meansymp <- data.frame(day=sort(unique(data2$day)))
data2$bodytemp[data2$bodytemp==0] <- NA
for (i in 1:nrow(meansymp)){
   meansymp$upA[i] <- mean(data2$up[data2$day==meansymp$day[i]&data2$QVres==1])
   meansymp$upB[i] <- mean(data2$up[data2$day==meansymp$day[i]&data2$QVres==2])
   meansymp$loA[i] <- mean(data2$lo[data2$day==meansymp$day[i]&data2$QVres==1])
   meansymp$loB[i] <- mean(data2$lo[data2$day==meansymp$day[i]&data2$QVres==2])
   meansymp$syA[i] <- mean(data2$sy[data2$day==meansymp$day[i]&data2$QVres==1])
   meansymp$syB[i] <- mean(data2$sy[data2$day==meansymp$day[i]&data2$QVres==2])
   meansymp$feverA[i] <- mean(data2$bodytemp[data2$day==meansymp$day[i]&data2$QVres==1],na.rm=T)
   meansymp$feverB[i] <- mean(data2$bodytemp[data2$day==meansymp$day[i]&data2$QVres==2],na.rm=T)
}


###                                     ###
###                 Plot                ###
###                                     ###

windows(height=14,width=12)
layout(matrix(1:8, ncol=2, byrow=TRUE), heights=c(2,2,2,1))
par(mar=c(1.5,3,1,1), oma=c(2,0,0,0))

# Flu A - Viral Shedding
plot(NA, type='l', axes=FALSE, xlim = c(-6,10), ylim=c(2,10),xlab="", ylab= "", font.lab=2, bty='n', cex.lab=1.15,lwd=0, mgp=c(1.5,1,0))
    lines(c(-6,10),c(log10(450),log10(450)),col="grey")
    lines(meandata$day[4:13], meandata$pcrA[4:13], lwd=2)
    axis(1, pos=2, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font =1, at=c(2,4,6,8,10),
        label=c(expression(10^2),expression(10^4),expression(10^6),expression(10^8),expression(10^10)),cex.axis=1.3,las=1, mgp=c(1.5,0.7,0))
    points(jitter(pcr2$day[pcr2$adult == 1 & pcr2$QVres == 1]), pcr2$logpcr[pcr2$adult == 1 & pcr2$QVres == 1], pch=1, cex=0.75)
    points(jitter(pcr2$day[pcr2$adult == 0 & pcr2$QVres == 1]), pcr2$logpcr[pcr2$adult == 0 & pcr2$QVres == 1], pch=3, cex=0.75)
    mtext("Viral load in NTS (copies/ml)", side = 2, line = 1.9, cex=0.7, font=2)
    text(-5,10,"A",font=2, cex=1.2)

# Flu B - Viral Shedding
plot(NA, type='l', axes=FALSE, xlim = c(-6,10), ylim=c(2,10),xlab="", ylab= "", font.lab=2, bty='n', cex.lab=1.3,lwd=0)
    lines(c(-6,10),c(log10(450),log10(450)),col="grey")
    lines(meandata$day[4:14], meandata$pcrB[4:14], lwd=2)
    axis(1, pos=2, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font =1, at=c(2,4,6,8,10),
        label=c(expression(10^2),expression(10^4),expression(10^6),expression(10^8),expression(10^10)),cex.axis=1.3,las=1, mgp=c(1.5,0.7,0))
    points(jitter(pcr2$day[pcr2$adult == 1 & pcr2$QVres == 2]), pcr2$logpcr[pcr2$adult == 1 & pcr2$QVres == 2], pch=1, cex=0.75)
    points(jitter(pcr2$day[pcr2$adult == 0 & pcr2$QVres == 2]), pcr2$logpcr[pcr2$adult == 0 & pcr2$QVres == 2], pch=3, cex=0.75)
    text (8.6,2.87, c("PCR detection limit"), cex =0.65)
    legend (6.4,9.9, c("child","adult"), pch=c(3,1), col=c("black","black"), cex=0.9, bty="n")
    text(-5,10,"B",font=2, cex=1.2)

# Flu A - TCID 50
plot(NA, axes=FALSE, ylim=c(0,6), xlim=c(-6,10), bty='n', xlab="", ylab="")
    lines(c(-6,10),c(0.3,0.3),col="grey") # detection limit
    lines(meandata$day[c(5:13)],meandata$TCIDA[c(5:13)], lwd=2)
    points(jitter(pcr2$day[pcr2$adult == 1 & pcr2$QVres == 1]), pcr2$TCID[pcr2$adult == 1 & pcr2$QVres == 1] , pch=1, cex=0.75)
    points(jitter(pcr2$day[pcr2$adult == 0 & pcr2$QVres == 1]), pcr2$TCID[pcr2$adult == 0 & pcr2$QVres == 1] , pch=3, cex=0.75)
    axis(1, pos=0, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10),, cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font=1, las=1, cex.axis=1.3, mgp=c(1.5,0.6,0))
    mtext(expression(bold(paste(log[10]," ",TCID[50], sep=""))), side = 2, line = 1.8, cex=0.7, font=2)
    text(-5,6,"C",font=2, cex=1.2)

# Flu B - TCID 50
plot(NA, axes=FALSE, ylim=c(0,6), xlim=c(-6,10), bty='n', xlab="", ylab="")
    lines(c(-6,10),c(0.3,0.3),col="grey") # detection limit
    lines(meandata$day[c(4:12)],meandata$TCIDB[c(4:12)], lwd=2)
    points(jitter(pcr2$day[pcr2$adult == 1 & pcr2$QVres == 2]), pcr2$TCID[pcr2$adult == 1 & pcr2$QVres == 2] , pch=1, cex=0.75)
    points(jitter(pcr2$day[pcr2$adult == 0 & pcr2$QVres == 2]), pcr2$TCID[pcr2$adult == 0 & pcr2$QVres == 2] , pch=3, cex=0.75)
    axis(1, pos=0, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10),, cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font=1, las=1, cex.axis=1.3, mgp=c(1.5,0.6,0))
    text(-5,6,"D",font=2, cex=1.2)
    legend (6.4,4.7, c("child","adult"), pch=c(3,1), col=c("black","black"), cex=0.9, bty="n")
    text (8.6,0.425, c("TCID detection limit"), cex =0.65 , bty="n")

# Flu A - Symptoms
plot(meansymp$day[4:14],meansymp$upA[4:14], axes=FALSE, type='l', ylim=c(0,0.8), xlim=c(-6,10), bty='n', xlab="",
         ylab="", font.lab=2, cex.lab=1.15, lwd=2, lty="dashed", mgp=c(1.5,1,0))
    axis(1, pos=0, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font=1,las=1, cex.axis=1.3, mgp=c(1.5,0.7,0))
    lines(meansymp$day[4:14], meansymp$loA[4:14], type='l', lwd=2,lty="dotdash")
    lines(meansymp$day[4:14], meansymp$syA[4:14], type='l', lwd=2)
    mtext("Symptom scores", side = 2, line = 1.9, cex=0.7, font=2)
    text(-5,0.8,"E",font=2, cex=1.2)

# Flu B - Symptoms
plot(meansymp$day[4:14],meansymp$upB[4:14], axes=FALSE, type='l', ylim=c(0,0.8), xlim=c(-6,10), bty='n', xlab="",
         ylab="", font.lab=2, cex.lab=1.15, lwd=2, lty="dashed", mgp=c(1.5,1,0))
    axis(1, pos=0, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font=1,las=1, cex.axis=1.3, mgp=c(1.5,0.7,0))
    lines(meansymp$day[4:14], meansymp$loB[4:14], type='l', lwd=2,lty="dotdash")
    lines(meansymp$day[4:14], meansymp$syB[4:14], type='l', col="black", lwd=2)
    legend (4.8,0.77, c("lower respir.","upper respir.","systemic"), lwd=1.7, lty=c("dotdash","dashed","solid"),
            col=c("black","black","black"), cex=0.9, bty="n")
    text(-5,0.8,"F",font=2, cex=1.2)

# Flu A - Temp
plot(meansymp$day[4:14], meansymp$feverA[4:14], axes=FALSE, type='l', ylim=c(36.35,38), xlim=c(-6,10), bty='n', xlab="",
         ylab="", font.lab=2, cex.lab=1.15, lwd=2, mgp=c(1.7,1,0))
    lines(c(-6,10),c(37.8,37.8),col="grey", lty= "dashed")
    axis(1, pos=36.35, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font=1,las=1, at=c(36.35, 37, 38),
        label=c("", 37, 38), cex.axis=1.3, mgp=c(1.5,0.7,0))
    mtext(expression(bold(paste("Temp (",degree,"C)", sep=""))), side = 2, line=1.8, cex= 0.7)
    text(-5,37.95,"G",font=2, cex=1.2)

# Flu B - Temp
plot(meansymp$day[4:14], meansymp$feverB[4:14], axes=FALSE, type='l', ylim=c(36.35, 38), xlim=c(-6,10), bty='n', xlab="",
         ylab="", font.lab=2, cex.lab=1.15, lwd=2,mgp=c(1.7,1,0))
    lines(c(-6,10),c(37.8,37.8),col="grey", lty= "dashed")
    axis(1, pos=36.35, lwd=1.5, font=1, at=c(-6,-4,-2,0,2,4,6,8,10), label=c(-6,-4,-2,0,2,4,6,8,10), cex.axis=1.3, mgp=c(1.5,0.6,0))
    axis(2, pos=-6, lwd=1.5, font=1,las=1, at=c(36.35, 37, 38),
        label=c("", 37, 38), cex.axis=1.3, mgp=c(1.5,0.7,0))
    legend (7.3,37.92, c(expression(paste("37.8",degree,"C", sep=""))), cex =0.9, bty="n")
    text(-5,37.95,"H",font=2, cex=1.2)

mtext("Time (days since symptom onset)", side = 1, at=0.28, line = 0.2, cex=0.7, font=2, outer=TRUE)
mtext("Time (days since symptom onset)", side = 1, at=0.78, line = 0.2, cex=0.7, font=2, outer=TRUE)

#
# End of script
#

