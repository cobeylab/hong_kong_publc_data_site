#
# Script to reproduce information in Appendix Figure 2 (adherence) from:
#
# Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, et al.
# Facemasks and hand hygiene to prevent influenza transmission 
# in households, a randomized trial.
# Annals of Internal Medicine, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Sep 08, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV3/"
source("http://www.hku.hk/bcowling/influenza/NPImain_scripts/Analyzed_hh.r")
analyzed <- hculture[hculture$member==0,c("hhID","analyzed")]

handday <- read.csv(paste(dir, "adherence_d.csv", sep=""))
handday <- merge(handday,housechar[1:2],by="hhID",all.x=TRUE)
handday$hygiene_day <- rep(NA,nrow(handday))
for (i in 1:nrow(handday)){
      handday$hygiene_day[i] <- min(handday$handrub_day[i],handday$soap_day[i],na.rm=TRUE)
}
handday$hygiene_day[handday$hygiene_day=="Inf"] <- NA
handday <- merge(handday,analyzed,all.x=TRUE)
handday <- handday[handday$analyzed==1,1:9]

#
# Calculate the mean and s.e. of daily frequency of hand rub use for index vs contact (hand group)
#

handuse <- handday[handday$intervention==3,]
tmp <- rep(NA, 10)
output <- data.frame(day=0:9, index.mean=tmp, index.se=tmp, contact.mean=tmp, contact.se=tmp)
for(i in 1:10){
  output$index.mean[i] <- mean(handuse$hygiene_day[!is.na(handuse$hygiene_day) & handuse$member==0 & handuse$day==(i-1)])
  output$index.se[i] <- sqrt(var(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member==0 & handuse$day==(i-1)]) / length(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member==0 & handuse$day==(i-1)]))
  output$contact.mean[i] <- mean(handuse$hygiene_day[!is.na(handuse$hygiene_day) & handuse$member!=0 & handuse$day==(i-1)])
  output$contact.se[i] <- sqrt(var(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member!=0 & handuse$day==(i-1)]) / length(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member!=0 & handuse$day==(i-1)]))
}
output1 <- output[1:8,]

#
# Calculate the mean and s.e. of daily frequency of hand rub use for index vs contact (maskhh group)
#
handuse <- handday[handday$intervention==4,]
tmp <- rep(NA, 10)
output <- data.frame(day=0:9, index.mean=tmp, index.se=tmp, contact.mean=tmp, contact.se=tmp)
for(i in 1:10){
  output$index.mean[i] <- mean(handuse$hygiene_day[!is.na(handuse$hygiene_day) & handuse$member==0 & handuse$day==(i-1)])
  output$index.se[i] <- sqrt(var(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member==0 & handuse$day==(i-1)]) / length(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member==0 & handuse$day==(i-1)]))
  output$contact.mean[i] <- mean(handuse$hygiene_day[!is.na(handuse$hygiene_day) & handuse$member!=0 & handuse$day==(i-1)])
  output$contact.se[i] <- sqrt(var(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member!=0 & handuse$day==(i-1)]) / length(handuse$hygiene_day[!is.na(handuse$hygiene_day) &
    handuse$member!=0 & handuse$day==(i-1)]))
}
output2 <- output[1:8,]

#
# Calculate the mean and s.e. of daily frequency of mask use for index vs contact
#

maskday <- read.csv(paste(dir, "adherence_d.csv", sep=""))
maskday <- merge(maskday,housechar[1:2],by="hhID",all.x=TRUE)
maskuse <- maskday[maskday$intervention==4,]
tmp <- rep(NA, 10)
output <- data.frame(day=0:9, index.mean=tmp, index.se=tmp,contact.mean=tmp, contact.se=tmp)
for(i in 1:10){
  output$index.mean[i] <- mean(maskuse$mask_day[!is.na(maskuse$mask_day) & maskuse$member==0 & maskuse$day==(i-1)])
  output$index.se[i] <- sqrt(var(maskuse$mask_day[!is.na(maskuse$mask_day) &
    maskuse$member==0 & maskuse$day==(i-1)]) / length(maskuse$mask_day[!is.na(maskuse$mask_day) &
    maskuse$member==0 & maskuse$day==(i-1)]))
  output$contact.mean[i] <- mean(maskuse$mask_day[!is.na(maskuse$mask_day) & maskuse$member!=0 & maskuse$day==(i-1)])
  output$contact.se[i] <- sqrt(var(maskuse$mask_day[!is.na(maskuse$mask_day) &
    maskuse$member!=0 & maskuse$day==(i-1)]) / length(maskuse$mask_day[!is.na(maskuse$mask_day) &
    maskuse$member!=0 & maskuse$day==(i-1)]))
}
output3 <- output[1:8,]


windows(width=6, height=12)
layout(matrix(1:3,ncol=1))

# plot daily frequency of hand rub use for index vs contact
par(mar=c(5,6,2,1), xaxs="i", yaxs="i")
plot(0, axes=FALSE, xlim=c(-0.3,7.3), ylim=c(1,4), xlab="", ylab="")

lines(type= "l", output1$day-0.07, 5-output1$index.mean,lty=1,lwd=1.2)
segments(output1$day-0.07, 5-output1$index.mean-1.96*output1$index.se,
  output1$day-0.07, 5-output1$index.mean+1.96*output1$index.se, lwd=1.2)
points(output1$day-0.07, 5-output1$index.mean, pch=16, cex=1.8)

lines(type= "l", output1$day+0.07, 5-output1$contact.mean, lty=2,lwd=1.2)
segments(output1$day+0.07, 5-output1$contact.mean-1.96*output1$contact.se,
  output1$day+0.07, 5-output1$contact.mean+1.96*output1$contact.se,lwd=1.2)
points(output1$day+0.07, 5-output1$contact.mean, pch=17,cex=1.8)

legend(4.5,4.0, c("Index", "Contact"), pch=c(16,17), lty=c(1,2), bty="n")
axis(1, at=0:7)
axis(2, at=1:4, labels=c("Never", "Sometimes", "Often", "Always"), las=1)
mtext("Day",side=1,line=3.5)

# plot daily frequency of hand rub use for index vs contact
par(mar=c(5,6,2,1), xaxs="i", yaxs="i")
plot(0, axes=FALSE, xlim=c(-0.3,7.3), ylim=c(1,4), xlab="", ylab="")

lines(type= "l", output2$day-0.07, 5-output2$index.mean,lty=1,lwd=1.2)
segments(output2$day-0.07, 5-output2$index.mean-1.96*output2$index.se,
  output2$day-0.07, 5-output2$index.mean+1.96*output2$index.se, lwd=1.2)
points(output2$day-0.07, 5-output2$index.mean, pch=16, cex=1.8)

lines(type= "l", output2$day+0.07, 5-output2$contact.mean,lty=2,lwd=1.2)
segments(output2$day+0.07, 5-output2$contact.mean-1.96*output2$contact.se,
  output2$day+0.07, 5-output2$contact.mean+1.96*output2$contact.se,lwd=1.2)
points(output2$day+0.07, 5-output2$contact.mean, pch=17,cex=1.8)

axis(1, at=0:7)
axis(2, at=1:4, labels=c("Never", "Sometimes", "Often", "Always"), las=1)
mtext("Day",side=1,line=3.5)

# plot daily frequency of mask use for index vs contact
par(mar=c(5,6,2,1), xaxs="i", yaxs="i")
plot(0, type="n", axes=FALSE, xlim=c(-0.3,7.3), ylim=c(1,4), xlab="", ylab="")

lines(type= "l", output3$day-0.07, 5-output3$index.mean,lty=1,lwd=1.2)
segments(output3$day-0.07, 5-output3$index.mean-1.96*output3$index.se,
  output3$day-0.07, 5-output3$index.mean+1.96*output3$index.se, lwd=1.2)
points(output3$day-0.07, 5-output3$index.mean, pch=16, cex=1.8)

lines(type= "l", output3$day+0.07, 5-output3$contact.mean,lty=2,lwd=1.2)
segments(output3$day+0.07, 5-output3$contact.mean-1.96*output3$contact.se,
  output3$day+0.07, 5-output3$contact.mean+1.96*output3$contact.se,lwd=1.2)
points(output3$day+0.07, 5-output3$contact.mean, pch=17,cex=1.8)

axis(1, at=0:7)
axis(2, at=1:4, labels=c("Never", "Sometimes", "Often", "Always"), las=1)
mtext("Day",side=1,line=3.5)


# End of script



