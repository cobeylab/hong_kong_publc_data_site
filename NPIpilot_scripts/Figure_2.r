#
# Script to reproduce information in Figure 2 (onset/randomization/intervention delays) from:
#
# Cowling BJ, Fung ROP, Cheng KY, Fang VJ, Chan KH, Seto WH, et al.
# Preliminary findings of a randomized trial of non-pharmaceutical
# interventions to prevent influenza transmission in households.
# PLoS ONE, 2008; 3(5):e2101.
# http://dx.doi.org/10.1371/journal.pone.0002101
#
# Last updated by Vicky Fang and Ben Cowling
# January 5, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

qv <- read.csv(paste(dir, "clinicdat_h.csv", sep=""), header=TRUE)
dat <- read.csv(paste(dir, "hchar_h.csv", sep=""), header=TRUE)
dat <- merge(dat,qv[,c(2,25)],by="hhID",all.x=TRUE)
dat$onsettime[dat$onsettime!=5] <- floor(dat$onsettime[dat$onsettime!=5]/2)
dat$onsettime[dat$onsettime==5] <- 3
names(dat)[9] <- "clinic_day"

# Function for extracting the delays for plotting
delayday <- function(dataset){
  temp1 <- as.matrix(dataset[dataset==0])
  d0 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==1])
  d1 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==2])
  d2 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==3])
  d3 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset==4])
  d4 <- dim(na.omit(temp1))[1]
  result <- c(d0,d1,d2,d3,d4)
  return (result) 
 }
#end

# For vertically-aligned plots
windows(width=8, height=12)
layout(matrix(1:3, nrow=3))
par(mar=c(5,5,2,2))


#
# 1st part - horizonal barplot of time from appearance of first symptoms to subject recruitment

da <- delayday(dat$clinic_day)
names(da)<-c(0,1,2,3,4)

barplot(rev(da), width = 1, horiz = TRUE, xlim=(c(0,80)), ylim=(c(0,5)), 
 xlab="", axisnames=FALSE, axes= FALSE,
 ylab="Delay (days)", cex.lab=1.5, cex.axis=1.5,
 space=0, las=1)
title("(a)",cex.main=1.5)
axis(1, at=0:8*10, 
     labels=c("0","10","20","30","40","50","60","70","80"), cex.axis=1.5)
axis(2, at=c(0.5,1.5,2.5,3.5,4.5), labels=c("4","3","2","1","0"), cex.axis=1.5,las=1)


#
# 2nd part - horizonal barplot of time from clinic recruitment to first home visit 

between <-  dat$v1_day
da <- delayday(between)
names(da)<- c(0,1,2,3,4)

barplot(rev(da), , width = 1, space=0, xlim=c(0,80), ylim =c(0,5), horiz = TRUE,
  xlab="", ylab="Delay (days)",las=1, axisnames=FALSE, axes= FALSE,
  cex.lab=1.5, cex.axis=1.5)
title("(b)",cex.main=1.5)
axis(1, at=0:8*10, 
     labels=c("0","10","20","30","40","50","60","70","80"), cex.axis=1.5)
axis(2, at=c(0.5,1.5,2.5,3.5,4.5), labels=c("4","3","2","1","0"), cex.axis=1.5,las=1)


#
# 3rd part - horizonal barplot of time from appearance of first symptoms to first home visit 

da <- delayday(dat$v1_day+dat$clinic_day)

names(da)<- c(0,1,2,3,4)
barplot(rev(da), , width = 1, space=0, xlim=c(0,80), ylim =c(0,5), horiz = TRUE,
  axisnames=FALSE, axes= FALSE, 
  xlab="Number of index cases", cex.axis=1.5, ylab="Delay (days)",las=1, cex.lab=1.5 )
title("(c)",cex.main=1.5)
axis(1, at=0:8*10,
     labels=c("0","10","20","30","40","50","60","70","80"), cex.axis=1.5)
axis(2, at=c(0.5,1.5,2.5,3.5,4.5), labels=c("4","3","2","1","0"), cex.axis=1.5,las=1)


# End of script
