#
# Script to reproduce information in Figure S1 (study recruitment vs local influenza activity) from:
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

plot.data <- read.csv(paste(dir, "qvchpqmrate.csv", sep=""), header=TRUE)

months <- rep(0,40)
months[1] <- 1
for(i in 2:40){
    if(plot.data$month[i]!=plot.data$month[i-1]) months[i] <- 1
}
labs <- c(sort(unique(months*1:40))[-1], 40)

# For vertically aligned plots
windows(width=8, height=7)
layout(matrix(1:3, ncol=1))

#
# 1st part - NPI study recruitment data by QV result

par(mar=c(3,5,3,1))
barplot(height=t(as.matrix(plot.data[1:40,9:10])), space=0,
 axes=FALSE, axisnames=FALSE, col=c(gray(0.2),0),
 xlim=c(0,40), ylim=c(0,12), xlab=" ", ylab="Recruitment rate",
 cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,3,6,9,12) ,  labels=c("0","3","6","9","12"), las=1, cex.axis=1.2)
legend(33, 11, cex=1,
  legend=c("QV -ve", "QV +ve"), fill=c(0, gray(0.2)), lty=NULL)
mtext("(a)", cex=1)

#
# 2nd part - Centre for Health Protection ILI sentinel surveillance data, GOPC and GP rates

par(mar=c(3,5,3,1))
plot(as.matrix(plot.data[1:40,18]), axes=FALSE, type="l",
 xlim=c(0,40), ylim=c(0,80), xlab=" ", ylab="ILI rate", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:4*20, labels=c("0","20","40","60","80"), las=1, cex.axis=1.2)
mtext("(b)", cex=1)

#
# 3rd part - Queen Mary Hospital laboratory surveillance data

par(mar=c(5,5,3,1))
barplot(height=t(as.matrix(plot.data[1:40,15:16])), space=0, axes=FALSE, axisnames=FALSE,
 xlim=c(0,40), ylim=c(0,0.3), xlab=" ", ylab="Isolation rate", cex.lab=1.5,
 col=c(gray(0.2),0))
legend(32.5, 0.3, cex=1,
 legend=c("Flu B +ve", "Flu A +ve"), fill=c(0, gray(0.2)), lty=NULL)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:3/10, labels=c("0%","10%","20%","30%"), las=1, cex.axis=1.2)

mtext("Jan", side=1, at=2.5, line=1, cex=1)
mtext("Feb", side=1, at=7, line=1, cex=1)
mtext("Mar", side=1, at=11, line=1, cex=1)
mtext("Apr", side=1, at=15.5, line=1, cex=1)
mtext("May", side=1, at=20, line=1, cex=1)
mtext("Jun", side=1, at=24, line=1, cex=1)
mtext("Jul", side=1, at=28.5, line=1, cex=1)
mtext("Aug", side=1, at=33, line=1, cex=1)
mtext("Sep", side=1, at=37, line=1, cex=1)

mtext("(c)", cex=1)

# End of script
