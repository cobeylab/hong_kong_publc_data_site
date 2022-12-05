#
# Script to reproduce information in Appendix Figure 1 (surveillance) from:
#
# Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, et al.
# Facemasks and hand hygiene to prevent influenza transmission 
# in households, a randomized trial.
# Annals of Internal Medicine, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Sep 08, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV3/"
plots <- read.csv(paste(dir, "QV_CHP_QM_rate.csv", sep=""))

months <- rep(0,39)
months[1] <- 1
for(i in 2:39){
    if(plots$month[i]!=plots$month[i-1]) months[i] <- 1
}
labs <- c(sort(unique(months*1:39))[-1], 40)

# for plots side by side

windows(width=7, height=8)
layout(matrix(1:3, ncol=1))
par(mar=c(3,6,3,1))

#
# 1st window - QV recruitment data

barplot(height=t(as.matrix(plots[1:39,5:7] )), space=0, axes=FALSE, axisnames=FALSE, col=c(gray(0.2), gray(0.8),0),
        xlim=c(0,39), ylim=c(0,210), xlab=" ", ylab="Participants Recruited, n", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=c(0,40,80,120,160,200) ,  labels=c("0","40","80","120","160","200"), las=1, cex.axis=1.2)
legend(33, 192, cex=1, bty="n",
       legend=c("QV -ve","QV Flu B +ve", "QV Flu A +ve"), fill=c(0, gray(0.8), gray(0.2)), lty=NULL)

#
# 2nd window - chp surveillance data, GOPC and GP rate

par(mar=c(3,6,3,1))

plot(as.matrix(plots[1:39,4]), axes=FALSE, type="l",xlim=c(0,39), ylim=c(0,80), xlab=" ",
     ylab="ILI Rate, n per 1000 consultations", cex.lab=1.5)
axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:4*20, labels=c("0","20","40","60","80"), las=1, cex.axis=1.2)

#
# 3rd window - queen mary Lab surveillance data

par(mar=c(5,6,3,1))

plots$qmfluarate <- plots$qmflua/plots$qmtotal
plots$qmflubrate <- plots$qmflub/plots$qmtotal
barplot(height=t(as.matrix(plots[1:39,11:12])), space=0, axes=FALSE, axisnames=FALSE,
        xlim=c(0,39), ylim=c(0,0.2), xlab=" ", ylab="Isolation Rate, %", cex.lab=1.5,
        col=c(gray(0.2),gray(0.8)))
legend(32.5, 0.2, cex=1, bty="n",
  legend=c("Flu B +ve", "Flu A +ve"), fill=c(gray(0.8),gray(0.2)), lty=NULL)

axis(1, pos=0, at=labs-1, labels=rep(NA, length(labs)), cex.axis=1.2)
axis(2, pos=0, at=0:2/10, labels=c("0","10","20"), las=1, cex.axis=1.2)

mtext("Jan", side=1, at=2.5, line=1, cex=1)
mtext("Feb", side=1, at=7, line=1, cex=1)
mtext("Mar", side=1, at=11, line=1, cex=1)
mtext("Apr", side=1, at=15.5, line=1, cex=1)
mtext("May", side=1, at=20, line=1, cex=1)
mtext("Jun", side=1, at=24, line=1, cex=1)
mtext("Jul", side=1, at=28.5, line=1, cex=1)
mtext("Aug", side=1, at=33, line=1, cex=1)
mtext("Sep", side=1, at=37, line=1, cex=1)

# End of script

