#
# R syntax to reproduce information in Appendix Figure 1 (viral load vs symptom onset time) from:
#
# Cheng CKY, Cowling BJ, Chan KH, Fang VJ, Seto WH, et al.
# Factors affecting QuickVue Influenza A+B rapid test performance in
# the community setting.
# Diagnostic Microbiology and Infectious Disease, 2009; 65(1): 35-41.
#
# Last updated by Calvin Cheng and Vicky Fang
# Sep 08, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV3/"

cdc <- read.csv(paste(dir,"QuickVue_data.csv",sep=""))
source ("http://www.hku.hk/bcowling/influenza/QuickVue_scripts/add_groups.r")

dat <-cdc[!is.na(cdc$qPCR),]

dat$qPCR[dat$qPCR > 0 & dat$qPCR <900] <- 1000
dat$logqPCR <- log(dat$qPCR,10)
dat$logqPCR[dat$logqPCR == -Inf] <- log(70,10)
dat$logqPCR[dat$logqPCR < log(900,10) ] <- log(70,10)

dat.a <- dat[dat$culture =="A",]
dat.b <- dat[dat$culture =="B",]

group1b <- dat.a[dat.a$onsettime ==1 & dat.a$isadult ==0,]
group2b <- dat.a[dat.a$onsettime ==2 & dat.a$isadult ==0,]
group3b <- dat.a[dat.a$onsettime ==3 & dat.a$isadult ==0,]
group4b <- dat.a[dat.a$onsettime ==4 & dat.a$isadult ==0,]
group5b <- dat.a[dat.a$onsettime ==5 & dat.a$isadult ==0,]

group1a <- dat.a[dat.a$onsettime ==1 & dat.a$isadult ==1,]
group2a <- dat.a[dat.a$onsettime ==2 & dat.a$isadult ==1,]
group3a <- dat.a[dat.a$onsettime ==3 & dat.a$isadult ==1,]
group4a <- dat.a[dat.a$onsettime ==4 & dat.a$isadult ==1,]
group5a <- dat.a[dat.a$onsettime ==5 & dat.a$isadult ==1,]


group1bb <- dat.b[dat.b$onsettime ==1 & dat.b$isadult ==0,]
group2bb <- dat.b[dat.b$onsettime ==2 & dat.b$isadult ==0,]
group3bb <- dat.b[dat.b$onsettime ==3 & dat.b$isadult ==0,]
group4bb <- dat.b[dat.b$onsettime ==4 & dat.b$isadult ==0,]
group5bb <- dat.b[dat.b$onsettime ==5 & dat.b$isadult ==0,]

group1ab <- dat.b[dat.b$onsettime ==1 & dat.b$isadult ==1,]
group2ab <- dat.b[dat.b$onsettime ==2 & dat.b$isadult ==1,]
group3ab <- dat.b[dat.b$onsettime ==3 & dat.b$isadult ==1,]
group4ab <- dat.b[dat.b$onsettime ==4 & dat.b$isadult ==1,]
group5ab <- dat.b[dat.b$onsettime ==5 & dat.b$isadult ==1,]


windows(width=22, height=18)
layout(matrix(1:2, nrow=2))

par(mar=c(4,7,1,4))
plot(0, type="n", axes=FALSE, xlim=c(1.4,11.5), ylim=c(2,10),xlab="", ylab="",las=1)
axis(1, c(2,4,6,8,10), labels=rep("", 5))
axis(2, pos=1, at=2:10, labels= c(expression(10^2),"",expression(10^4),"",expression(10^6),"",expression(10^8),"",expression(10^10)
), las=1)

points(jitter(rep(2, length(group1a$logqPCR)), factor=1), group1a$logqPCR, pch=19)
points(jitter(rep(4, length(group2a$logqPCR)), factor=1), group2a$logqPCR, pch=19)
points(jitter(rep(6, length(group3a$logqPCR)), factor=1), group3a$logqPCR, pch=19)
points(jitter(rep(8, length(group4a$logqPCR)), factor=1), group4a$logqPCR, pch=19)
points(jitter(rep(10, length(group5a$logqPCR)), factor=1), group5a$logqPCR, pch=19)

points(jitter(rep(2.25, length(group1b$logqPCR)), factor=1), group1b$logqPCR, pch=1)
points(jitter(rep(4.25, length(group2b$logqPCR)), factor=1), group2b$logqPCR, pch=1)
points(jitter(rep(6.25, length(group3b$logqPCR)), factor=1), group3b$logqPCR, pch=1)
points(jitter(rep(8.25, length(group4b$logqPCR)), factor=1), group4b$logqPCR, pch=1)
points(jitter(rep(10.25, length(group5b$logqPCR)), factor=1), group5b$logqPCR, pch=1)

mtext("Viral load (copies/ml)",side=2,at=6,line=4)
mtext("Symptom onset",side=1,at=0.5, line=1)
mtext("0-12hrs" ,side=1,at=2, line=1)
mtext("13-24hrs" ,side=1,at=4, line=1)
mtext("25-36hrs" ,side=1,at=6, line=1)
mtext("37-48hrs" ,side=1,at=8, line=1)
mtext(">48hrs",side=1,at=10, line=1)

lower.lim <- log(900,10)
lines(c(0,11.5),c(lower.lim,lower.lim),col=gray(0.5))
mtext("RT-qPCR
detection limit", side=1, at=11.5, line=-4)
mtext("(a)",side=3)


#for median and IQR of each group#

mediangp <-vector(length=6)
lowiqrgp<-vector(length=6)
upiqrgp <-vector(length=6)
mediangp[1] <- median(group1a$logqPCR)
mediangp[2] <- median(group1b$logqPCR)
mediangp[3] <- median(group2a$logqPCR)
mediangp[4] <- median(group2b$logqPCR)
mediangp[5] <- median(group3a$logqPCR)
mediangp[6] <- median(group4a$logqPCR)

lowiqrgp[1] <-quantile (group1a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[1] <-quantile (group1a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[2] <-quantile (group1b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[2] <-quantile (group1b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[3] <-quantile (group2a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[3] <-quantile (group2a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[4] <-quantile (group2b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[4] <-quantile (group2b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[5] <-quantile (group3a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[5] <-quantile (group3a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[6] <-quantile (group4a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[6] <-quantile (group4a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]

lines(c(2.55,2.65), c(mediangp[1],mediangp[1]))
lines(c(2.6,2.6),c(lowiqrgp[1],upiqrgp[1]))

lines(c(2.85,2.95), c(mediangp[2],mediangp[2]))
lines(c(2.9,2.9),c(lowiqrgp[2],upiqrgp[2]),lty=2)

lines(c(4.55,4.65), c(mediangp[3],mediangp[3]))
lines(c(4.6,4.6),c(lowiqrgp[3],upiqrgp[3]))

lines(c(4.85,4.95), c(mediangp[4],mediangp[4]))
lines(c(4.9,4.9),c(lowiqrgp[4],upiqrgp[4]),lty=2)

lines(c(6.55,6.65), c(mediangp[5],mediangp[5]))
lines(c(6.6,6.6),c(lowiqrgp[5],upiqrgp[5]))

lines(c(8.55,8.65), c(mediangp[6],mediangp[6]))
lines(c(8.6,8.6),c(lowiqrgp[6],upiqrgp[6]))


######################################################

par(mar=c(4,7,1,4))
plot(0, type="n", axes=FALSE, xlim=c(1.4,11.5), ylim=c(2,10),xlab="", ylab="",las=1)
axis(1, c(2,4,6,8,10), labels=rep("", 5))
axis(2, pos=1, at=2:10, labels= c(expression(10^2),"",expression(10^4),"",expression(10^6),"",expression(10^8),"",expression(10^10)
), las=1)

points(jitter(rep(2, length(group1ab$logqPCR)), factor=1), group1ab$logqPCR, pch=19)
points(jitter(rep(4, length(group2ab$logqPCR)), factor=1), group2ab$logqPCR, pch=19)
points(jitter(rep(6, length(group3ab$logqPCR)), factor=1), group3ab$logqPCR, pch=19)
points(jitter(rep(8, length(group4ab$logqPCR)), factor=1), group4ab$logqPCR, pch=19)
points(jitter(rep(10, length(group5ab$logqPCR)), factor=1), group5ab$logqPCR, pch=19)

points(jitter(rep(2.25, length(group1bb$logqPCR)), factor=1), group1bb$logqPCR, pch=1)
points(jitter(rep(4.25, length(group2bb$logqPCR)), factor=1), group2bb$logqPCR, pch=1)
points(jitter(rep(6.25, length(group3bb$logqPCR)), factor=1), group3bb$logqPCR, pch=1)
points(jitter(rep(8.25, length(group4bb$logqPCR)), factor=1), group4bb$logqPCR, pch=1)
points(jitter(rep(10.25, length(group5bb$logqPCR)), factor=1), group5bb$logqPCR, pch=1)

mtext("Viral load (copies/ml)",side=2,at=6,line=4)
mtext("Symptom onset",side=1,at=0.5, line=1)
mtext("0-12hrs" ,side=1,at=2, line=1)
mtext("13-24hrs" ,side=1,at=4, line=1)
mtext("25-36hrs" ,side=1,at=6, line=1)
mtext("37-48hrs" ,side=1,at=8, line=1)
mtext(">48hrs",side=1,at=10, line=1)

lower.lim <- log(900,10)
lines(c(0,11.5),c(lower.lim,lower.lim),col=gray(0.5))
mtext("RT-qPCR
detection limit", side=1, at=11.5, line=-4)
mtext("(b)",side=3)

#for median and IQR of each group#
mediangp <-vector(length=3)
lowiqrgp<-vector(length=3)
upiqrgp <-vector(length=3)
mediangp[1] <- median(group2ab$logqPCR)
mediangp[2] <- median(group2bb$logqPCR)
mediangp[3] <- median(group4ab$logqPCR)

lowiqrgp[1] <-quantile (group2ab$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[1] <-quantile (group2ab$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[2] <-quantile (group2bb$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[2] <-quantile (group2bb$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[3] <-quantile (group4ab$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[3] <-quantile (group4ab$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]

lines(c(4.55,4.65), c(mediangp[1],mediangp[1]))
lines(c(4.6,4.6),c(lowiqrgp[1],upiqrgp[1]))

lines(c(4.85,4.95), c(mediangp[2],mediangp[2]))
lines(c(4.9,4.9),c(lowiqrgp[2],upiqrgp[2]),lty=2)

lines(c(8.55,8.65), c(mediangp[3],mediangp[3]))
lines(c(8.6,8.6),c(lowiqrgp[3],upiqrgp[3]))


# END OF SCRIPTS
