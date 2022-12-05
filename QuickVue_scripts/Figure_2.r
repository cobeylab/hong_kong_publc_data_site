#
# R syntax to reproduce information in Figure 2 (viral load vs culture/QV results) from:
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


#seperate into 3 groups, quickvue culture --, A-+, B-+, A++, B++
group0a <- dat[dat$QVres == 1 & dat$culture==0,]
group0b <- dat[dat$QVres == 2 & dat$culture==0,]
group1 <- dat[dat$QVres == 3 & dat$culture==0,]
group2 <- dat[dat$QVres == 3 & dat$culture!=0,]
group3 <- dat[dat$QVres != 3 & dat$culture!=0,]

group2a <- group2[group2$culture =="A",]
group2b <- group2[group2$culture =="B",]
group3a <- group3[group3$culture =="A",]
group3b <- group3[group3$culture =="B",]

groupa <-rbind (group2a,group3a)
groupb <-rbind (group2b,group3b)

lower.lim <- log(900,10)

# Plot 

windows(width=10, height=7)
par(mar=c(5.5,4,1,2))
plot(0, type="n", axes=FALSE, xlim=c(-1,11), ylim=c(2,10),xlab="", ylab="",las=1)
axis(1,c(0.5, 3, 4.5, 7, 8.5),labels=rep("", 5), cex.axis=1.5)         
axis(2, pos=0, at=2:10, labels= c(expression(10^2),"",expression(10^4),"",expression(10^6),"",expression(10^8),"",expression(10^10)
), las=1)
mtext("Viral load (copies/ml)",side=2,line=1)

# for more jitter in group 1
tmp <- group1$logqPCR
for (i in 1:15){
  inc_fac <- 0.1
  tmp[i] <- tmp[i] + inc_fac
  inc_fac <- inc_fac +0.5
}

points(jitter(rep(0.5, length(group1$logqPCR)), factor=10), pch=16,tmp)
points(jitter(rep(3, length(group2a$logqPCR)), factor=1), pch=16,group2a$logqPCR)
points(jitter(rep(4.5, length(group3a$logqPCR)), factor=1), pch=16,group3a$logqPCR)
points(jitter(rep(7, length(group2b$logqPCR)), factor=1), pch=16,group2b$logqPCR)
points(jitter(rep(8.5, length(group3b$logqPCR)), factor=1), pch=16,group3b$logqPCR)

mtext("QuickVue result:",side=1,at=-1.5, line=1.0)
mtext("Culture result:",side=1,at=-1.7, line=2.5)
mtext("Neg" ,side=1,at=0.5, line=1)
mtext("Neg" ,side=1,at=3, line=1)
mtext("Pos A" ,side=1,at=4.5, line=1)
mtext("Neg" ,side=1,at=7, line=1)
mtext("Pos B",side=1,at=8.5, line=1)

mtext("Negative",  side=1, at=2-1.5, line=2.5)
mtext("Influenza A", side=1, at=6.75-marg, line=2.5)
mtext("Influenza B", side=1, at=10.75-marg, line=2.5)

mtext("Sample size: ",side=1,at=-1.7, line=4)
mtext(dim(group1)[1] ,side=1,at=0.5, line=4)
mtext(dim(group2a)[1] ,side=1,at=3, line=4)
mtext(dim(group3a)[1] ,side=1,at=4.5, line=4)
mtext(dim(group2b)[1] ,side=1,at=7, line=4)
mtext(dim(group3b)[1],side=1,at=8.5, line=4)

lines(c(0,15),c(lower.lim,lower.lim),col=gray(0.5))
mtext("RT-qPCR
detection limit", side=1, at=13.5-marg, line=-6)

#for median and IQR of each group

mediangp <-vector(length=4)
lowiqrgp<-vector(length=4)
upiqrgp <-vector(length=4)
mediangp[1] <- median(group2a$logqPCR)
mediangp[2] <- median(group2b$logqPCR)
mediangp[3] <- median(group3a$logqPCR)
mediangp[4] <- median(group3b$logqPCR)

lowiqrgp[1] <-quantile (group2a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[1] <-quantile (group2a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[2] <-quantile (group2b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[2] <-quantile (group2b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[3] <-quantile (group3a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[3] <-quantile (group3a$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]
lowiqrgp[4] <-quantile (group3b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[2]
upiqrgp[4] <-quantile (group3b$logqPCR, c(0, 0.25, 0.5, 0.75, 1))[4]

lines(c(3.35,3.45), c(mediangp[1],mediangp[1]),col="red")
lines(c(3.4,3.4),c(lowiqrgp[1],upiqrgp[1]),col="red")

lines(c(4.85,4.95), c(mediangp[3],mediangp[3]),col="red")
lines(c(4.9,4.9),c(lowiqrgp[3],upiqrgp[3]),col="red")

lines(c(7.35,7.45), c(mediangp[2],mediangp[2]),col="red")
lines(c(7.4,7.4),c(lowiqrgp[2],upiqrgp[2]),col="red")

lines(c(8.85,8.95), c(mediangp[4],mediangp[4]),col="red")
lines(c(8.9,8.9),c(lowiqrgp[4],upiqrgp[4]),col="red")

############################################################
#090119 add number of qPCR -ve cases in each group
negsize <- vector(length=7)
negsize[1] <- sum(group1$qPCR==0)
negsize[2] <- sum(group0a$qPCR==0)
negsize[3] <- sum(group0b$qPCR==0)
negsize[4] <- sum(group2a$qPCR==0)
negsize[5] <- sum(group3a$qPCR==0)
negsize[6] <- sum(group2b$qPCR==0)
negsize[7] <- sum(group3b$qPCR==0)


# END OF SCRIPTS

