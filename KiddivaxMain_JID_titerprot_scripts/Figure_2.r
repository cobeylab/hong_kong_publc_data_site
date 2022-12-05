#
# R syntax to reproduce Figure 2 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_3.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_4.r")

# Format data~~~~~~~
dat.vic <-Pred.B.titer
dat.vic$E.risk <- exp(dat.vic$postvax.b * vic.vic$coef[1] )
denom <- exp(log2(5)* vic.vic$coef[1] )

mrow<-361

#plot                          
windows(width=10.5,height=4.3)
par (oma=c(0,2,0,0),mar=c(4,3,1,1))

layout(matrix(1:2, nrow=1, byrow=T))

# B TIV 6:8

dat <-dat.vic[ dat.vic$age%in%5:8 & dat.vic$TIV==1,]
dat$protect <- 1-(dat$E.risk/denom)

decile <- matrix(rep(NA,5*mrow), nrow=5)
for ( i in 1:mrow) {
                          decile[,i] <- quantile(dat$protect[dat$t==i-1], c((0:4)* 0.25), na.rm=T)
                          }

plot(0,0,xlim=c(0,30*9), ylim=c(0,1), type="n", axes=F, xlab="",ylab="")
axis(1, at=(0:9)*30, labels=0:9)
axis(2, at=(0:4)*0.25, labels=(0:4)*25, las=1)
mtext("Protection (%)", side=2, line=3.5)
mtext("Months from vaccination", side=1, line=3)
mtext("6-8y", side=3, line=0)
mtext("A)", side=2, line=3.5, at=1.07, las=1)

for ( i in 31:(9*30+1)) {
                          segments(i, decile[1,i], i, decile[5,i], col=rgb(0.2,0.2,0.2, alpha=0.1))
                          segments(i, decile[2,i], i, decile[4,i], col=rgb(0.2,0.2,0.2, alpha=0.1))
                          points(i, decile[3,i], pch=16, col=rgb(0.2,0.2,0.2, alpha=0.2))
                          }
                          
text((9*30), decile[2,(9*30)]-0, "25th pct", adj=1, cex=0.8)
text((9*30), decile[3,(9*30)]-0.03, "50th pct", adj=1, cex=0.8)
text((9*30), decile[4,(9*30)]-0.03, "75th pct", adj=1, cex=0.8)
text((9*30), decile[5,(9*30)]-0.03, "100th pct", adj=1, cex=0.8)

# B TIV 9:18

dat <-dat.vic[dat.vic$age%in%9:18 & dat.vic$TIV==1,]
dat$protect <- 1-(dat$E.risk/denom)

decile <- matrix(rep(NA,5*mrow), nrow=5)
for ( i in 31:mrow) {
                          decile[,i] <- quantile(dat$protect[dat$t==i-1], c((0:4)* 0.25), na.rm=T)
                          }

plot(0,0,xlim=c(0,30*9), ylim=c(0,1), type="n", axes=F, xlab="",ylab="")
axis(1, at=(0:9)*30, labels=0:9)
axis(2, at=(0:4)*0.25, labels=(0:4)*25, las=1)
#mtext("Protection (%)", side=2, line=3.5)
mtext("Months from vaccination", side=1, line=3)
mtext("9-17y", side=3, line=0)
mtext("B)", side=2, line=3.5, at=1.07, las=1)

for ( i in 31:(9*30+1)) { segments(i, decile[1,i], i, decile[5,i], col=rgb(0.2,0.2,0.2, alpha=0.1))
                          segments(i, decile[2,i], i, decile[4,i], col=rgb(0.2,0.2,0.2, alpha=0.1))
                          points(i, decile[3,i], pch=16, col=rgb(0.2,0.2,0.2, alpha=0.2))
                          }
                          
text((9*30), decile[2,(9*30)]-0.03, "25th pct", adj=1, cex=0.8)
text((9*30), decile[3,(9*30)]-0.03, "50th pct", adj=1, cex=0.8)
text((9*30), decile[4,(9*30)]-0.03, "75th pct", adj=1, cex=0.8)
text((9*30), decile[5,(9*30)]-0.03, "100th pct", adj=1, cex=0.8)


#
# End of script.
#


