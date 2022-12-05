#
# R syntax to reproduce information for Figure 1 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Humoral antibody response after receipt of inactivated seasonal
# influenza vaccinations one year apart in children
# PIDJ, 2012.
#
# Last updated by Ng S, Fang VJ and Cowling BJ.
# February 15, 2013
#

source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_PIDJ_antibody_scripts/dataframe.r")

dat <- sero
dat$group <- paste(dat$pilot.intervention,dat$intervention,sep="")

dat.sH1 <- dat[!is.na(dat$prevax.sH1) & !is.na(dat$postvax.sH1),]
dat.sH3 <- dat[!is.na(dat$prevax.sH3) & !is.na(dat$postvax.sH3),]
dat.pH1 <- dat[!is.na(dat$prevax.pH1) & !is.na(dat$postvax.pH1),]
dat.B.Brisbane <- dat[!is.na(dat$prevax.B.Brisbane) & !is.na(dat$postvax.B.Brisbane),]

#values for y axis
y <- c(5,10,20,40,80,160,320,640,1280,2560,5120,10240,20480)
logy <- log(y,1.15)

#color
c1 <- gray(c(rep(0.6,8),rep(0.2,16))); c2 <- "black"

#plot
windows(height=12, width=7)
layout(matrix(1:4, nrow=4, byrow=TRUE))
par(mar=c(2,3,1,0), oma=c(0,1,0,0), mgp=c(3,1.5,0))

#sH1

plot(-2,-2, axes=FALSE, xlim=c(0, 4), ylim=c(0, 80), xlab="",ylab="")
text(0.35+0:3,1.5,"pre",cex=1)
text(0.65+0:3,1.5,"post",cex=1)
segments(0.5,76,1.5,76, lty=1)
segments(0.5,76,0.5,75, lty=1)
segments(1.5,76,1.5,75, lty=1)
rect(0.72, 71,1.29,80,  col = 'white', border="white")
text(1, 76, "(fold rise) p<0.01", cex=1, bg="white")

segments(0.65,70,1.65,70, lty=1)
segments(0.65,69,0.65,70, lty=1)
segments(1.65,70,1.65,69, lty=1)
rect(1.03, 67, 1.27,71,  col = 'white', border="white")
text(1.15, 70, paste("p=", round(wilcox.test(dat.sH1$postvax.sH1[dat.sH1$group=="TIVTIV"], dat.sH1$postvax.sH1[dat.sH1$group=="placeboTIV"])
$p.value,2) , sep=""), cex=1)

segments(0.35,6,1.35,6, lty=1)
segments(0.35,6,0.35,7, lty=1)
segments(1.35,6,1.35,7, lty=1)
rect(0.73, 3, 0.98,8,  col = 'white', border="white")
#text(0.85, 6, paste("p=", round(wilcox.test(dat.sH1$prevax.sH1[dat.sH1$group=="TIVTIV"], dat.sH1$prevax.sH1[dat.sH1$group=="placeboTIV"])
#$p.value,2) , sep=""), cex=1)
text(0.85, 6, "p<0.01", cex=1)

segments(2.5,74,3.5,74, lty=1)
segments(3.5,74,3.5,73, lty=1)
segments(2.5,74,2.5,73, lty=1)
rect(2.72, 69, 3.29,78,  col = 'white', border="white")
text(3, 74, paste("(fold rise) p=", round(wilcox.test(dat.sH1$postvax.sH1[dat.sH1$group=="TIVplacebo"]/dat.sH1$prevax.sH1[dat.sH1$group=="TIVplacebo"], dat.sH1$postvax.sH1[dat.sH1$group=="placeboplacebo"]/dat.sH1$prevax.sH1[dat.sH1$group=="placeboplacebo"])$p.value,2) , sep=""), cex=1)

segments(0.5+0:3,0,0.5+0:3,65, lty=2)
segments(2.65,67,3.65,67, lty=1)
segments(2.65,67,2.65,66, lty=1)
segments(3.65,67,3.65,66, lty=1)
rect(3.03, 65, 3.27,69,  col = 'white', border="white")
text(3.15, 67, paste("p=", round(wilcox.test(dat.sH1$postvax.sH1[dat.sH1$group=="TIVplacebo"], dat.sH1$postvax.sH1[dat.sH1$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")

segments(2.35,6,3.35,6, lty=1)
segments(2.35,6,2.35,7, lty=1)
segments(3.35,6,3.35,7, lty=1)
rect(2.78, 3, 3.03,8,  col = 'white', border="white")
text(2.9, 6, paste("p=", round(wilcox.test(dat.sH1$prevax.sH1[dat.sH1$group=="TIVplacebo"], dat.sH1$prevax.sH1[dat.sH1$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")

#11
points( jitter(rep(0.35, length(dat.sH1$hhID[dat.sH1$group=="TIVTIV"])),factor=7),
        log(dat.sH1$prevax.sH1[dat.sH1$group=="TIVTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="TIVTIV"]])

points( jitter(rep(0.65, length(dat.sH1$hhID[dat.sH1$group=="TIVTIV"])),factor=7/(0.65/0.35)),
log(dat.sH1$postvax.sH1[dat.sH1$group=="TIVTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="TIVTIV"]])

#10
points( jitter(rep(1.35, length(dat.sH1$hhID[dat.sH1$group=="placeboTIV"])),factor=7/(1.35/0.35)),
        log(dat.sH1$prevax.sH1[dat.sH1$group=="placeboTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="placeboTIV"]])

points( jitter(rep(1.65, length(dat.sH1$hhID[dat.sH1$group=="placeboTIV"])),factor=7/(1.65/0.35)),
log(dat.sH1$postvax.sH1[dat.sH1$group=="placeboTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="placeboTIV"]])

#01
points( jitter(rep(2.35, length(dat.sH1$hhID[dat.sH1$group=="TIVplacebo"])),factor=7/(2.35/0.35)),
        log(dat.sH1$prevax.sH1[dat.sH1$group=="TIVplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="TIVplacebo"]])

points( jitter(rep(2.65, length(dat.sH1$hhID[dat.sH1$group=="TIVplacebo"])),factor=7/(2.65/0.35)),
log(dat.sH1$postvax.sH1[dat.sH1$group=="TIVplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="TIVplacebo"]])

#00
points( jitter(rep(3.35, length(dat.sH1$hhID[dat.sH1$group=="placeboplacebo"])),factor=7/(3.35/0.35)),
        log(dat.sH1$prevax.sH1[dat.sH1$group=="placeboplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="placeboplacebo"]])

points( jitter(rep(3.65, length(dat.sH1$hhID[dat.sH1$group=="placeboplacebo"])),factor=7/(3.65/0.35)),
log(dat.sH1$postvax.sH1[dat.sH1$group=="placeboplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH1$age[dat.sH1$group=="placeboplacebo"]])

axis(1,pos=0,at=c(0,4),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c(0.5,1.5,2.5,3.5), label=c("TIV 2009-10\nTIV 2008-09", "TIV 2009-10\nPlacebo 2008-09", "Placebo 2009-10\nTIV 2008-09", "Placebo 2009-10\nPlacebo 2008-09"), cex.axis = 1)
    #axis(2, pos=0,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","","1:20","","1:80","","1:320","","1:1280","","1:5120","","1:20480"), pos=0, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Seasonal A(H1N1)", side = 3, line=-1, cex=1, font=2)
mtext("Antibody titers",side=2,line=2)

#sH3

plot(-2,-2, axes=FALSE, xlim=c(0, 4), ylim=c(0, 80), xlab="",ylab="")
text(0.35+0:3,1.5,"pre",cex=1)
text(0.65+0:3,1.5,"post",cex=1)
segments(0.5,76,1.5,76, lty=1)
segments(0.5,76,0.5,75, lty=1)
segments(1.5,76,1.5,75, lty=1)
rect(0.72, 71, 1.29,80,  col = 'white', border="white")
#text(0.8, 76, paste("(fold rise)\np=", round(wilcox.test(dat.sH3$postvax.sH3[dat.sH3$group=="placeboTIV"]/dat.sH3$prevax.sH3[dat.sH3$group=="placeboTIV"], dat.sH3$postvax.sH3[dat.sH3$group=="TIVTIV"]/dat.sH3$prevax.sH3[dat.sH3$group=="TIVTIV"])
#$p.value,2) , sep=""), cex=1)
text(1, 76, "(fold rise) p<0.01", cex=1, bg="white")

segments(0.65,70,1.65,70, lty=1)
segments(0.65,69,0.65,70, lty=1)
segments(1.65,70,1.65,69, lty=1)
rect(1.03, 67, 1.27,71,  col = 'white', border="white")
#text(1.15, 70, paste("p=", round(wilcox.test(dat.sH3$postvax.sH3[dat.sH3$group=="TIVTIV"], dat.sH3$postvax.sH3[dat.sH3$group=="placeboTIV"])
#$p.value,2) , sep=""), cex=1)
text(1.15, 70, "p<0.01", cex=1)

segments(0.35,6,1.35,6, lty=1)
segments(0.35,6,0.35,7, lty=1)
segments(1.35,6,1.35,7, lty=1)
rect(0.73, 3, 0.98,8,  col = 'white', border="white")
text(0.85, 6, paste("p=", round(wilcox.test(dat.sH3$prevax.sH3[dat.sH3$group=="TIVTIV"], dat.sH3$prevax.sH3[dat.sH3$group=="placeboTIV"])
$p.value,2) , sep=""), cex=1)
#text(0.85, 6, "p<0.01", cex=1)

segments(2.5,74,3.5,74, lty=1)
segments(3.5,74,3.5,73, lty=1)
segments(2.5,74,2.5,73, lty=1)
rect(2.72, 69, 3.29,78,  col = 'white', border="white")
text(3, 74, paste("(fold rise) p=", round(wilcox.test(dat.sH3$postvax.sH3[dat.sH3$group=="TIVplacebo"]/dat.sH3$prevax.sH3[dat.sH3$group=="TIVplacebo"], dat.sH3$postvax.sH3[dat.sH3$group=="placeboplacebo"]/dat.sH3$prevax.sH3[dat.sH3$group=="placeboplacebo"])$p.value,2) , sep=""), cex=1)

segments(0.5+0:3,0,0.5+0:3,65, lty=2)
segments(2.65,67,3.65,67, lty=1)
segments(2.65,67,2.65,66, lty=1)
segments(3.65,67,3.65,66, lty=1)
rect(3.03, 65, 3.27,69,  col = 'white', border="white")
#text(3.15, 67, paste("p=", round(wilcox.test(dat.sH3$postvax.sH3[dat.sH3$group=="TIVplacebo"], dat.sH3$postvax.sH3[dat.sH3$group=="placeboplacebo"])
#$p.value,2) , sep=""), cex=1, bg="white")
text(3.15, 67, "p<0.01", cex=1, bg="white")

segments(2.35,6,3.35,6, lty=1)
segments(2.35,6,2.35,7, lty=1)
segments(3.35,6,3.35,7, lty=1)
rect(2.78, 3, 3.03,8,  col = 'white', border="white")
text(2.9, 6, paste("p=", round(wilcox.test(dat.sH3$prevax.sH3[dat.sH3$group=="TIVplacebo"], dat.sH3$prevax.sH3[dat.sH3$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")

#11
points( jitter(rep(0.35, length(dat.sH3$hhID[dat.sH3$group=="TIVTIV"])),factor=7),
        log(dat.sH3$prevax.sH3[dat.sH3$group=="TIVTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="TIVTIV"]])

points( jitter(rep(0.65, length(dat.sH3$hhID[dat.sH3$group=="TIVTIV"])),factor=7/(0.65/0.35)),
log(dat.sH3$postvax.sH3[dat.sH3$group=="TIVTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="TIVTIV"]])

#10
points( jitter(rep(1.35, length(dat.sH3$hhID[dat.sH3$group=="placeboTIV"])),factor=7/(1.35/0.35)),
        log(dat.sH3$prevax.sH3[dat.sH3$group=="placeboTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="placeboTIV"]])

points( jitter(rep(1.65, length(dat.sH3$hhID[dat.sH3$group=="placeboTIV"])),factor=7/(1.65/0.35)),
log(dat.sH3$postvax.sH3[dat.sH3$group=="placeboTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="placeboTIV"]])

#01
points( jitter(rep(2.35, length(dat.sH3$hhID[dat.sH3$group=="TIVplacebo"])),factor=7/(2.35/0.35)),
        log(dat.sH3$prevax.sH3[dat.sH3$group=="TIVplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="TIVplacebo"]])

points( jitter(rep(2.65, length(dat.sH3$hhID[dat.sH3$group=="TIVplacebo"])),factor=7/(2.65/0.35)),
log(dat.sH3$postvax.sH3[dat.sH3$group=="TIVplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="TIVplacebo"]])

#00
points( jitter(rep(3.35, length(dat.sH3$hhID[dat.sH3$group=="placeboplacebo"])),factor=7/(3.35/0.35)),
        log(dat.sH3$prevax.sH3[dat.sH3$group=="placeboplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="placeboplacebo"]])

points( jitter(rep(3.65, length(dat.sH3$hhID[dat.sH3$group=="placeboplacebo"])),factor=7/(3.65/0.35)),
log(dat.sH3$postvax.sH3[dat.sH3$group=="placeboplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.sH3$age[dat.sH3$group=="placeboplacebo"]])

axis(1,pos=0,at=c(0,4),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c(0.5,1.5,2.5,3.5), label=c("TIV 2009-10\nTIV 2008-09", "TIV 2009-10\nPlacebo 2008-09", "Placebo 2009-10\nTIV 2008-09", "Placebo 2009-10\nPlacebo 2008-09"), cex.axis = 1)
    #axis(2, pos=0,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","","1:20","","1:80","","1:320","","1:1280","","1:5120","","1:20480"), pos=0, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Seasonal A(H3N2)", side = 3, line=-1, cex=1, font=2)
mtext("Antibody titers",side=2,line=2)


#B.Brisbane
plot(-2,-2, axes=FALSE, xlim=c(0, 4), ylim=c(0, 80), xlab="",ylab="")
text(0.35+0:3,1.5,"pre",cex=1)
text(0.65+0:3,1.5,"post",cex=1)
segments(0.5,76,1.5,76, lty=1)
segments(0.5,76,0.5,75, lty=1)
segments(1.5,76,1.5,75, lty=1)
rect(0.72, 71, 1.29,80,  col = 'white', border="white")
text(1, 76, paste("(fold rise) p=", round(wilcox.test(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="placeboTIV"]/dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="placeboTIV"], dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="TIVTIV"]/dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="TIVTIV"])
$p.value,2) , sep=""), cex=1)
#text(0.8, 76, "(fold rise)\np<0.01", cex=1, bg="white")

segments(0.65,70,1.65,70, lty=1)
segments(0.65,69,0.65,70, lty=1)
segments(1.65,70,1.65,69, lty=1)
rect(1.03, 67, 1.27,71,  col = 'white', border="white")
text(1.15, 70, paste("p=", round(wilcox.test(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="TIVTIV"], dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="placeboTIV"])
$p.value,2) , sep=""), cex=1)
#text(1.15, 70, "p<0.01", cex=1)

segments(0.35,6,1.35,6, lty=1)
segments(0.35,6,0.35,7, lty=1)
segments(1.35,6,1.35,7, lty=1)
rect(0.73, 3, 0.98,8,  col = 'white', border="white")
text(0.85, 6, paste("p=", round(wilcox.test(dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="TIVTIV"], dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="placeboTIV"])
$p.value,2) , sep=""), cex=1)
#text(0.85, 6, "p<0.01", cex=1)

segments(2.5,74,3.5,74, lty=1)
segments(3.5,74,3.5,73, lty=1)
segments(2.5,74,2.5,73, lty=1)
rect(2.72, 69, 3.29,78,  col = 'white', border="white")
text(3, 74, paste("(fold rise) p=", round(wilcox.test(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="TIVplacebo"]/dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="TIVplacebo"], dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="placeboplacebo"]/dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="placeboplacebo"])$p.value,2) , sep=""), cex=1)

segments(0.5+0:3,0,0.5+0:3,65, lty=2)
segments(2.65,67,3.65,67, lty=1)
segments(2.65,67,2.65,66, lty=1)
segments(3.65,67,3.65,66, lty=1)
rect(3.03, 65, 3.27,69,  col = 'white', border="white")
text(3.15, 67, paste("p=", round(wilcox.test(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="TIVplacebo"], dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")
#text(3.15, 67, "p<0.01", cex=1, bg="white")

segments(2.35,6,3.35,6, lty=1)
segments(2.35,6,2.35,7, lty=1)
segments(3.35,6,3.35,7, lty=1)
rect(2.78, 3, 3.03,8,  col = 'white', border="white")
text(2.9, 6, paste("p=", round(wilcox.test(dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="TIVplacebo"], dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")

#11
points( jitter(rep(0.35, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="TIVTIV"])),factor=7),
        log(dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="TIVTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="TIVTIV"]])

points( jitter(rep(0.65, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="TIVTIV"])),factor=7/(0.65/0.35)),
log(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="TIVTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="TIVTIV"]])

#10
points( jitter(rep(1.35, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="placeboTIV"])),factor=7/(1.35/0.35)),
        log(dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="placeboTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="placeboTIV"]])

points( jitter(rep(1.65, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="placeboTIV"])),factor=7/(1.65/0.35)),
log(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="placeboTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="placeboTIV"]])

#01
points( jitter(rep(2.35, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="TIVplacebo"])),factor=7/(2.35/0.35)),
        log(dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="TIVplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="TIVplacebo"]])

points( jitter(rep(2.65, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="TIVplacebo"])),factor=7/(2.65/0.35)),
log(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="TIVplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="TIVplacebo"]])

#00
points( jitter(rep(3.35, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="placeboplacebo"])),factor=7/(3.35/0.35)),
        log(dat.B.Brisbane$prevax.B.Brisbane[dat.B.Brisbane$group=="placeboplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="placeboplacebo"]])

points( jitter(rep(3.65, length(dat.B.Brisbane$hhID[dat.B.Brisbane$group=="placeboplacebo"])),factor=7/(3.65/0.35)),
log(dat.B.Brisbane$postvax.B.Brisbane[dat.B.Brisbane$group=="placeboplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.B.Brisbane$age[dat.B.Brisbane$group=="placeboplacebo"]])

axis(1,pos=0,at=c(0,4),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c(0.5,1.5,2.5,3.5), label=c("TIV 2009-10\nTIV 2008-09", "TIV 2009-10\nPlacebo 2008-09", "Placebo 2009-10\nTIV 2008-09", "Placebo 2009-10\nPlacebo 2008-09"), cex.axis = 1)
    #axis(2, pos=0,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","","1:20","","1:80","","1:320","","1:1280","","1:5120","","1:20480"), pos=0, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Seasonal B/Brisbane", side = 3, line=-1, cex=1, font=2)
mtext("Antibody titers",side=2,line=2)


#pH1
plot(-2,-2, axes=FALSE, xlim=c(0, 4), ylim=c(0, 80), xlab="",ylab="")
text(0.35+0:3,1.5,"pre",cex=1)
text(0.65+0:3,1.5,"post",cex=1)
segments(0.5,76,1.5,76, lty=1)
segments(0.5,76,0.5,75, lty=1)
segments(1.5,76,1.5,75, lty=1)
rect(0.72, 71, 1.29,80,  col = 'white', border="white")
text(1, 76, paste("(fold rise) p=", round(wilcox.test(dat.pH1$postvax.pH1[dat.pH1$group=="placeboTIV"]/dat.pH1$prevax.pH1[dat.pH1$group=="placeboTIV"], dat.pH1$postvax.pH1[dat.pH1$group=="TIVTIV"]/dat.pH1$prevax.pH1[dat.pH1$group=="TIVTIV"])
$p.value,2) , sep=""), cex=1)
#text(0.8, 76, "(fold rise)\np<0.01", cex=1, bg="white")

segments(0.65,70,1.65,70, lty=1)
segments(0.65,69,0.65,70, lty=1)
segments(1.65,70,1.65,69, lty=1)
rect(1.03, 67, 1.27,71,  col = 'white', border="white")
text(1.15, 70, paste("p=", round(wilcox.test(dat.pH1$postvax.pH1[dat.pH1$group=="TIVTIV"], dat.pH1$postvax.pH1[dat.pH1$group=="placeboTIV"])
$p.value,2) , sep=""), cex=1)
#text(1.15, 70, "p<0.01", cex=1)

segments(0.35,6,1.35,6, lty=1)
segments(0.35,6,0.35,7, lty=1)
segments(1.35,6,1.35,7, lty=1)
rect(0.73, 3, 0.98,8,  col = 'white', border="white")
text(0.85, 6, paste("p=", round(wilcox.test(dat.pH1$prevax.pH1[dat.pH1$group=="TIVTIV"], dat.pH1$prevax.pH1[dat.pH1$group=="placeboTIV"])
$p.value,2) , sep=""), cex=1)
#text(0.85, 6, "p<0.01", cex=1)

segments(2.5,74,3.5,74, lty=1)
segments(3.5,74,3.5,73, lty=1)
segments(2.5,74,2.5,73, lty=1)
rect(2.72, 69, 3.29,78,  col = 'white', border="white")
text(3, 74, paste("(fold rise) p=", round(wilcox.test(dat.pH1$postvax.pH1[dat.pH1$group=="TIVplacebo"]/dat.pH1$prevax.pH1[dat.pH1$group=="TIVplacebo"], dat.pH1$postvax.pH1[dat.pH1$group=="placeboplacebo"]/dat.pH1$prevax.pH1[dat.pH1$group=="placeboplacebo"])$p.value,2) , sep=""), cex=1)


segments(0.5+0:3,0,0.5+0:3,65, lty=2)
segments(2.65,67,3.65,67, lty=1)
segments(2.65,67,2.65,66, lty=1)
segments(3.65,67,3.65,66, lty=1)
rect(3.03, 65, 3.27,69,  col = 'white', border="white")
text(3.15, 67, paste("p=", round(wilcox.test(dat.pH1$postvax.pH1[dat.pH1$group=="TIVplacebo"], dat.pH1$postvax.pH1[dat.pH1$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")
#text(3.15, 67, "p<0.01", cex=1, bg="white")

segments(2.35,6,3.35,6, lty=1)
segments(2.35,6,2.35,7, lty=1)
segments(3.35,6,3.35,7, lty=1)
rect(2.78, 3, 3.03,8,  col = 'white', border="white")
text(2.9, 6, paste("p=", round(wilcox.test(dat.pH1$prevax.pH1[dat.pH1$group=="TIVplacebo"], dat.pH1$prevax.pH1[dat.pH1$group=="placeboplacebo"])
$p.value,2) , sep=""), cex=1, bg="white")

#11
points( jitter(rep(0.35, length(dat.pH1$hhID[dat.pH1$group=="TIVTIV"])),factor=7),
        log(dat.pH1$prevax.pH1[dat.pH1$group=="TIVTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="TIVTIV"]])

points( jitter(rep(0.65, length(dat.pH1$hhID[dat.pH1$group=="TIVTIV"])),factor=7/(0.65/0.35)),
log(dat.pH1$postvax.pH1[dat.pH1$group=="TIVTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="TIVTIV"]])

#10
points( jitter(rep(1.35, length(dat.pH1$hhID[dat.pH1$group=="placeboTIV"])),factor=7/(1.35/0.35)),
        log(dat.pH1$prevax.pH1[dat.pH1$group=="placeboTIV"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="placeboTIV"]])

points( jitter(rep(1.65, length(dat.pH1$hhID[dat.pH1$group=="placeboTIV"])),factor=7/(1.65/0.35)),
log(dat.pH1$postvax.pH1[dat.pH1$group=="placeboTIV"],1.15),
pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="placeboTIV"]])

#01
points( jitter(rep(2.35, length(dat.pH1$hhID[dat.pH1$group=="TIVplacebo"])),factor=7/(2.35/0.35)),
        log(dat.pH1$prevax.pH1[dat.pH1$group=="TIVplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="TIVplacebo"]])

points( jitter(rep(2.65, length(dat.pH1$hhID[dat.pH1$group=="TIVplacebo"])),factor=7/(2.65/0.35)),
log(dat.pH1$postvax.pH1[dat.pH1$group=="TIVplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="TIVplacebo"]])

#00
points( jitter(rep(3.35, length(dat.pH1$hhID[dat.pH1$group=="placeboplacebo"])),factor=7/(3.35/0.35)),
        log(dat.pH1$prevax.pH1[dat.pH1$group=="placeboplacebo"],1.15),
        pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="placeboplacebo"]])

points( jitter(rep(3.65, length(dat.pH1$hhID[dat.pH1$group=="placeboplacebo"])),factor=7/(3.65/0.35)),
log(dat.pH1$postvax.pH1[dat.pH1$group=="placeboplacebo"],1.15),
pch = 1, cex = 1.5, col=c1[dat.pH1$age[dat.pH1$group=="placeboplacebo"]])

axis(1,pos=0,at=c(0,4),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c(0.5,1.5,2.5,3.5), label=c("TIV 2009-10\nTIV 2008-09", "TIV 2009-10\nPlacebo 2008-09", "Placebo 2009-10\nTIV 2008-09", "Placebo 2009-10\nPlacebo 2008-09"), cex.axis = 1)
    axis(2, at= c(logy), labels=c("<1:10","","1:20","","1:80","","1:320","","1:1280","","1:5120","","1:20480"), pos=0, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Pandemic A(H1N1)", side = 3, line=-1, cex=1, font=2)
mtext("Antibody titers",side=2,line=2)


#
# End of script.
#


