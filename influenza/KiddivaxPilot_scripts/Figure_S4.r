#
# R syntax to reproduce information for Figure S4 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 13, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe.r")

#values for y axis
y <-c(5,10,20,40,80,160,320,640,1280)
logy <-log(y,1.15)

###

windows(width=8, height=4)
par(mar=c(1,1,0,1), oma=c(0,0.5,0,0))
c1 <-"black"

plot(-2,-2, axes=FALSE, xlim=c(0, 4.5), ylim=c(0, 55), xlab="",ylab="")
points( jitter(rep(1, length(all$hhID[(all$sh1w==0 & all$sh3w==0 & all$flo.w==0)])),factor=4),
        log(all$ph1.mid.titer[(all$sh1w==0 & all$sh3w==0 & all$flo.w==0)],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(2, length(all$hhID[(all$sh1w==1 |all$swab.sh1w)])),factor=2),
        log(all$ph1.mid.titer[(all$sh1w==1 |all$swab.sh1w )],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(3, length(all$hhID[(all$sh3w==1 |all$swab.sh3w )])),factor=1.5),
        log(all$ph1.mid.titer[(all$sh3w==1 |all$swab.sh3w)],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(4, length(all$hhID[(all$flo.w==1)])),factor=1),
        log(all$ph1.mid.titer[(all$flo.w==1)],1.15),
        pch = 1, cex = 1, col=c1)

axis(1,pos=5,at=c(0.5,4.5),col.ticks="white",labels=F)
axis(1, pos=5, lwd=1, at=c(1, 2,3,4), label=c("No seasonal\ninfection", "Lab-confirmed\nA/H1N1", "Lab-confirmed\nA/H3N2", "Lab-confirmed\nB"),padj=0.3)
    axis(2, pos=0.5,at=c(5,55),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280"), pos=0.5, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Antibody titers to pandemic A/H1N1", side = 2, line =-1,adj=0.7)

lines(c(1.3,1.3),c(quantile(log(all$ph1.mid.titer[(all$sh1w==0 & all$sh3w==0)],1.15),0.25, na.rm=T),quantile(log(all$ph1.mid.titer[(all$sh1w==0 & all$sh3w==0)],1.15),0.75, na.rm=T)),lty=1,col="gray")
points(1.3,median(log(all$ph1.mid.titer[(all$sh1w==0 & all$sh3w==0)],1.15),na.rm=T),pch="-",bg="black",col="black")

lines(c(2.3,2.3),c(quantile(log(all$ph1.mid.titer[(all$sh1w==1 |all$swab.sh1w==1)],1.15),0.25, na.rm=T),quantile(log(all$ph1.mid.titer[(all$sh1w==1 |all$swab.sh1w==1 )],1.15),0.75, na.rm=T)),lty=1,col="black")
points(2.3,median(log(all$ph1.mid.titer[(all$sh1w==1 |all$swab.sh1w==1 )],1.15), na.rm=T),pch="-",bg="black",col="black")

lines(c(3.3,3.3),c(quantile(log(all$ph1.mid.titer[( all$sh3w==1 |all$swab.sh3w==1)],1.15),0.25, na.rm=T),quantile(log(all$ph1.mid.titer[(all$sh3w==1 |all$swab.sh3w==1)],1.15),0.75, na.rm=T)),lty=1,col="black")
points(3.3,median(log(all$ph1.mid.titer[(all$sh3w==1 |all$swab.sh3w==1)],1.15), na.rm=T),pch="-",bg="black",col="black")

lines(c(4.3,4.3),c(quantile(log(all$ph1.mid.titer[(all$flo.w==1)],1.15),0.25, na.rm=T),quantile(log(all$ph1.mid.titer[(all$flo.w==1)],1.15),0.75, na.rm=T)),lty=1,col="gray")
points(4.3,median(log(all$ph1.mid.titer[(all$flo.w==1)],1.15),na.rm=T),pch="-",bg="black",col="black")

lines(c(1.3,2.3),c(45,45))
lines(c(1.3,1.3),c(45,44))
lines(c(2.3,2.3),c(45,44))

points(1.8,45,pch=15,col="white",cex=7)
text(1.8,45,
paste("p=",round(
wilcox.test(
all$ph1.mid.titer[(all$sh1w==1 | all$swab.sh1w==1)],
all$ph1.mid.titer[(all$sh1w==0 & all$sh3w==0)])$p.value,2),sep=""),cex=0.9)

lines(c(1.3,3.3),c(50,50))
lines(c(1.3,1.3),c(50,49))
lines(c(3.3,3.3),c(50,49))

points(2.3,50,pch=15,col="white",cex=7)
text(2.3,50,
paste("p=",round(
wilcox.test(
all$ph1.mid.titer[(all$sh1w==3 | all$swab.sh3w==1)],
all$ph1.mid.titer[(all$sh1w==0 & all$sh3w==0)])$p.value,2),sep=""),cex=0.9)

lines(c(1.3,4.3),c(55,55))
lines(c(1.3,1.3),c(55,54))
lines(c(4.3,4.3),c(55,54))

points(2.8,55,pch=15,col="white",cex=7)
text(2.8,55,"p=1.00",cex=0.9)

#
# End of script.
#


