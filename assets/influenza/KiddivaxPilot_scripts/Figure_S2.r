#
# R syntax to reproduce information for Figure S2 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 13, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe.r")

c1<-"gray"
c2<-"black"

#x marker for intervention group 1==TIV, 0==placebo
index.pre$gp[index.pre$intervention=="TIV"] <-0
index.pre$gp[index.pre$intervention=="placebo"] <-1

#values for y axis
y <-c(5,10,20,40,80,160,320,640,1280,2560,5120,10240,20480)
logy <-log(y,1.15)

############################################################################
# Plot col 1: postvac Titers Vs winter infection by intervention gp        #
############################################################################


## Pre-post vac seasonal Titers -------------------------------------------------------------------------
windows(height=6, width=7)
layout(matrix(1:4, ncol=2, byrow=TRUE))
par(mar=c(2,4,1.5,0), oma=c(0,2,0,0))

#sh1

plot(-2,-2, axes=FALSE, xlim=c(-0.5, 3.5), ylim=c(0, 75), xlab="",ylab="")

lines(c(0.5,2.5),c(75,75))
lines(c(0.5,0.5),c(75,74))
lines(c(2.5,2.5),c(75,74))
points(1.5,75,pch=22,col="white",bg="white",cex=9)
text(1.5,75,"p<0.01")

points( jitter(rep(2, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$sh1.pre[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(2.5, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$sh1.postv[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c2)

points( jitter(rep(0, length(index.pre$gp[index.pre$gp==0])),factor=6),
        log(index.pre$sh1.pre[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(0.5, length(index.pre$gp[index.pre$gp==0])),factor=10),
        log(index.pre$sh1.postv[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c2)

axis(1,pos=0,at=c(-0.5,3.5),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c( 0.3, 2.3), label=c("TIV", "Placebo"), cex.axis = 1)
    axis(2, pos=-0.5,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560","1:5120","1:10240","1:20480"),
         pos=-0.5, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))


mtext("Seasonal A/H1N1", side = 3, line=0, cex=1, font=2)
mtext("Antibody titers",side=2,line=3.5)

text(0.25,5,"p<0.01",cex=1)

text(2.25,5,
paste("p=",round(
wilcox.test(
index.pre$sh1.pre[index.pre$gp==1],
index.pre$sh1.postv[index.pre$gp==1])$p.value,2)),cex=1)


lines(c(2.8,2.8),c(quantile(log(index.pre$sh1.pre[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$sh1.pre[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(2.8,median(log(index.pre$sh1.pre[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="gray",col="gray",cex=1)

lines(c(3.0,3.0),c(quantile(log(index.pre$sh1.postv[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$sh1.postv[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c2)
points(3.0,median(log(index.pre$sh1.postv[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="black",cex=1)

lines(c(0.8,0.8),c(quantile(log(index.pre$sh1.pre[index.pre$gp==0],1.15),0.25,na.rm=T),quantile(log(index.pre$sh1.pre[index.pre$gp==0],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(0.8,median(log(index.pre$sh1.pre[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="gray",col="gray",cex=1)

lines(c(1.0,1.0),c(quantile(log(index.pre$sh1.postv[index.pre$gp==0],1.15),0.25, na.rm=T),quantile(log(index.pre$sh1.postv[index.pre$gp==0],1.15),0.75, na.rm=T)),lty=1,col=c2)
points(1.0,median(log(index.pre$sh1.postv[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="black",cex=1)

lines(c(0,0.5),c(8,8))
lines(c(2,2.5),c(8,8))
lines(c(0,0),c(8,9))
lines(c(2,2),c(8,9))
lines(c(0.5,0.5),c(8,9))
lines(c(2.5,2.5),c(8,9))



#sH3  ----------------------------------------------------------------------------

plot(-2,-2, axes=FALSE, xlim=c(-0.5, 3.5), ylim=c(0, 75), xlab="",ylab="")

lines(c(0.5,2.5),c(75,75))
lines(c(0.5,0.5),c(75,74))
lines(c(2.5,2.5),c(75,74))
points(1.5,75,pch=22,col="white",bg="white",cex=9)
text(1.5,75,"p<0.01")

points( jitter(rep(2, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$sh3.pre[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(2.5, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$sh3.postv[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c2)

points( jitter(rep(0, length(index.pre$gp[index.pre$gp==0])),factor=6),
        log(index.pre$sh3.pre[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(0.5, length(index.pre$gp[index.pre$gp==0])),factor=10),
        log(index.pre$sh3.postv[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c2)

axis(1,pos=0,at=c(-0.5,3.5),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c( 0.3, 2.3), label=c("TIV", "Placebo"), cex.axis = 1)
    axis(2, pos=-0.5,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560","1:5120","1:10240","1:20480"), pos=-0.5, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Seasonal A/H3N2", side = 3, line = 0, cex=1, font=2)

text(0.25,5,"p<0.01",cex=1)

text(2.25,5,
paste("p=",round(
wilcox.test(
index.pre$sh3.pre[index.pre$gp==1],
index.pre$sh3.postv[index.pre$gp==1])$p.value,2)),cex=1)


lines(c(2.8,2.8),c(quantile(log(index.pre$sh3.pre[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$sh3.pre[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(2.8,median(log(index.pre$sh3.pre[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="gray",col="gray")

lines(c(3.0,3.0),c(quantile(log(index.pre$sh3.postv[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$sh3.postv[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c2)
points(3.0,median(log(index.pre$sh3.postv[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="black")

lines(c(0.8,0.8),c(quantile(log(index.pre$sh3.pre[index.pre$gp==0],1.15),0.25,na.rm=T),quantile(log(index.pre$sh3.pre[index.pre$gp==0],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(0.8,median(log(index.pre$sh3.pre[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="gray",col="gray")

lines(c(1.0,1.0),c(quantile(log(index.pre$sh3.postv[index.pre$gp==0],1.15),0.25,na.rm=T),quantile(log(index.pre$sh3.postv[index.pre$gp==0],1.15),0.75,na.rm=T)),lty=1,col=c2)
points(1.0,median(log(index.pre$sh3.postv[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="black")


lines(c(0,0.5),c(8,8))
lines(c(2,2.5),c(8,8))
lines(c(0,0),c(8,9))
lines(c(2,2),c(8,9))
lines(c(0.5,0.5),c(8,9))
lines(c(2.5,2.5),c(8,9))

#b

plot(-2,-2, axes=FALSE, xlim=c(-0.5, 3.5), ylim=c(0, 75), xlab="",ylab="")

lines(c(0.5,2.5),c(75,75))
lines(c(0.5,0.5),c(75,74))
lines(c(2.5,2.5),c(75,74))
points(1.5,75,pch=22,col="white",bg="white",cex=9)
text(1.5,75,"p<0.01")

points( jitter(rep(2, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$FluB.Florida.pre[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(2.5, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$FluB.Florida.postv[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c2)

points( jitter(rep(0, length(index.pre$gp[index.pre$gp==0])),factor=6),
        log(index.pre$FluB.Florida.pre[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(0.5, length(index.pre$gp[index.pre$gp==0])),factor=10),
        log(index.pre$FluB.Florida.postv[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c2)

axis(1,pos=0,at=c(-0.5,3.5),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c( 0.3, 2.3), label=c("TIV", "Placebo"), cex.axis = 1)
    axis(2, pos=-0.5,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560","1:5120","1:10240","1:20480"), pos=-0.5, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Seasonal B", side = 3, line = 0, cex=1, font=2)
mtext("Antibody titers",side=2,line=3.5)

text(0.25,5,"p<0.01",cex=1)

text(2.25,5,
paste("p=",round(
wilcox.test(
index.pre$FluB.Florida.pre[index.pre$gp==1],
index.pre$FluB.Florida.postv[index.pre$gp==1])$p.value,2)),cex=1)


lines(c(2.8,2.8),c(quantile(log(index.pre$FluB.Florida.pre[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$FluB.Florida.pre[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(2.8,median(log(index.pre$FluB.Florida.pre[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="gray",col="gray",cex=1)

lines(c(3.0,3.0),c(quantile(log(index.pre$FluB.Florida.postv[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$FluB.Florida.postv[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c2)
points(3.0,median(log(index.pre$FluB.Florida.postv[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="black",cex=1)

lines(c(0.8,0.8),c(quantile(log(index.pre$FluB.Florida.pre[index.pre$gp==0],1.15),0.25,na.rm=T),quantile(log(index.pre$FluB.Florida.pre[index.pre$gp==0],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(0.8,median(log(index.pre$FluB.Florida.pre[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="gray",col="gray",cex=1)

lines(c(1.0,1.0),c(quantile(log(index.pre$FluB.Florida.postv[index.pre$gp==0],1.15),0.25, na.rm=T),quantile(log(index.pre$FluB.Florida.postv[index.pre$gp==0],1.15),0.75, na.rm=T)),lty=1,col=c2)
points(1.0,median(log(index.pre$FluB.Florida.postv[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="black",cex=1)

lines(c(0,0.5),c(8,8))
lines(c(2,2.5),c(8,8))
lines(c(0,0),c(8,9))
lines(c(2,2),c(8,9))
lines(c(0.5,0.5),c(8,9))
lines(c(2.5,2.5),c(8,9))


#ph1  ----------------------------------------------------------------------------

plot(-2,-2, axes=FALSE, xlim=c(-0.5, 3.5), ylim=c(0, 75), xlab="",ylab="")
points( jitter(rep(2, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$ph1.pre[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(2.5, length(index.pre$gp[index.pre$gp==1])),factor=2),
        log(index.pre$ph1.postv[index.pre$gp==1],1.15),
        pch = 1, cex = 1, col=c2)

points( jitter(rep(0, length(index.pre$gp[index.pre$gp==0])),factor=6),
        log(index.pre$ph1.pre[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c1)

points( jitter(rep(0.5, length(index.pre$gp[index.pre$gp==0])),factor=10),
        log(index.pre$ph1.postv[index.pre$gp==0],1.15),
        pch = 1, cex = 1, col=c2)

axis(1,pos=0,at=c(-0.5,3.5),col.ticks="white",labels=F)
axis(1, pos=0, lwd=1, font=1, at=c( 0.3, 2.3), label=c("TIV", "Placebo"), cex.axis = 1)
    axis(2, pos=-0.5,at=c(0,75),col.ticks="white",labels=F)
    axis(2, at= c(logy), labels=c("<1:10","1:10","1:20","1:40","1:80","1:160","1:320","1:640","1:1280","1:2560","1:5120","1:10240","1:20480"), pos=-0.5, lwd=1, cex.axis=1,las=1, mgp=c(1.5,0.7,0))

mtext("Pandemic A/H1N1", side = 3, line = 0, cex=1, font=2)

text(0.25,5,"p<0.01",cex=1)

text(2.25,5,
paste("p=",round(
wilcox.test(
index.pre$ph1.pre[index.pre$gp==1],
index.pre$ph1.postv[index.pre$gp==1])$p.value,2)),cex=1)


lines(c(2.8,2.8),c(quantile(log(index.pre$ph1.pre[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$ph1.pre[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(2.8,median(log(index.pre$ph1.pre[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="gray",col="gray")

lines(c(3.0,3.0),c(quantile(log(index.pre$ph1.postv[index.pre$gp==1],1.15),0.25,na.rm=T),quantile(log(index.pre$ph1.postv[index.pre$gp==1],1.15),0.75,na.rm=T)),lty=1,col=c2)
points(3.0,median(log(index.pre$ph1.postv[index.pre$gp==1],1.15),na.rm=T),pch="-",bg="black")

lines(c(0.8,0.8),c(quantile(log(index.pre$ph1.pre[index.pre$gp==0],1.15),0.25,na.rm=T),quantile(log(index.pre$ph1.pre[index.pre$gp==0],1.15),0.75,na.rm=T)),lty=1,col=c1)
points(0.8,median(log(index.pre$ph1.pre[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="gray",col="gray")

lines(c(1.0,1.0),c(quantile(log(index.pre$ph1.postv[index.pre$gp==0],1.15),0.25,na.rm=T),quantile(log(index.pre$ph1.postv[index.pre$gp==0],1.15),0.75,na.rm=T)),lty=1,col=c2)
points(1.0,median(log(index.pre$ph1.postv[index.pre$gp==0],1.15),na.rm=T),pch="-",bg="black")

lines(c(0.5,2.5),c(30,30))
lines(c(0.5,0.5),c(30,29))
lines(c(2.5,2.5),c(30,29))
points(1.5,31,pch=22,col="white",bg="white",cex=9)
text(1.5,31,"p<0.01")

lines(c(0,0.5),c(8,8))
lines(c(2,2.5),c(8,8))
lines(c(0,0),c(8,9))
lines(c(2,2),c(8,9))
lines(c(0.5,0.5),c(8,9))
lines(c(2.5,2.5),c(8,9))


#
# End of script.
#





