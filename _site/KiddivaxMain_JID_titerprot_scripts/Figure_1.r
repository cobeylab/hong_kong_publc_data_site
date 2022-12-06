#
# R syntax to reproduce Figure 1 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

# Create a figure showing the timeline of vaccination, blood taking, completion and flu activitiy of the study
# with cumulataive % vaccinated and completed study increasing with calendar time
#

dir <- "http://sph.hku.hk/data/KiddivaxMainV2/"
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_1.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_2.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_3.r")

kdata <- read.csv(paste(dir, "serology.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))
fluproxy <- read.csv("http://www.hku.hk/bcowling/influenza/influenza_proxy_1998to2013.csv")

time0 <- "31/07/2009"
fluproxy <- fluproxy[614:679,]
fluproxy$ti <- as.numeric(as.Date(fluproxy$Week.ending,format="%d/%m/%Y")-as.Date(time0,format="%d/%m/%Y"))-1
swab$ti <-as.numeric(as.Date(as.character(swab$date), format="%Y%m%d")-as.Date(time0,format="%d/%m/%Y"))-1

d <-c(31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31)

pcr<-matrix(rep(NA, 4*length(d)), nrow=4, dimnames=list(c("pdm","Vic","Yam","H3"),c(0,cumsum(d)[1:length(d)-1])))

for ( i in 1:length(d)) { pcr[1,i] <-length(unique(swab$hhID[swab$member==0 & swab$ti<=as.numeric(as.character(colnames(pcr)[i+1])) & swab$Swine.H1=="P"]))
                          pcr[2,i] <-length(unique(swab$hhID[swab$member==0 & swab$ti<=as.numeric(as.character(colnames(pcr)[i+1])) & swab$FluB.subtype=="Victoria"]))
                          pcr[3,i] <-length(unique(swab$hhID[swab$member==0 & swab$ti<=as.numeric(as.character(colnames(pcr)[i+1])) & swab$FluB.subtype=="Yamagata"]))
                          pcr[4,i] <-length(unique(swab$hhID[swab$member==0 & swab$ti<=as.numeric(as.character(colnames(pcr)[i+1])) & swab$H3=="P"]))
                          }

pcr[1,2:length(d)] <-pcr[1,2:length(d)] - pcr[1,1:length(d)-1]
pcr[2,2:length(d)] <-pcr[2,2:length(d)] - pcr[2,1:length(d)-1]
pcr[3,2:length(d)] <-pcr[3,2:length(d)] - pcr[3,1:length(d)-1]
pcr[4,2:length(d)] <-pcr[4,2:length(d)] - pcr[4,1:length(d)-1]

pcr[,length(d)] <-0

proxy <- data.frame(ti=fluproxy$ti,proxy=fluproxy$B.proxy)
s.proxy <-data.frame(ti=1:max(proxy$ti), s.proxy=predict(smooth.spline(proxy$ti, proxy$proxy, spar=0.5),1:max(proxy$ti))$y)

proxy.p <- data.frame(ti=fluproxy$ti,proxy.p=(1/r)*fluproxy$A.H1N1pdm.proxy)
s.proxy.p <-data.frame(ti=1:max(proxy.p$ti), s.proxy.p=predict(smooth.spline(proxy.p$ti, proxy.p$proxy.p, spar=0.6),1:max(proxy.p$ti))$y)    #run rescale script first  to get r !!


dat <-kdata[kdata$member==0 & !is.na(kdata$postvax.B.Brisbane) & !is.na(kdata$post.season.B.Brisbane),]
dat$followup <-as.numeric(as.Date(dat$post.season, format="%d/%m/%Y")-as.Date(dat$postvax, format="%d/%m/%Y"))
dat$followup2 <-as.numeric(as.Date(dat$mid.season, format="%d/%m/%Y")-as.Date(dat$postvax, format="%d/%m/%Y"))
dat$calpre <-as.numeric(as.Date(dat$prevax, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calpostv <-as.numeric(as.Date(dat$postvax, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calposts<-as.numeric(as.Date(dat$post.season, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))


pre <-data.frame(table(dat$calpre))
pre$cumfreq <-cumsum (pre$Freq) / max(cumsum(pre$Freq))
pre$t<-as.numeric(as.character(pre$Var1))

postv <-data.frame(table(dat$calpostv))
postv$cumfreq <-cumsum (postv$Freq) / max(cumsum(postv$Freq))
postv$t<-as.numeric(as.character(postv$Var1))

posts <-data.frame(table(dat$calposts))
posts$cumfreq <-cumsum (posts$Freq) / max(cumsum(posts$Freq))
posts$t <-as.numeric(as.character(posts$Var1))

# Plot parameters  -----------------------------------------------------------------------------------------------------------------------
mycol <-c(rgb(0.1,0.1,0.1,alpha=0.1),rgb(0.1,0.1,0.1,alpha=0.6),rgb(0.1,0.1,0.1,alpha=0.03) )

x<-78:500

d <-c(31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31)
m <-c("A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D")

# Plot 

windows(width=12,height=11.5)
par(mar=c(5,5,0,2), oma=c(0,0,0,0), family="sans")
layout(matrix(1:6, nrow=3, byrow=T), widths=c(6,6,6,6,6,6), heights=c(3.5,3.5,3.5) )

# Timeline and isolate rate plot ~~~~~~~~~

plot(0,0,type="n",xlab="",ylab="", xlim=c(0,518+50),ylim=c(0,1.1), axes=F)

polygon(c(rep(postv$t, each=2), max(postv$t)+2, sort(c(rep(posts$t, each=2), max(posts$t)+2), decreasing=T), min(postv$t)),
        c(0,rep(postv$cumfreq, each=2),sort(c(0,rep(posts$cumfreq, each=2)), decreasing=T),0)
,border=F, col=mycol[3])

barplot(pcr/50, width=d, col=c(colors()[490], gray(0.5), colors()[584],colors()[404]), add=T, offset=0, axes=F, space=0, border="white", names.arg=rep("",length(d)))
legend(270,0.6,c("A(H3N2)","B(Yamagata)","B(Victoria)","A(H1N1)pdm09"), fill=c( colors()[404], colors()[584],gray(0.5),colors()[490]), border=NA, bty="n",y.intersp=0.7, x.intersp=0.2 )

lines( c(rep(pre$t, each=2), max(pre$t)+2), c(0,rep(pre$cumfreq, each=2)), col=mycol[2], lty=2)
lines( c(rep(postv$t, each=2), max(postv$t)+2), c(0,rep(postv$cumfreq, each=2)), col=mycol[2], lty=2)
lines( c(rep(posts$t, each=2), max(posts$t)+2), c(0,rep(posts$cumfreq, each=2)), col=mycol[2], lty=2)


d <-c(31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31)
m <-c("A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D")


axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=0, tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=0, tcl=-0.6)

axis(4,at=0:10*0.1, labels=0:10*10, pos=518, las=1, col=gray(0.5), col.axis=gray(0.5))
axis(2, at=0:5*(1/5), labels=(0:5)*50*(1/5), pos=0,las=1)
#text(cumsum(d)-14, -0.03, m, adj=0.5)

mtext("No. PCR-confirmed influenza", side=2, line=3)
mtext("Cumulative % of subjects", side=4, line=-1, col=mycol[2])
text(90, 0.9, "Received\nTIV/placebo", col=mycol[2], cex=1)
text(230, 0.9, "Post-\nvaccination\nserum", col=mycol[2], cex=1)
text(450, 0.9, "Completed\nfollow-up\n(serum)", col=mycol[2], cex=1)

segments(cumsum(d)[8], 3*1/5, cumsum(d)[9]+10,3*1/5, col=mycol[2], lty=2)
segments(cumsum(d)[8], 3*1/5+0.01, cumsum(d)[8],3*1/5-0.01, col=mycol[2])
segments(cumsum(d)[9]+10, 3*1/5+0.01, cumsum(d)[9]+10,3*1/5-0.01, col=mycol[2])
text(cumsum(d)[8]+20, 3*1/5+0.08, "mid-study\nserum", col=mycol[2], cex=1)

mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=0.5, font=1,cex=1, side=1)
mtext("A)", side=2, at=1.08, line=3, las=1)

# Plot proxy ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot(0,0,type="n",xlab="",ylab="", xlim=c(0,518+8),ylim=c(0,0.9), axes=F)

lines(s.proxy$s.proxy[80:500]/10~s.proxy$ti[80:500],col=gray(0.3))
lines(s.proxy.p$s.proxy.p[98:500]/10~s.proxy.p$ti[98:500], col=colors()[566])

d <-c(31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31)
m <-c("A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=0, tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=0, tcl=-0.6)
axis(2, at=0:4*(1/5), labels=(0:4)*(1/5)*10/1000, pos=0,las=1, cex.axis=1)

mtext("Proxy influenza activity", side=2, line=3)
text(238, 0.47, "B",col=gray(0.3))
text(380, 0.25, "A(H1N1)pdm09", col=colors()[566])
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=0.5, font=1,cex=1, side=1)
mtext("B)", side=2, at=0.90, line=3, las=1)

# plot estimated titer change ~~~~~~~~~~

plot(0,0,type="n",xlab="",ylab="", xlim=c(0,518+8),ylim=c(log2(2.5), log2(10240)), axes=F, main="")

points(datr.p1$calpostv, jitter(log2(datr.p1$postvax.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.p1$calmid, jitter(log2(datr.p1$mid.season.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.p1$calposts, jitter(log2(datr.p1$post.season.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.p0$calpostv, jitter(log2(datr.p0$postvax.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.p0$calmid, jitter(log2(datr.p0$mid.season.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.p0$calposts, jitter(log2(datr.p0$post.season.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)

lines(x, model.p0$coef[1]+x*model.p0$coef[2])
lines(x, model.p1$coef[1]+x*model.p1$coef[2] , col="red")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=log2(2.5), tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=log2(2.5), tcl=-0.6)
axis(2,at=c(log2(5*2^(-1:11))), labels=c("","<10",5*2^(1:11)), las=1, pos=0)
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=0.5, side=1, font=1)

mtext("C)", font=1, side=2, line=3, at=log2(10240)+0.8, las=1)
mtext("A(H1N1)pdm09", side=3, font=1, line=-1.2,  at=540/2)
mtext("HAI titer", side=2, line=3)

#Plot for B (Victoria)--------------------------------------------------------------------------------------------------------------------------------------

plot(0,0,type="n",xlab="",ylab="", xlim=c(0,518+8),ylim=c(log2(2.5), log2(10240)), axes=F, main="")

points(datr.b1$calpostv[datr.b1$age%in%9:18], jitter(log2(datr.b1$postvax.B.Brisbane[datr.b1$age%in%9:18]),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calmid[datr.b1$age%in%9:18], jitter(log2(datr.b1$mid.season.B.Brisbane[datr.b1$age%in%9:18]),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calposts[datr.b1$age%in%9:18], jitter(log2(datr.b1$post.season.B.Brisbane[datr.b1$age%in%9:18]),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.b1$calpostv[datr.b1$age%in%6:8], jitter(log2(datr.b1$postvax.B.Brisbane[datr.b1$age%in%6:8]),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calmid[datr.b1$age%in%6:8], jitter(log2(datr.b1$mid.season.B.Brisbane[datr.b1$age%in%6:8]),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calposts[datr.b1$age%in%6:8], jitter(log2(datr.b1$post.season.B.Brisbane[datr.b1$age%in%6:8]),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.b0$calpostv, jitter(log2(datr.b0$postvax.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.b0$calmid, jitter(log2(datr.b0$mid.season.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.b0$calposts, jitter(log2(datr.b0$post.season.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)

lines(x, model.b0$coef[1]+x*model.b0$coef[2])
lines(x, model.b1.age$coef[1]+x*model.b1.age$coef[2]+model.b1.age$coef[3]+model.b1.age$coef[4]*x , col="red", lty=2)
lines(x, model.b1.age$coef[1]+x*model.b1.age$coef[2] , col="red", lty=6)
text(350,log2(80),"9-17y", col="red")
text(350,log2(50),"6-8y", col="red")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=log2(2.5), tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=log2(2.5), tcl=-0.6)
axis(2,at=c(log2(5*2^(-1:11))), labels=c("","<10",5*2^(1:11)), las=1, pos=0)
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=0.5, side=1, font=1)
mtext("B(Victoria-lineage)", side=3, font=1, line=-1.2, at=540/2)
mtext("D)", font=1, side=2, line=3, at=log2(10240)+0.8, las=1)

# Plot HIP curve ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta <-cox.td.A.m$coef
se <-sqrt(sum(vcov(cox.td.A.m)))
est <-data.frame(titer=log2(5*2^(0:9)),proxy=2, mu=NA, ci1=NA, ci2=NA)
est$mu <-exp(beta[1] * est$titer + beta[2]*est$proxy  )
est$ci1 <-exp(beta[1] * est$titer + beta[2]*est$proxy  - 1.96*se)
est$ci2 <-exp(beta[1] * est$titer + beta[2]*est$proxy + 1.96*se)
est$d <-(est$mu[1]-est$mu)/est$mu[1]
est$d.ci1 <-(est$mu[1]-est$ci2)/est$mu[1]
est$d.ci2 <-(est$mu[1]-est$ci1)/est$mu[1]

plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(log2(5),log2(2560)), ylim=c(0,1), type="n")
axis(1,at=log2(5*2^(0:9)), labels=5*2^(0:9))
axis(2,at=0.2*0:5, labels=0.2*0:5, las=1, pos=log2(5))

lines(est$titer,est$d)
lines(est$titer,est$d.ci1, lty=2)
lines(est$titer,est$d.ci2, lty=2)

mtext("HAI titer", side=1, line=4, at=(log2(80)+log2(160))/2)
mtext("Relative risk reduction", side=2, line=3)
mtext("A(H1N1)pdm09",side=3,line=-1.2, font=1, at=(log2(5)+log2(2560))/2)
mtext("E)", font=1, side=2, line=3, at=1.08, las=1)

#plot for flu B  (TIV group)
beta <-vic.vic$coef
se <-sqrt(sum(vcov(cox.td.B)))

est <-data.frame(titer=log2(5*2^(0:9)),proxy=2, mu=NA, ci1=NA, ci2=NA)
est$mu <-exp(beta[1] * est$titer + beta[2]*est$proxy  )
est$ci1 <-exp(beta[1] * est$titer + beta[2]*est$proxy  - 1.96*se)
est$ci2 <-exp(beta[1] * est$titer + beta[2]*est$proxy + 1.96*se)
est$d <-(est$mu[1]-est$mu)/est$mu[1]
est$d.ci1 <-(est$mu[1]-est$ci2)/est$mu[1]
est$d.ci2 <-(est$mu[1]-est$ci1)/est$mu[1]

plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(log2(5),log2(2560)), ylim=c(0,1),type="n")
axis(1,at=log2(5*2^(0:9)), labels=5*2^(0:9))
axis(2,at=0.2*0:5, labels=0.2*0:5, las=1, pos=log2(5))

lines(est$titer,est$d)
lines(est$titer,est$d.ci1, lty=2)
lines(est$titer,est$d.ci2, lty=2)

mtext("HAI titer", side=1, line=4, at=(log2(80)+log2(160))/2)
mtext("B(Victoria-lineage)",side=3,line=-1.2, font=1, at=(log2(5)+log2(2560))/2)
mtext("F)", font=1, side=2, line=3, at=1.08, las=1)


#
# End of script.
#


