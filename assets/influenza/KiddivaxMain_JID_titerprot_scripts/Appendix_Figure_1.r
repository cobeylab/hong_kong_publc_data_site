#
# R syntax to reproduce Appendix Figure 1 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

#
# Statistical model to evaluate the correlation between HI titer and subsequent protection
# Cox model taking into account : Background risk of infection (Lab * ILI), waning HI antibodies + beta*titer(t0)
# To smooth the proxy line, spline(x,y,xout=1:max(t)).

#
#
# Read in main results
dir <- "http://sph.hku.hk/data/KiddivaxMainV2/"
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_1.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_2.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_3.r")

#
# Pandemic A(H1N1)  Model
#

# Read in data (Received from Vicky on 2012-12-12)
kdata <- read.csv(paste(dir, "serology.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))
fluproxy <- read.csv("http://www.hku.hk/bcowling/influenza/influenza_proxy_1998to2013.csv")

time0 <- "01/10/2009"
fluproxy <- fluproxy[614:679,]
fluproxy$ti <- as.numeric(as.Date(fluproxy$Week.ending,format="%d/%m/%Y")-as.Date(time0,format="%d/%m/%Y"))-1
proxy <- data.frame(ti=ili$ti,proxy=fluproxy$A.H1N1pdm.proxy)
# create smoothed proxy for background risk
s.proxy <-data.frame(ti=1:max(proxy$ti), s.proxy=predict(smooth.spline(proxy$ti, proxy$proxy, spar=0.6),1:max(proxy$ti))$y)

kdata <- kdata[kdata$member==0,]
kdata <-merge(kdata,random, by="hhID",all.x=T)

swab <- swab[swab$Swine.H1=="P",]; swab <- swab[order(swab$hhID,swab$member,swab$date),]
swab$mark <- 0
for(i in 2:nrow(swab)){
   if(swab$hhID[i]==swab$hhID[i-1]&swab$member[i]==swab$member[i-1]) swab$mark[i] <- 1
}
swab <- swab[swab$mark==0&swab$member==0,]

# merge in PCR type and date
kdata$swab.pH1 <- 1*(kdata$hhID%in%swab$hhID[swab$Swine.H1=="P"])
kdata <- kdata[!is.na(kdata$postvax.pH1)&kdata$postvax!="",]
kdata2 <- merge(kdata,swab[swab$Swine.H1=="P",c(1,3)],all.x=T)
kdata2$t0 <- as.numeric(as.Date(as.character(kdata2$postvax),format="%d/%m/%Y")-as.Date(as.character("1/10/2009"),format="%d/%m/%Y"))+1
kdata2$T <- as.numeric(as.Date(as.character(kdata2$date),format="%Y%m%d")-as.Date(as.character("1/10/2009"),format="%d/%m/%Y"))+1
kdata2$T[is.na(kdata2$T)] <- as.numeric(as.Date(as.character(kdata2$end.date[is.na(kdata2$T)]),format="%d/%m/%Y")-as.Date(as.character("1/10/2009"),format="%d/%m/%Y"))+1

kdata3 <- kdata2[kdata2$t0<kdata2$T,]
kdata3.p1 <-kdata3[kdata3$intervention=="TIV",]
kdata3.p0 <-kdata3[kdata3$intervention=="placebo",]

mrow <- kdata3.p1$T[1]-kdata3.p1$t0[1]+1
kdata.td.p1 <- data.frame(hhID=rep(kdata3.p1$hhID[1],mrow),TIV=1*(kdata3.p1$intervention[1]=="TIV"),postvax.P=rep(log2(kdata3.p1$postvax.pH1[1]),mrow),t=(1:mrow)-1,
                       event=c(rep(0,mrow-1),kdata3.p1$swab.pH1[1]))
kdata.td.p1$proxy <- s.proxy$s.proxy[kdata3.p1$t0[1]:kdata3.p1$T[1]]
kdata.td2.p1 <- kdata.td.p1[nrow(kdata.td.p1),]

for(j in 2:nrow(kdata3.p1)){
   mrow <- kdata3.p1$T[j]-kdata3.p1$t0[j]+1
   temp <- data.frame(hhID=rep(kdata3.p1$hhID[j],mrow),TIV=1*(kdata3.p1$intervention[j]=="TIV"),postvax.P=rep(log2(kdata3.p1$postvax.pH1[j]),mrow),t=(1:mrow)-1,
                      event=c(rep(0,mrow-1),kdata3.p1$swab.pH1[j]))
   temp$proxy <- s.proxy$s.proxy[kdata3.p1$t0[j]:kdata3.p1$T[j]]
   temp2 <- temp[nrow(temp),]
   kdata.td.p1 <- rbind(kdata.td.p1,temp)
   kdata.td2.p1 <- rbind(kdata.td2.p1,temp2)
}

mrow <- kdata3.p0$T[1]-kdata3.p0$t0[1]+1
kdata.td.p0 <- data.frame(hhID=rep(kdata3.p0$hhID[1],mrow),TIV=1*(kdata3.p0$intervention[1]=="TIV"),postvax.P=rep(log2(kdata3.p0$postvax.pH1[1]),mrow),t=(1:mrow)-1,
                       event=c(rep(0,mrow-1),kdata3.p0$swab.pH1[1]))
kdata.td.p0$proxy <- s.proxy$s.proxy[kdata3.p0$t0[1]:kdata3.p0$T[1]]
kdata.td2.p0 <- kdata.td.p0[nrow(kdata.td.p0),]

for(j in 2:nrow(kdata3.p0)){
   mrow <- kdata3.p0$T[j]-kdata3.p0$t0[j]+1
   temp <- data.frame(hhID=rep(kdata3.p0$hhID[j],mrow),TIV=1*(kdata3.p0$intervention[j]=="TIV"),postvax.P=rep(log2(kdata3.p0$postvax.pH1[j]),mrow),t=(1:mrow)-1,
                      event=c(rep(0,mrow-1),kdata3.p0$swab.pH1[j]))
   temp$proxy <- s.proxy$s.proxy[kdata3.p0$t0[j]:kdata3.p0$T[j]]
   temp2 <- temp[nrow(temp),]
   kdata.td.p0 <- rbind(kdata.td.p0,temp)
   kdata.td2.p0 <- rbind(kdata.td2.p0,temp2)
}

kdata.td.p <-rbind(kdata.td.p0, kdata.td.p1)

cox.td.A <- coxph(Surv(t,event)~postvax.P+proxy,data=kdata.td.p)

#---------------------------------------------------------------------------------------------------------------------------------

#
# Influenza B  Model
#
# Read in data (Received from Vicky on 2012-12-12)
kdata <- read.csv(paste(dir, "serology.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))
fluproxy <- read.csv("http://www.hku.hk/bcowling/influenza/influenza_proxy_1998to2013.csv")

time0 <- "01/10/2009"
fluproxy <- fluproxy[614:679,]
fluproxy$ti <- as.numeric(as.Date(fluproxy$Week.ending,format="%d/%m/%Y")-as.Date(time0,format="%d/%m/%Y"))-1
proxy <- data.frame(ti=ili$ti,proxy=(fluproxy$B.proxy))
# create smoothed proxy for background risk
s.proxy <-data.frame(ti=1:max(proxy$ti), s.proxy=predict(smooth.spline(proxy$ti, proxy$proxy, spar=0.5),1:max(proxy$ti))$y)

kdata <- kdata[kdata$member==0,]
kdata <-merge(kdata,random, by="hhID",all.x=T)

swab <- swab[swab$FluB=="P",]; swab <- swab[order(swab$hhID,swab$member,swab$date),]
swab$mark <- 0
for(i in 2:nrow(swab)){
   if(swab$hhID[i]==swab$hhID[i-1]&swab$member[i]==swab$member[i-1]) swab$mark[i] <- 1
}
swab <- swab[swab$mark==0&swab$member==0,]

# merge in PCR type and date
kdata$swab.Victoria <- 1*(kdata$hhID%in%swab$hhID[swab$FluB.subtype=="Victoria"])
kdata$swab.Yamagata <- 1*(kdata$hhID%in%swab$hhID[swab$FluB.subtype=="Yamagata"])
kdata <- kdata[!is.na(kdata$postvax.B.Brisbane)&kdata$postvax!="",]
kdata2 <- merge(kdata,swab[swab$FluB.subtype=="Victoria"|swab$FluB.subtype=="Yamagata",c(1,3)],all.x=T)
kdata2$t0 <- as.numeric(as.Date(as.character(kdata2$postvax),format="%d/%m/%Y")-as.Date(as.character("1/10/2009"),format="%d/%m/%Y"))+1
kdata2$T <- as.numeric(as.Date(as.character(kdata2$date),format="%Y%m%d")-as.Date(as.character("1/10/2009"),format="%d/%m/%Y"))+1
kdata2$T[is.na(kdata2$T)] <- as.numeric(as.Date(as.character(kdata2$end.date[is.na(kdata2$T)]),format="%d/%m/%Y")-as.Date(as.character("1/10/2009"),format="%d/%m/%Y"))+1

kdata3 <- kdata2[kdata2$t0<kdata2$T,]
kdata3.b1 <-kdata3[kdata3$intervention=="TIV",]
kdata3.b0 <-kdata3[kdata3$intervention=="placebo",]

mrow <- kdata3.b1$T[1]-kdata3.b1$t0[1]+1
kdata.td.b1 <- data.frame(hhID=rep(kdata3.b1$hhID[1],mrow),TIV=1*(kdata3.b1$intervention[1]=="TIV"),postvax.b=rep(log2(kdata3.b1$postvax.B.Brisbane[1]),mrow),t=(1:mrow)-1,
                       event=c(rep(0,mrow-1),1*(kdata3.b1$swab.Victoria[1]==1 |kdata3.b1$swab.Yamagata[1]==1 )))
kdata.td.b1$proxy <- s.proxy$s.proxy[kdata3.b1$t0[1]:kdata3.b1$T[1]]
kdata.td2.b1 <- kdata.td.b1[nrow(kdata.td.b1),]

for(j in 2:nrow(kdata3.b1)){
   mrow <- kdata3.b1$T[j]-kdata3.b1$t0[j]+1
   temp <- data.frame(hhID=rep(kdata3.b1$hhID[j],mrow),TIV=1*(kdata3.b1$intervention[j]=="TIV"),postvax.b=rep(log2(kdata3.b1$postvax.B.Brisbane[j]),mrow),t=(1:mrow)-1,
                      event=c(rep(0,mrow-1),1*(kdata3.b1$swab.Victoria[j]==1 |kdata3.b1$swab.Yamagata[j]==1 )) )
   temp$proxy <- s.proxy$s.proxy[kdata3.b1$t0[j]:kdata3.b1$T[j]]
   temp2 <- temp[nrow(temp),]
   kdata.td.b1 <- rbind(kdata.td.b1,temp)
   kdata.td2.b1 <- rbind(kdata.td2.b1,temp2)
}

mrow <- kdata3.b0$T[1]-kdata3.b0$t0[1]+1
kdata.td.b0 <- data.frame(hhID=rep(kdata3.b0$hhID[1],mrow),TIV=1*(kdata3.b0$intervention[1]=="TIV"),postvax.b=rep(log2(kdata3.b0$postvax.B.Brisbane[1]),mrow),t=(1:mrow)-1,
                       event=c(rep(0,mrow-1),1*(kdata3.b0$swab.Victoria[1]==1 |kdata3.b0$swab.Yamagata[1]==1 ))     )
kdata.td.b0$proxy <- s.proxy$s.proxy[kdata3.b0$t0[1]:kdata3.b0$T[1]]
kdata.td2.b0 <- kdata.td.b0[nrow(kdata.td.b0),]

for(j in 2:nrow(kdata3.b0)){
   mrow <- kdata3.b0$T[j]-kdata3.b0$t0[j]+1
   temp <- data.frame(hhID=rep(kdata3.b0$hhID[j],mrow),TIV=1*(kdata3.b0$intervention[j]=="TIV"),postvax.b=rep(log2(kdata3.b0$postvax.B.Brisbane[j]),mrow),t=(1:mrow)-1,
                      event=c(rep(0,mrow-1),1*(kdata3.b0$swab.Victoria[j]==1 |kdata3.b0$swab.Yamagata[j]==1 ))  )
   temp$proxy <- s.proxy$s.proxy[kdata3.b0$t0[j]:kdata3.b0$T[j]]
   temp2 <- temp[nrow(temp),]
   kdata.td.b0 <- rbind(kdata.td.b0,temp)
   kdata.td2.b0 <- rbind(kdata.td2.b0,temp2)
}

kdata.td.b <-rbind(kdata.td.b0, kdata.td.b1)
cox.td.B <- coxph(Surv(t,event)~postvax.b+proxy,data=kdata.td.b)

# plot

windows(width=11.6,height=12)
layout(matrix(1:6, nrow=3, byrow=T), widths=rep(5.8,6), heights=rep(4,6))
par(mar=c(4,0,0,0), oma=c(0,0,0,0), family="sans")
                           
#  Waning sen {upper 2 panels) 

x<-78:500

d <-c(31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31)
m <-c("A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D")

#### exclude 30%
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_5.r")

#ph1
plot(0,0,type="n",xlab="",ylab="", xlim=c(0-65,540),ylim=c(log2(2.5)-0.6, log2(10240)), axes=F, main="")

points(datr.p1$calpostv, jitter(log2(datr.p1$postvax.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.p1$calmid, jitter(log2(datr.p1$mid.season.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.p1$calposts, jitter(log2(datr.p1$post.season.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.p0$calpostv, jitter(log2(datr.p0$postvax.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.p0$calmid, jitter(log2(datr.p0$mid.season.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.p0$calposts, jitter(log2(datr.p0$post.season.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)

lines(x, model.p0.30$coef[1]+x*model.p0.30$coef[2])
lines(x, model.p1.30$coef[1]+x*model.p1.30$coef[2] , col="red")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=log2(2.5), tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=log2(2.5), tcl=-0.6)
axis(2,at=c(log2(5*2^(-1:11))), labels=c("","<10",5*2^(1:11)), las=1, pos=0)
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=-0.7, font=1,cex=1, side=1)
mtext("A)", side=2, font=1, line=-1.4, at=log2(10240), las=1)
mtext("A(H1N1)pdm09", side=3, font=1, line=-1.2, at=540/2)
mtext("HAI titer", side=2, line=-1.5)

# Victoria

plot(0,0,type="n",xlab="",ylab="", xlim=c(0-65,540),ylim=c(log2(2.5)-0.6, log2(10240)), axes=F, main="")

points(datr.b1$calpostv, jitter(log2(datr.b1$postvax.B.Brisbane),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calmid, jitter(log2(datr.b1$mid.season.B.Brisbane),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calposts, jitter(log2(datr.b1$post.season.B.Brisbane),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.b0$calpostv, jitter(log2(datr.b0$postvax.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.b0$calmid, jitter(log2(datr.b0$mid.season.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.b0$calposts, jitter(log2(datr.b0$post.season.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)

lines(x, model.b0.30$coef[1]+x*model.b0.30$coef[2])
lines(x, model.b1.30.age$coef[1]+x*model.b1.30.age$coef[2]+model.b1.30.age$coef[3]+x*model.b1.30.age$coef[4] , col="red", lty=2)
lines(x, model.b1.30.age$coef[1]+x*model.b1.30.age$coef[2] , col="red", lty=6)
text(350,log2(80),"9-17y", cex=1, col="red")
text(350,log2(50),"6-8y",cex=1, col="red")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=log2(2.5), tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=log2(2.5), tcl=-0.6)
axis(2,at=c(log2(5*2^(-1:11))), labels=c("","<10",5*2^(1:11)), las=1, pos=0)
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=-0.7, font=1,cex=1, side=1)
mtext("B)", side=2, font=1, line=-1.4, at=log2(10240), las=1)
mtext("B(Victoria-lineage)", side=3, font=1, line=-1.2, at=540/2)


################################### exclude pcr only ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_6.r")

#ph1
plot(0,0,type="n",xlab="",ylab="", xlim=c(0-65,540),ylim=c(log2(2.5)-0.6, log2(10240)), axes=F, main="")

points(datr.p1$calpostv, jitter(log2(datr.p1$postvax.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.p1$calmid, jitter(log2(datr.p1$mid.season.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.p1$calposts, jitter(log2(datr.p1$post.season.pH1),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.p0$calpostv, jitter(log2(datr.p0$postvax.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.p0$calmid, jitter(log2(datr.p0$mid.season.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.p0$calposts, jitter(log2(datr.p0$post.season.pH1),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)

lines(x, model.p0.xpcr$coef[1]+x*model.p0.xpcr$coef[2])
lines(x, model.p1.xpcr$coef[1]+x*model.p1.xpcr$coef[2] , col="red")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=log2(2.5), tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=log2(2.5), tcl=-0.6)
axis(2,at=c(log2(5*2^(-1:11))), labels=c("","<10",5*2^(1:11)), las=1, pos=0)
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=-0.7, font=1,cex=1, side=1)
mtext("C)", side=2, font=1, line=-1.4, at=log2(10240), las=1)
mtext("A(H1N1)pdm09", side=3, font=1, line=-1.4, at=540/2)
mtext("HAI titer", side=2, line=-1.5)

#Victoria
plot(0,0,type="n",xlab="",ylab="", xlim=c(0-65,540),ylim=c(log2(2.5)-0.6, log2(10240)), axes=F, main="")

points(datr.b1$calpostv, jitter(log2(datr.b1$postvax.B.Brisbane),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calmid, jitter(log2(datr.b1$mid.season.B.Brisbane),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)
points(datr.b1$calposts, jitter(log2(datr.b1$post.season.B.Brisbane),2), col=rgb(0.8,0,0.2, alpha=0.4), cex=0.8, pch=16)

points(datr.b0$calpostv, jitter(log2(datr.b0$postvax.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.b0$calmid, jitter(log2(datr.b0$mid.season.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)
points(datr.b0$calposts, jitter(log2(datr.b0$post.season.B.Brisbane),2), col=rgb(0.1,0.1,0.1, alpha=0.4), cex=0.8, pch=16)

lines(x, model.b0.xpcr$coef[1]+x*model.b0.xpcr$coef[2])
lines(x, model.b1.xpcr.age$coef[1]+x*model.b1.xpcr.age$coef[2]+model.b1.xpcr.age$coef[3]+x*model.b1.xpcr.age$coef[4] , col="red", lty=2)
lines(x, model.b1.xpcr.age$coef[1]+x*model.b1.xpcr.age$coef[2] , col="red", lty=6)
text(350,log2(55),"9-17y", cex=1, col="red")
text(350,log2(24),"6-8y",cex=1, col="red")

axis(1,at=c(0,cumsum(d)), labels=rep("",18), pos=log2(2.5), tcl=-0.3)
axis(1,at=c(cumsum(d)[c(1,5,9,13,17)]), labels=rep("",5), pos=log2(2.5), tcl=-0.6)
axis(2,at=c(log2(5*2^(-1:11))), labels=c("","<10",5*2^(1:11)), las=1, pos=0)
mtext(c("1 Sep 2009","1 Jan 2010", "1 May 2010", "1 Sep 2010", "1 Jan 2011") ,at=cumsum(d)[c(1,5,9,13,17)], line=-0.7, font=1,cex=1, side=1)
mtext("D)", side=2, font=1, line=-1.4, at=log2(10240), las=1)
mtext("B(Victoria-lineage)", side=3, font=1, line=-1.2, at=540/2)

#                       
# HIP sen (lower panel) 

#plot for pH1
beta <-cox.td.A$coef
se <-sqrt(sum(vcov(cox.td.A)))

est <-data.frame(titer=log2(5*2^(0:9)),proxy=2, mu=NA, ci1=NA, ci2=NA)
est$mu <-exp(beta[1] * est$titer + beta[2]*est$proxy  )
est$ci1 <-exp(beta[1] * est$titer + beta[2]*est$proxy  - 1.96*se)
est$ci2 <-exp(beta[1] * est$titer + beta[2]*est$proxy + 1.96*se)
est$d <-(est$mu[1]-est$mu)/est$mu[1]
est$d.ci1 <-(est$mu[1]-est$ci2)/est$mu[1]
est$d.ci2 <-(est$mu[1]-est$ci1)/est$mu[1]

plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(log2(5)-1,log2(2560)), ylim=c(0,1))
axis(1,at=log2(5*2^(0:9)), labels=5*2^(0:9))
axis(2,at=0.2*0:5, labels=0.2*0:5, las=1, pos=log2(5))

lines(est$titer,est$d)
lines(est$titer,est$d.ci1, lty=2)
lines(est$titer,est$d.ci2, lty=2)

mtext("HAI titer", side=1, line=3, at=(log2(80)+log2(160))/2)
mtext("Relative risk reduction", side=2, line=-1.5)
mtext("A(H1N1)pdm09",side=3,line=-1.2, font=1, at=(log2(80)+log2(160))/2)
mtext("E)", side=2, line=-1.5, at=1, las=1)

#plot for pH1 (main)
beta <-cox.td.A.m $coef
se <-sqrt(sum(vcov(cox.td.A.m)))

est.pm <-data.frame(titer=log2(5*2^(0:9)),proxy=2, mu=NA, ci1=NA, ci2=NA)
est.pm$mu <-exp(beta[1] * est.pm$titer + beta[2]*est.pm$proxy  )
est.pm$ci1 <-exp(beta[1] * est.pm$titer + beta[2]*est.pm$proxy  - 1.96*se)
est.pm$ci2 <-exp(beta[1] * est.pm$titer + beta[2]*est.pm$proxy + 1.96*se)
est.pm$d <-(est.pm$mu[1]-est.pm$mu)/est.pm$mu[1]
est.pm$d.ci1 <-(est.pm$mu[1]-est.pm$ci2)/est.pm$mu[1]
est.pm$d.ci2 <-(est.pm$mu[1]-est.pm$ci1)/est.pm$mu[1]

lines(est.pm$titer,est.pm$d, col=gray(0.7))
lines(est.pm$titer,est.pm$d.ci1, lty=2, col=gray(0.7))
lines(est.pm$titer,est.pm$d.ci2, lty=2, col=gray(0.7))

mtext("Relative risk reduction", side=2, line=2)

#plot for flu B  (TIV group)
beta <-cox.td.B$coef
se <-sqrt(sum(vcov(cox.td.B)))

est <-data.frame(titer=log2(5*2^(0:9)),proxy=2, mu=NA, ci1=NA, ci2=NA)
est$mu <-exp(beta[1] * est$titer + beta[2]*est$proxy  )
est$ci1 <-exp(beta[1] * est$titer + beta[2]*est$proxy  - 1.96*se)
est$ci2 <-exp(beta[1] * est$titer + beta[2]*est$proxy + 1.96*se)
est$d <-(est$mu[1]-est$mu)/est$mu[1]
est$d.ci1 <-(est$mu[1]-est$ci2)/est$mu[1]
est$d.ci2 <-(est$mu[1]-est$ci1)/est$mu[1]

plot(0,0,xlab="",ylab="", main="", axes=F, xlim=c(log2(5)-1,log2(2560)), ylim=c(0,1))
axis(1,at=log2(5*2^(0:9)), labels=5*2^(0:9))
axis(2,at=0.2*0:5, labels=0.2*0:5, las=1, pos=log2(5))

lines(est$titer,est$d)
lines(est$titer,est$d.ci1, lty=2)
lines(est$titer,est$d.ci2, lty=2)

#plot for flu B  (main)
beta <-vic.vic$coef
se <-sqrt(sum(vcov(vic.vic)))
                                                       
est.bm <-data.frame(titer=log2(5*2^(0:9)),proxy=2, mu=NA, ci1=NA, ci2=NA)
est.bm$mu <-exp(beta[1] * est.bm$titer + beta[2]*est.bm$proxy  )
est.bm$ci1 <-exp(beta[1] * est.bm$titer + beta[2]*est.bm$proxy - 1.96*se)
est.bm$ci2 <-exp(beta[1] * est.bm$titer + beta[2]*est.bm$proxy + 1.96*se)
est.bm$d <-(est.bm$mu[1]-est.bm$mu)/est.bm$mu[1]
est.bm$d.ci1 <-(est.bm$mu[1]-est.bm$ci2)/est.bm$mu[1]
est.bm$d.ci2 <-(est.bm$mu[1]-est.bm$ci1)/est.bm$mu[1]

lines(est.bm$titer,est.bm$d, col=gray(0.7))
lines(est.bm$titer,est.bm$d.ci1, lty=2, col=gray(0.7))
lines(est.bm$titer,est.bm$d.ci2, lty=2, col=gray(0.7))

mtext("HAI titer", side=1, line=3, at=(log2(80)+log2(160))/2)
mtext("B(Victoria-lineage)",side=3,line=-1.2, font=1, at=(log2(80)+log2(160))/2)
mtext("F)", side=2, line=-1.5, at=1, las=1)


#
# End of script.
#



