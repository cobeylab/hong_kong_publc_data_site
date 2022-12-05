#
# R syntax to generate information used in:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

### pH1 Model
dir <- "http://sph.hku.hk/data/KiddivaxMainV2/"
kdata <- read.csv(paste(dir, "serology.csv", sep=""))
demog <- read.csv(paste(dir, "demographic.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))

time0 <-"31/07/2009"
dat <-kdata[kdata$member==0 & !is.na(kdata$postvax.pH1) & !is.na(kdata$post.season.pH1),]
dat$followup <-as.numeric(as.Date(dat$post.season, format="%d/%m/%Y")-as.Date(dat$postvax, format="%d/%m/%Y"))
dat$followup2 <-as.numeric(as.Date(dat$mid.season, format="%d/%m/%Y")-as.Date(dat$postvax, format="%d/%m/%Y"))
dat$calpre <-as.numeric(as.Date(dat$prevax, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calpostv <-as.numeric(as.Date(dat$postvax, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calposts<-as.numeric(as.Date(dat$post.season, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calmid<-as.numeric(as.Date(dat$mid.season, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$GMTR <-dat$post.season.pH1 / dat$postvax.pH1
dat$GMTR2 <-dat$mid.season.pH1 / dat$postvax.pH1
dat$GMTR3 <-dat$post.season.pH1 / dat$mid.season.pH1

dat<-merge(dat, random, by="hhID", all.x=T)
dat<-merge(dat, demog, by=c("hhID", "member"), all.x=T)

dat2 <-dat
swab2 <-swab[swab$member==0,]

datr.p1  <-dat2[dat2$intervention=="TIV" & !(dat2$hhID %in% swab2$hhID[swab2$Swine.H1=="P"]),]
datr.p0 <-dat2[dat2$intervention=="placebo" & !(dat2$hhID %in% swab2$hhID[swab2$Swine.H1=="P"]),]

datr2a.p1 <-data.frame(hhID=datr.p1$hhID, TIV=1*(datr.p1$intervention=="TIV"), caltime=datr.p1$calpostv ,titer=log2(datr.p1$postvax.pH1))
datr2b.p1 <-data.frame(hhID=datr.p1$hhID, TIV=1*(datr.p1$intervention=="TIV"), caltime=datr.p1$calmid ,titer=log2(datr.p1$mid.season.pH1))
datr2c.p1 <-data.frame(hhID=datr.p1$hhID, TIV=1*(datr.p1$intervention=="TIV"), caltime=datr.p1$calposts ,titer=log2(datr.p1$post.season.pH1))

datr2a.p0 <-data.frame(hhID=datr.p0$hhID, TIV=1*(datr.p0$intervention=="TIV"), caltime=datr.p0$calpostv ,titer=log2(datr.p0$postvax.pH1))
datr2b.p0 <-data.frame(hhID=datr.p0$hhID, TIV=1*(datr.p0$intervention=="TIV"), caltime=datr.p0$calmid ,titer=log2(datr.p0$mid.season.pH1))
datr2c.p0 <-data.frame(hhID=datr.p0$hhID, TIV=1*(datr.p0$intervention=="TIV"), caltime=datr.p0$calposts ,titer=log2(datr.p0$post.season.pH1))

datr2.p1 <-rbind(datr2a.p1,datr2b.p1,datr2c.p1)
datr2.p0 <-rbind(datr2a.p0,datr2b.p0,datr2c.p0)

require(gee)
model.p1.xpcr <-gee(titer~caltime,data=datr2.p1,family=gaussian, id=hhID)
model.p0.xpcr <-gee(titer~caltime,data=datr2.p0,family=gaussian, id=hhID)
model.pi.xpcr <- gee(titer~caltime+TIV+caltime*TIV,data=rbind(datr2.p1, datr2.p0),family=gaussian, id=hhID)

# B model -----------------------------------------------------------------------------------------------------------------------------------------
kdata <- read.csv(paste(dir, "serology.csv", sep=""))
demog <- read.csv(paste(dir, "demographic.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))

time0 <-"31/07/2009"
dat <-kdata[kdata$member==0 & !is.na(kdata$postvax.B.Brisbane) & !is.na(kdata$post.season.B.Brisbane),]
dat$followup <-as.numeric(as.Date(dat$post.season, format="%d/%m/%Y")-as.Date(dat$postvax, format="%d/%m/%Y"))
dat$followup2 <-as.numeric(as.Date(dat$mid.season, format="%d/%m/%Y")-as.Date(dat$postvax, format="%d/%m/%Y"))
dat$calpre <-as.numeric(as.Date(dat$prevax, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calpostv <-as.numeric(as.Date(dat$postvax, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calposts<-as.numeric(as.Date(dat$post.season, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$calmid<-as.numeric(as.Date(dat$mid.season, format="%d/%m/%Y")-as.Date(time0, format="%d/%m/%Y"))
dat$GMTR <-dat$post.season.B.Brisbane / dat$postvax.B.Brisbane
dat$GMTR2 <-dat$mid.season.B.Brisbane / dat$postvax.B.Brisbane
dat$GMTR3 <-dat$post.season.B.Brisbane / dat$mid.season.B.Brisbane

dat<-merge(dat, random, by="hhID", all.x=T)
dat<-merge(dat, demog, by=c("hhID", "member"), all.x=T)

dat2 <-dat
swab2 <-swab[swab$member==0,]

datr.b1 <-dat2[ !(dat2$hhID %in% swab2$hhID[swab2$FluB=="P"]) & dat2$intervention=="TIV",]
datr.b0 <-dat2[ !(dat2$hhID %in% swab2$hhID[swab2$FluB=="P"]) & dat2$intervention=="placebo",]

datr2a.b1 <-data.frame(hhID=datr.b1$hhID, TIV=1*(datr.b1$intervention=="TIV"), age=datr.b1$age, caltime=datr.b1$calpostv ,titer=log2(datr.b1$postvax.B.Brisbane))
datr2b.b1 <-data.frame(hhID=datr.b1$hhID, TIV=1*(datr.b1$intervention=="TIV"), age=datr.b1$age, caltime=datr.b1$calmid ,titer=log2(datr.b1$mid.season.B.Brisbane))
datr2c.b1 <-data.frame(hhID=datr.b1$hhID, TIV=1*(datr.b1$intervention=="TIV"), age=datr.b1$age, caltime=datr.b1$calposts ,titer=log2(datr.b1$post.season.B.Brisbane))

datr2a.b0 <-data.frame(hhID=datr.b0$hhID, TIV=1*(datr.b0$intervention=="TIV"), age=datr.b0$age, caltime=datr.b0$calpostv ,titer=log2(datr.b0$postvax.B.Brisbane))
datr2b.b0 <-data.frame(hhID=datr.b0$hhID, TIV=1*(datr.b0$intervention=="TIV"), age=datr.b0$age, caltime=datr.b0$calmid ,titer=log2(datr.b0$mid.season.B.Brisbane))
datr2c.b0 <-data.frame(hhID=datr.b0$hhID, TIV=1*(datr.b0$intervention=="TIV"), age=datr.b0$age, caltime=datr.b0$calposts ,titer=log2(datr.b0$post.season.B.Brisbane))

datr2.b1 <-rbind(datr2a.b1,datr2b.b1,datr2c.b1)
datr2.b0 <-rbind(datr2a.b0,datr2b.b0,datr2c.b0)

require(gee)
model.b1.xpcr <-gee(titer~caltime,data=datr2.b1,family=gaussian, id=hhID)
model.b0.xpcr <-gee(titer~caltime,data=datr2.b0,family=gaussian, id=hhID)
model.b1.xpcr.age <-gee(titer~caltime+(age<9)+(age<9)*caltime,data=datr2.b1,family=gaussian, id=hhID)
model.bi.xpcr <- gee(titer~caltime+TIV+caltime*TIV,data=rbind(datr2.b1, datr2.b0),family=gaussian, id=hhID)

n.xpcr <- c(nrow(datr.p1), nrow(datr.p0), nrow(datr.b1[datr.b1$age<9,]), nrow(datr.b1[datr.b1$age>=9,]), nrow(datr.b0))

#
# End of script.
#


