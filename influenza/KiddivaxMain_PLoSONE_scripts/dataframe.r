#
# R syntax to reformat raw data for:
#
# Ng S, Ip DKM, Fang VJ, et al.
# The effect of age and recent influenza vaccination history on the immunogenicity and efficacy 
# of 2009-10 seasonal trivalent inactivated influenza vaccination in children
# PLoS ONE 2013 (in press).
#
# Last updated by Ng S, Fang VJ and Cowling BJ.
# March 8, 2013
#

require(chron)
dir1 <- "http://sph.hku.hk/data/KiddivaxPilotV2/"
random1 <- read.csv(paste(dir1, "randomcode_h.csv", sep="")); names(random1) <- c("pilot.hhID","pilot.intervention");random1$pilot.member <- 0
dir2 <- "http://sph.hku.hk/data/KiddivaxMainV2/"
chronic2 <- read.csv(paste(dir2, "chronic_condition.csv", sep=""))
random2 <- read.csv(paste(dir2, "randomcode.csv", sep=""))
swab <- read.csv(paste(dir2, "swab.csv", sep=""))
demog <- read.csv(paste(dir2, "demographic.csv", sep=""))
symp <- read.csv(paste(dir2, "symptom_d.csv", sep=""))
sero <- read.csv(paste(dir2, "serology.csv", sep=""))

# Define RT-PCR (swab) confirmed infections
swab <- merge(swab,sero[c(1,2,5,6)],by=c("hhID","member"),all.x=T)
swab$mark <- 1*(as.Date(as.character(swab$date),format="%Y%m%d")-as.Date(as.character(swab$start.date),format="%d/%m/%Y")>=0&
                as.Date(as.character(swab$date),format="%Y%m%d")-as.Date(as.character(swab$end.date),format="%d/%m/%Y")<=0)
swab <- swab[swab$mark==1,1:11]

swab.pH1 <- unique(swab[swab$Swine.H1=="P" ,c("hhID","member")]);
swab.sH3 <- unique(swab[swab$H3=="P",c("hhID","member")]); swab.B <- unique(swab[swab$FluB=="P",c("hhID","member")])
swab.pH1$swab.pH1 <-1; swab.sH3$swab.sH3 <-1; swab.B$swab.B <-1

# Define serology-confirmed infections (>=4 fold rise in antibody titers)
sero$pH1 <- 1*(sero$post.season.pH1/sero$postvax.pH1>=4)
sero$sH3 <- 1*(sero$post.season.sH3/sero$postvax.sH3>=4)
sero$B <- 1*(sero$post.season.B.Brisbane/sero$postvax.B.Brisbane>=4)

# construct the data frame
index.pre <- merge(demog,random,by="hhID",all.x=T) 
index.pre$age[index.pre$age<6] <- 6

# merge in serology results
index.pre <- merge(index.pre,sero,by=c("hhID","member"),all.x=T)

# vac0910.pre=1 if intervention=TIV
index.pre$vac0910.pre[index.pre$intervention=="TIV"]<-1
index.pre$vac0910.pre[index.pre$intervention=="placebo"]<-0

# merge in ARI/ILI
symp <- symp[symp$hhID!=2654,]  # hhID=2654 should be dropout, input error??
symp <- merge(symp,sero[c(1,2,5,6)],by=c("hhID","member"),all.x=T)
symp$mark <- 1*(as.Date(as.character(symp$date),format="%d/%m/%Y")-as.Date(as.character(symp$start.date),format="%d/%m/%Y")>=0&
                as.Date(as.character(symp$date),format="%d/%m/%Y")-as.Date(as.character(symp$end.date),format="%d/%m/%Y")<=0)
symp <- symp[symp$mark==1,1:14]

ARI <- unique(symp[symp$ARI==1,1:2]); ARI$ARI <- 1
ILI <- unique(symp[symp$ILI==1,1:2]); ILI$ILI <- 1
index.pre <- merge(index.pre,ARI, by=c("hhID","member"),all.x=T); index.pre <- merge(index.pre,ILI, by=c("hhID","member"),all.x=T)
index.pre$ARI[is.na(index.pre$ARI)] <- index.pre$ILI[is.na(index.pre$ILI)] <- 0

# merge swab A and B +ve
index.pre <- merge(index.pre,swab.pH1,by=c("hhID","member"),all.x=T)
index.pre <- merge(index.pre,swab.sH3,by=c("hhID","member"),all.x=T)
index.pre <- merge(index.pre,swab.B,by=c("hhID","member"),all.x=T)
index.pre$swab.pH1[is.na(index.pre$swab.pH1)] <- 0
index.pre$swab.sH3[is.na(index.pre$swab.sH3)] <- 0
index.pre$swab.B[is.na(index.pre$swab.B)] <- 0

##
index.pre <- merge(index.pre,chronic2,by=c("hhID","member"),all.x=T)
index.pre <- merge(index.pre,random1,by=c("pilot.hhID","pilot.member"),all.x=T)
index.pre$chron2 <-1*(index.pre$aller.inflam.immun==1|index.pre$canc==1|index.pre$cardvas==1|index.pre$endo==1|index.pre$resp==1)

dat <- index.pre
dat$vac0809[dat$pilot.intervention=="TIV"&!is.na(dat$pilot.intervention)] <-1       # Correct vac0809 by pilot intervention group
dat$vac0809[dat$pilot.intervention=="placebo"&!is.na(dat$pilot.intervention)] <-0

dat$prevax.sH1.c <-log(dat$prevax.sH1/2.5,2)
dat$prevax.sH3.c <-log(dat$prevax.sH3/2.5,2)
dat$prevax.pH1.c <-log(dat$prevax.pH1/2.5,2)
dat$prevax.B.Brisbane.c <-log(dat$prevax.B.Brisbane/2.5,2)
dat$prevax.B.Floride.c <-log(dat$prevax.B.Floride/2.5,2)

dat$postvax.sH1.c <-log(dat$postvax.sH1/2.5,2)
dat$postvax.sH3.c <-log(dat$postvax.sH3/2.5,2)
dat$postvax.pH1.c <-log(dat$postvax.pH1/2.5,2)
dat$postvax.B.Brisbane.c <-log(dat$postvax.B.Brisbane/2.5,2)
dat$postvax.B.Floride.c <-log(dat$postvax.B.Floride/2.5,2)

dat$post.season.sH1.c <-log(dat$post.season.sH1/2.5,2)
dat$post.season.sH3.c <-log(dat$post.season.sH3/2.5,2)
dat$post.season.pH1.c <-log(dat$post.season.pH1/2.5,2)
dat$post.season.B.Brisbane.c <-log(dat$post.season.B.Brisbane/2.5,2)
dat$post.season.B.Floride.c <-log(dat$post.season.B.Floride/2.5,2)

#
# End of script.
#




