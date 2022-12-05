#
# R syntax to reproduce Appendix Table 1 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

dir <- "http://sph.hku.hk/data/KiddivaxMainV2/"
kdata <- read.csv(paste(dir, "serology.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))
demog <- read.csv(paste(dir, "demographic.csv", sep=""))
chron <- read.csv(paste(dir, "chronic_condition.csv", sep=""))

##
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

dat1 <-merge(kdata3, demog,  by=c("hhID", "member"), all.x=T)
dat2 <-merge(dat1,chron,  by=c("hhID", "member"), all.x=T)

m <-matrix(rep(NA,7*11), ncol=7, dimnames=list(c("age<8", "age>8",NA, "male", "female", NA,"v0809-yes", "v0809-no", NA,"chron-yes", "chron-no"),c("n1","%",NA,"n0","%",NA,"p")))

m[1,c(1,4)]<-table(dat2$age<8, dat2$intervention)[2,2:1]
m[2,c(1,4)]<-table(dat2$age<8, dat2$intervention)[1,2:1]
m[1,7] <-round(chisq.test(table(dat2$age<8, dat2$intervention))$p.value,2)

m[4,c(1,4)]<-table(dat2$male, dat2$intervention)[2,2:1]
m[5,c(1,4)]<-table(dat2$male, dat2$intervention)[1,2:1]
m[4,7] <-round(chisq.test(table(dat2$male, dat2$intervention))$p.value,2)

m[7,c(1,4)]<-table(dat2$vac0809, dat2$intervention)[2,2:1]
m[8,c(1,4)]<-table(dat2$vac0809, dat2$intervention)[1,2:1]
m[7,7] <-round(chisq.test(table(dat2$vac0809, dat2$intervention))$p.value,2)

m[10,c(1,4)]<-table(dat2$asthma==1 | dat2$aller.inflam.immun==1 | dat2$cardvas==1 | dat2$resp==1 | dat2$endo==1 | dat2$canc==1 | dat2$haem==1 | dat2$ren==1, dat2$intervention)[2,2:1]
m[11,c(1,4)]<-table(dat2$asthma==1 | dat2$aller.inflam.immun==1 | dat2$cardvas==1 | dat2$resp==1 | dat2$endo==1 | dat2$canc==1 | dat2$haem==1 | dat2$ren==1, dat2$intervention)[1,2:1]
m[10,7] <-round(chisq.test(table(dat2$asthma==1 | dat2$aller.inflam.immun==1 | dat2$cardvas==1 | dat2$resp==1 | dat2$endo==1 | dat2$canc==1 | dat2$haem==1 | dat2$ren==1, dat2$intervention))$p.value,2)

m[,2] <-round(m[,1]/nrow(dat2[dat2$intervention=="TIV",]),2)
m[,5] <-round(m[,4]/nrow(dat2[dat2$intervention=="placebo",]),2)

m2<-rbind(matrix(c(nrow(dat2[dat2$randgp==1,]),NA,NA,nrow(dat2[dat2$randgp==0,]),NA,NA,NA),nrow=1),m)
m2

#
# End of script.
#
