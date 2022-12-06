#
# R syntax to reformat raw data for:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Humoral antibody response after receipt of inactivated seasonal
# influenza vaccinations one year apart in children
# PIDJ, 2012.
#
# Last updated by Ng S, Fang VJ and Cowling BJ.
# February 15, 2013
#

dir1 <- "http://sph.hku.hk/data/KiddivaxPilotV2/"
dir2 <- "http://sph.hku.hk/data/KiddivaxMainV1/"
sero1 <- read.csv(paste(dir1, "serology_m.csv", sep=""))
random1 <- read.csv(paste(dir1, "randomcode_h.csv", sep=""))
symp1 <- read.csv(paste(dir1, "symptom_d.csv", sep=""))
sero2 <- read.csv(paste(dir2, "serology.csv", sep=""))
random2 <- read.csv(paste(dir2, "randomcode.csv", sep=""))
demog2 <- read.csv(paste(dir2, "demographic.csv", sep=""))

id <- sero2$pilot.hhID[sero2$member==0&sero2$pilot.member==0&!is.na(sero2$pilot.member)]
sero1 <- sero1[sero1$hhID%in%id&sero1$member==0,]; sero1 <- merge(sero1,random1,all.x=T)
sero2 <- sero2[sero2$pilot.hhID%in%id&sero2$member==0,]; sero2 <- merge(sero2,random2,all.x=T)

#format data
sero2$aftervactime <- as.numeric(as.Date(sero2$postvax,"%d/%m/%Y")-as.Date(sero2$prevax,"%d/%m/%Y"))
sero2$aftervactime[!is.na(sero2$vaxdate)] <- as.numeric(as.Date(sero2$postvax[!is.na(sero2$vaxdate)],"%d/%m/%Y")-as.Date(sero2$vaxdate[!is.na(sero2$vaxdate)],"%d/%m/%Y"))
sero2[sero2$aftervactime>80,c("postvax.pH1","postvax.sH1","postvax.sH3","postvax.B.Brisbane","postvax.B.Floride")] <- NA

names(sero1) <- paste("pilot.",names(sero1),sep="")
sero <- merge(sero1,sero2)
sero <- merge(sero,demog2[c("hhID","age","male")],all.x=T)

# add serology confirmed infections (>=4 fold rise)
sero$pilot.ph1s <- 1*(as.numeric(as.character(sero$pilot.ph1.post.titer))/as.numeric(as.character(sero$pilot.ph1.mid.titer))>=4)
sero$pilot.sh1w <- 1*(as.numeric(as.character(sero$pilot.sh1.mids))/as.numeric(as.character(sero$pilot.sh1.postv))>=4)
sero$pilot.sh1s <- 1*(as.numeric(as.character(sero$pilot.sh1.posts))/as.numeric(as.character(sero$pilot.sh1.mids))>=4)
sero$pilot.sh3w <- 1*(as.numeric(as.character(sero$pilot.sh3.mids))/as.numeric(as.character(sero$pilot.sh3.postv))>=4)
sero$pilot.sh3s <- 1*(as.numeric(as.character(sero$pilot.sh3.posts))/as.numeric(as.character(sero$pilot.sh3.mids))>=4)
sero$pilot.ph1s[sero$pilot.hhID==9211] <- 0 # adjust for cross-reactivity
sero$pilot.sh1w[sero$pilot.hhID==9164] <- 1 # adjust from RT-PCR result
sero$pilot.sh1 <- 1*(sero$pilot.sh1w==1|sero$pilot.sh1s==1); sero$pilot.sh3 <- 1*(sero$pilot.sh3w==1|sero$pilot.sh3s==1)

symp1$ARI <- 1*((1*(symp1$bodytemp>=37.8)+symp1$headache+symp1$sthroat+symp1$cough+symp1$pmuscle+symp1$phlegm+symp1$rnose)>=2)
symp1$ILI <- 1*(symp1$bodytemp>=37.8&(symp1$cough==1|symp1$sthroat==1))
sero$pilot.ari <- 1*(sero$pilot.hhID%in%unique(symp1$hhID[symp1$ARI==1&symp1$member==0])); sero$pilot.ari[sero$pilot.hhID%in%c(9135,9170,9180)] <- NA
sero$pilot.ili <- 1*(sero$pilot.hhID%in%unique(symp1$hhID[symp1$ILI==1&symp1$member==0])); sero$pilot.ili[sero$pilot.hhID%in%c(9135,9170,9180)] <- NA

#
# End of script.
#

