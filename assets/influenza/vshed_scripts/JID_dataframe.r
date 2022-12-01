
#
# R syntax to reproduce information for data set from:
#
# Lau LLH, Cowling BJ, Fang VJ, Chan KH, Lau EHY, et al.
# Viral shedding and clinical illness in naturally acquired influenza virus infections
# JID, 2010 (in press).
#
# Last updated by Vicky Fang and Lincoln Lau
# Jan 30, 2010

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"

hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))                    
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""))      
symp <- read.csv(paste(dir, "symptomday_d.csv", sep=""))

# Households followed up and the family members

mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
hculture <- data.frame(hhID = baseflu$hhID, member = baseflu$member)

hc.temp <- reshape(hc[1:4], timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="qPCR")
hculture <- merge(hculture,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hculture) <- c("hhID","member","V1","V2","V3")

contact <- hculture[hculture$member>0,]

## Define lab-confirmed infection

for (i in 1:nrow(contact)){
    if ( (contact$V1[i] >0 & !is.na(contact$V1[i])) | 
         (contact$V2[i] >0 & !is.na(contact$V2[i])) | 
         (contact$V3[i] >0 & !is.na(contact$V3[i])) )   contact$infect[i] <- 1
	  else contact$infect[i] <- 0
}

# Define the subject with the first swabs collected were PCR+
contact$pos1 <- 1*( (!is.na(contact$V1)&contact$V1>0) | (is.na(contact$V1)&!is.na(contact$V2)&contact$V2>0) |
                    (is.na(contact$V1)&is.na(contact$V2)&!is.na(contact$V3)&contact$V3>0))

# Define ARI
symp <- data.frame(lapply(symp,function(x,...){x[is.na(x)] <- 0 ; x}))
symp$fever <- 1*(symp$bodytemp>=37.8)
symp$ARI <- 1*((symp$fever+symp$headache+symp$sthroat+symp$cough+symp$pmuscle+symp$rnose+symp$phlegm)>=2)
contact <- merge(contact,symp[symp$day==0,c("hhID","member","ARI")],by=c("hhID","member"),all.x=T)

# 59 subjects retained in the analysis
contact <- contact[contact$infect==1,]
contact$analyse <- 1*(contact$pos1==0&contact$ARI==0)
# table(contact$analyse)

# 44 subjects for plot (exclude those without ARI onset)
symp2 <- reshape(symp[c(1:3,12)],timevar="day", idvar=c("hhID","member"), direction="wide", v.names="ARI")
symp2$onset <- 1*((symp2$ARI.0+symp2$ARI.1+symp2$ARI.2+symp2$ARI.3+symp2$ARI.4+
                  symp2$ARI.5+symp2$ARI.6+symp2$ARI.7+symp2$ARI.8+symp2$ARI.9)>0)

contact <- merge(contact,symp2[c("hhID","member","onset")],by=c("hhID","member"),all.x=T)
contact$plot <- 1*(contact$analyse==1&contact$onset==1)
# table(contact$onset)

# make a new data frame for the 44 subjects
# ARI onset

symp2$any_onset <- 0
symp2$any_onset[symp2$ARI.0 == 1] <- 0
symp2$any_onset[symp2$ARI.0 == 0 & symp2$ARI.1 == 1] <- 1
symp2$any_onset[symp2$ARI.0 == 0 & symp2$ARI.1 == 0 & symp2$ARI.2 == 1] <- 2
symp2$any_onset[symp2$ARI.0 == 0 & symp2$ARI.1 == 0 & symp2$ARI.2 == 0
          &symp2$ARI.3 == 1 ] <- 3
symp2$any_onset[symp2$ARI.0 == 0 & symp2$ARI.1 == 0 & symp2$ARI.2 == 0
          &symp2$ARI.3 == 0 & symp2$ARI.4 == 1] <- 4
symp2$any_onset[symp2$ARI.0 == 0 & symp2$ARI.1 == 0 & symp2$ARI.2 == 0
          &symp2$ARI.3 == 0 & symp2$ARI.4 == 0 & symp2$ARI.5 == 1] <- 5
symp2$any_onset[symp2$ARI.0 == 0 & symp2$ARI.1 == 0 & symp2$ARI.2 == 0
          &symp2$ARI.3 == 0 & symp2$ARI.4 == 0 & symp2$ARI.5 == 0 & symp2$ARI.6 == 1] <- 6

plotdata <- merge(contact[contact$plot==1,1:5],symp2[c(1,2,14)], by=c("hhID", "member"),all.x=T)

hc2 <- reshape(hc[c(1:3,5)], timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="q_culture") 
names(hc2) <- c("hhID","member","T1","T2","T3")
plotdata <- merge(plotdata,hc2,by=c("hhID","member"),all.x=T)

plotdata <- plotdata[order(plotdata$hhID,plotdata$member),]

#
# End of script
#

