#
# Script to reproduce information in analyzed households (n=259) from:
#
# Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, et al.
# Facemasks and hand hygiene to prevent influenza transmission 
# in households, a randomized trial.
# Annals of Internal Medicine, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Sep 08, 2009

# Analyzed households: no co-index, index confirmed +ve at baseline

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV3/"

hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
housechar <- read.csv(paste(dir, "hchar_h.csv", sep=""))                       
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""))      

# Households followed up and the family members

mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
hculture <- data.frame(hhID = baseflu$hhID, member = baseflu$member)

hc.temp <- reshape(hc, timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="PCR")
hculture <- merge(hculture,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hculture) <- c("hhID","member","V1","V2","V3")

## exd_index: none of V0/V1 culture is A/B; Contact exclusion: V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (!is.na(hculture$V1[i]) & hculture$V1[i]==0)|(is.na(hculture$V1[i])) )) {hculture$exd_index[i]<-1}
       else {hculture$exd_index[i]<-0}
     if(hculture$member[i]!=0 &  !is.na(hculture$V1[i]) & hculture$V1[i]!=0)
       {hculture$d_contact[i]<-1}
       else {hculture$d_contact[i]<-0}
     if(hculture$member[i]!=0 & ( (!is.na(hculture$V1[i]) & hculture$V1[i]!=0) | is.na(hculture$V1[i]) ))
          {hculture$exd_contact[i]=1}
	  else{hculture$exd_contact[i]=0}
}

d_contactid <- unique(hculture$hhID[hculture$d_contact==1]) # for excluding co-index households
d_contact <- data.frame(hhID=d_contactid)
d_contact$d_contact <- 1

exd_index <- hculture[hculture$member==0,c(1,6)]

dim(hculture)
hculture <- merge(hculture[,-6], exd_index)
hculture <- merge(hculture[,-6], d_contact,all.x=TRUE)
hculture$d_contact[is.na(hculture$d_contact)] <- 0
dim(hculture)
hculture <- hculture[order(hculture$hhID,hculture$member),]

hculture$analyzed <- 1*(hculture$exd_index==0&hculture$d_contact==0)

# End of script
