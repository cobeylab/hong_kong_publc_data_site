#
# R syntax to reproduce information for Table 2 from:
#
# Cowling BJ, Ip DKM, Fang VJ, et al.
# Modes of transmission of influenza B virus in households
# PLoS ONE, 2014 (in press).
#
# Last updated by ang VJ and Cowling BJ.
# August 20, 2014
#

source("http://www.hku.hk/bcowling/influenza/NPI_PLoSOne_motB_scripts/dataframe.r")

demog1 <- read.csv("http://sph.hku.hk/data/HongKongNPIstudyV4/adherence_m.csv")
hchar1 <- read.csv("http://sph.hku.hk/data/HongKongNPIstudyV4/hchar_h.csv")
demog2 <- read.csv("http://sph.hku.hk/data/HongKongNPIstudy2009V1/demog_m.csv")
hchar2 <- read.csv("http://sph.hku.hk/data/HongKongNPIstudy2009V1/hchar_h.csv")

names(demog1)[5] <- "vaccine"; names(demog2)[5] <- "vaccine"
hchar <- rbind(hchar1[c("hhID","intervention","familysize")],hchar2[c("hhID","intervention","familysize")])
id <- unique(hkdata$hhID)
cdata <- rbind(demog1[c("hhID","member","age","male","vaccine")],demog2[c("hhID","member","age","male","vaccine")])
cdata <- cdata[cdata$hhID%in%id,]
hchar <- hchar[hchar$hhID%in%id,]

tab <- matrix(rep(NA, 14*6), nrow=14, dimnames=list(c("Index","Iage<=5","Iage6-15","Iage>16","male","hhsize","contact","age<=5","age6-15","age16-30","age31-50","age>50","male","vaccine"),c("Control","%","HH","%","Mask+HH","%")),byrow=T)

tab[1,c(1,3,5)] <- table(hchar$intervention)
index <- merge(cdata[cdata$member==0,],hchar,by="hhID")
index$agegp <- cut(index$age,c(0,5,15,100))
tab[2:4,c(1,3,5)] <- table(index$agegp,index$intervention)
tab[5,c(1,3,5)] <- table(index$male,index$intervention)[2,]
tab[6,1] <- median(index$familysize[index$intervention==1])
tab[6,3] <- median(index$familysize[index$intervention==3])
tab[6,5] <- median(index$familysize[index$intervention==4])
for(i in 1:3){tab[2:5,i*2] <- round(tab[2:5,i*2-1]/tab[1,i*2-1],2)}

contact <- merge(cdata[cdata$member!=0,],hchar,by="hhID")
tab[7,c(1,3,5)] <- table(contact$intervention)
contact$agegp <- cut(contact$age,c(-0.1,5,15,30,50,100))
tab[8:12,c(1,3,5)] <- table(contact$agegp,contact$intervention)
tab[13,c(1,3,5)] <- table(contact$male,contact$intervention)[2,]
tab[14,c(1,3,5)] <- table(contact$vaccine,contact$intervention)[2,]
for(i in 1:3){tab[8:14,i*2] <- round(tab[8:14,i*2-1]/tab[7,i*2-1],2)}

#
# End of script.
#
