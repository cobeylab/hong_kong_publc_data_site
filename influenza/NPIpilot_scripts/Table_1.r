#
# Script to reproduce information in Table 1 (subject characteristics) from:
#
# Cowling BJ, Fung ROP, Cheng KY, Fang VJ, Chan KH, Seto WH, et al.
# Preliminary findings of a randomized trial of non-pharmaceutical
# interventions to prevent influenza transmission in households.
# PLoS ONE, 2008; 3(5):e2101.
# http://dx.doi.org/10.1371/journal.pone.0002101
#
# Last updated by Vicky Fang and Ben Cowling
# January 5, 2009


dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

q1data <- read.csv(paste(dir, "clinicdat_h.csv", sep=""), header=TRUE)
intervention <- read.csv(paste(dir, "randomarm_198.csv", sep=""), header=TRUE)

table1p1 <- matrix(rep(NA,180),ncol=12,byrow=FALSE)

####################################### Part 1-1: For randomized index subjects (n=198) #####################################################

for(i in 1:nrow(q1data)){
   if(!is.na(q1data$bodytemp[i])&q1data$bodytemp[i]>=38) q1data$fever[i] <- 1
   else if(!is.na(q1data$bodytemp[i])&q1data$bodytemp[i]<38) q1data$fever[i] <- 0
   else q1data$fever[i] <- NA
}
q1 <- q1data[,c(1,3,4,28,7:25)]

tab1 <- merge(intervention,q1,by="scrID",all.x=TRUE)

tab1_c <- tab1[tab1$intervention==1,]
tab1_m <- tab1[tab1$intervention==2,]
tab1_h <- tab1[tab1$intervention==3,]

####### Age group #######
table1p1[1:4,1] <- c(dim(tab1_c[tab1_c$age<=15,])[1], dim(tab1_c[tab1_c$age<=30&tab1_c$age>15,])[1],
                     dim(tab1_c[tab1_c$age<=50&tab1_c$age>30,])[1], dim(tab1_c[tab1_c$age>50,])[1])

table1p1[1:4,5] <- c(dim(tab1_m[tab1_m$age<=15,])[1], dim(tab1_m[tab1_m$age<=30&tab1_m$age>15,])[1],
		     dim(tab1_m[tab1_m$age<=50&tab1_m$age>30,])[1], dim(tab1_m[tab1_m$age>50,])[1])

table1p1[1:4,9] <- c(dim(tab1_h[tab1_h$age<=15,])[1], dim(tab1_h[tab1_h$age<=30&tab1_h$age>15,])[1],
		     dim(tab1_h[tab1_h$age<=50&tab1_h$age>30,])[1], dim(tab1_h[tab1_h$age>50,])[1])

####### Sex #######
table1p1[5,1] <- dim(tab1_c[tab1_c$sex=="m",])[1]
table1p1[5,5] <- dim(tab1_m[tab1_m$sex=="m",])[1]
table1p1[5,9] <- dim(tab1_h[tab1_h$sex=="m",])[1]

####### Symptoms #######
table1p1[6:12,1] <- colSums(tab1_c[,c(12,9,22,6,7,11,20)],na.rm=TRUE)
table1p1[6:12,5] <- colSums(tab1_m[,c(12,9,22,6,7,11,20)],na.rm=TRUE)
table1p1[6:12,9] <- colSums(tab1_h[,c(12,9,22,6,7,11,20)],na.rm=TRUE)

####### Delay from symptom onset to randomization #######

table1p1[13:15,1] <- c(dim(tab1_c[tab1_c$onsettime<=2&!is.na(tab1_c$onsettime),])[1], 
                       dim(tab1_c[tab1_c$onsettime>2&tab1_c$onsettime<=4&!is.na(tab1_c$onsettime),])[1],
		       dim(tab1_c[tab1_c$onsettime>4&tab1_c$onsettime!=9&!is.na(tab1_c$onsettime),])[1])

table1p1[13:15,5] <- c(dim(tab1_m[tab1_m$onsettime<=2&!is.na(tab1_m$onsettime),])[1], 
                       dim(tab1_m[tab1_m$onsettime>2&tab1_m$onsettime<=4&!is.na(tab1_m$onsettime),])[1],
		       dim(tab1_m[tab1_m$onsettime>4&tab1_m$onsettime!=9&!is.na(tab1_m$onsettime),])[1])

table1p1[13:15,9] <- c(dim(tab1_h[tab1_h$onsettime<=2&!is.na(tab1_h$onsettime),])[1], 
                       dim(tab1_h[tab1_h$onsettime>2&tab1_h$onsettime<=4&!is.na(tab1_h$onsettime),])[1],
		       dim(tab1_h[tab1_h$onsettime>4&tab1_h$onsettime!=9&!is.na(tab1_h$onsettime),])[1])

table1p1[,2] <- round(table1p1[,1]/nrow(tab1_c),2)
table1p1[,6] <- round(table1p1[,5]/nrow(tab1_m),2)
table1p1[,10] <- round(table1p1[,9]/nrow(tab1_h),2)

####################################### Part 1-2: For followed up index subjects (n=128) ###################################################

hchar <-  read.csv(paste(dir, "hchar_h.csv", sep=""), header=TRUE)

hchar <- hchar[,1:2] # hhID, intervention
tab1.2 <- merge(tab1,hchar,by=c("hhID","intervention")) # Only include questionnaire records for 128 randomized index subjects

tab1_c <- tab1.2[tab1.2$intervention==1,]
tab1_m <- tab1.2[tab1.2$intervention==2,]
tab1_h <- tab1.2[tab1.2$intervention==3,]

####### Age group #######
table1p1[1:4,3] <- c(dim(tab1_c[tab1_c$age<=15,])[1], dim(tab1_c[tab1_c$age<=30&tab1_c$age>15,])[1],
                     dim(tab1_c[tab1_c$age<=50&tab1_c$age>30,])[1], dim(tab1_c[tab1_c$age>50,])[1])

table1p1[1:4,7] <- c(dim(tab1_m[tab1_m$age<=15,])[1], dim(tab1_m[tab1_m$age<=30&tab1_m$age>15,])[1],
		     dim(tab1_m[tab1_m$age<=50&tab1_m$age>30,])[1], dim(tab1_m[tab1_m$age>50,])[1])

table1p1[1:4,11] <- c(dim(tab1_h[tab1_h$age<=15,])[1], dim(tab1_h[tab1_h$age<=30&tab1_h$age>15,])[1],
		     dim(tab1_h[tab1_h$age<=50&tab1_h$age>30,])[1], dim(tab1_h[tab1_h$age>50,])[1])

####### Sex #######
table1p1[5,3] <- dim(tab1_c[tab1_c$sex=="m",])[1]
table1p1[5,7] <- dim(tab1_m[tab1_m$sex=="m",])[1]
table1p1[5,11] <- dim(tab1_h[tab1_h$sex=="m",])[1]

####### Symptoms #######
table1p1[6:12,3] <- colSums(tab1_c[,c(12,9,22,6,7,11,20)],na.rm=TRUE)
table1p1[6:12,7] <- colSums(tab1_m[,c(12,9,22,6,7,11,20)],na.rm=TRUE)
table1p1[6:12,11] <- colSums(tab1_h[,c(12,9,22,6,7,11,20)],na.rm=TRUE)

####### Delay from symptom onset to randomization #######

table1p1[13:15,3] <- c(dim(tab1_c[tab1_c$onsettime<=2&!is.na(tab1_c$onsettime),])[1], 
                       dim(tab1_c[tab1_c$onsettime>2&tab1_c$onsettime<=4&!is.na(tab1_c$onsettime),])[1],
		       dim(tab1_c[tab1_c$onsettime>4&tab1_c$onsettime!=9&!is.na(tab1_c$onsettime),])[1])

table1p1[13:15,7] <- c(dim(tab1_m[tab1_m$onsettime<=2&!is.na(tab1_m$onsettime),])[1], 
                       dim(tab1_m[tab1_m$onsettime>2&tab1_m$onsettime<=4&!is.na(tab1_m$onsettime),])[1],
		       dim(tab1_m[tab1_m$onsettime>4&tab1_m$onsettime!=9&!is.na(tab1_m$onsettime),])[1])

table1p1[13:15,11] <- c(dim(tab1_h[tab1_h$onsettime<=2&!is.na(tab1_h$onsettime),])[1], 
                       dim(tab1_h[tab1_h$onsettime>2&tab1_h$onsettime<=4&!is.na(tab1_h$onsettime),])[1],
		       dim(tab1_h[tab1_h$onsettime>4&tab1_h$onsettime!=9&!is.na(tab1_h$onsettime),])[1])

table1p1[,4] <- round(table1p1[,3]/nrow(tab1_c),2)
table1p1[,8] <- round(table1p1[,7]/nrow(tab1_m),2)
table1p1[,12] <- round(table1p1[,11]/nrow(tab1_h),2)

rownames(table1p1) <- c("2-15yr","16-30yr","31-50yr","50+yr","men","cough","rnose","tired","fever","headache","sthroat","pmuscle","0-24hr","24-48hr","48+hr")
colnames(table1p1) <- c("ramdom_c","%","follow_c","%","ramdom_m","%","follow_m","%","ramdom_h","%","follow_h","%")
table1p1

####################################### Part 2: For contacts in 128 followed-up households ################################################

demog <- read.csv(paste(dir, "adherence_m.csv", sep=""), header=TRUE)

demog <- demog[,1:5] # hhID, member, age, sex
tab1.2 <- merge(demog,hchar,by="hhID")
contact <- tab1.2[tab1.2$member>0,]

ct_c <- contact[contact$intervention==1,]
ct_m <- contact[contact$intervention==2,]
ct_h <- contact[contact$intervention==3,]

table1p2 <- matrix(rep(NA,36),ncol=6,byrow=FALSE)

#### Age group ####

table1p2[1:4,1] <- c(length(na.exclude(ct_c$age[ct_c$age<=15])), length(na.exclude(ct_c$age[ct_c$age<=30&ct_c$age>15])),
		     length(na.exclude(ct_c$age[ct_c$age<=50&ct_c$age>30])), length(na.exclude(ct_c$age[ct_c$age>50])))

table1p2[1:4,3] <- c(length(na.exclude(ct_m$age[ct_m$age<=15])), length(na.exclude(ct_m$age[ct_m$age<=30&ct_m$age>15])),
		     length(na.exclude(ct_m$age[ct_m$age<=50&ct_m$age>30])), length(na.exclude(ct_m$age[ct_m$age>50])))

table1p2[1:4,5] <- c(length(na.exclude(ct_h$age[ct_h$age<=15])), length(na.exclude(ct_h$age[ct_h$age<=30&ct_h$age>15])),
		     length(na.exclude(ct_h$age[ct_h$age<=50&ct_h$age>30])), length(na.exclude(ct_h$age[ct_h$age>50])))

#### Sex ####

table1p2[5,1] <- length(ct_c$sex[!is.na(ct_c$sex)&ct_c$sex=="m"])
table1p2[5,3] <- length(ct_m$sex[!is.na(ct_m$sex)&ct_m$sex=="m"])
table1p2[5,5] <- length(ct_h$sex[!is.na(ct_h$sex)&ct_h$sex=="m"])

#### Vaccine history ####

table1p2[6,1] <- sum(ct_c$vaccine, na.rm = TRUE)
table1p2[6,3] <- sum(ct_m$vaccine, na.rm = TRUE)
table1p2[6,5] <- sum(ct_h$vaccine, na.rm = TRUE)

table1p2[,2] <- round(table1p2[,1]/nrow(ct_c),2)
table1p2[,4] <- round(table1p2[,3]/nrow(ct_m),2)
table1p2[,6] <- round(table1p2[,5]/nrow(ct_h),2)

rownames(table1p2) <- c("0-15yr","16-30yr","31-50yr","50+yr","men","vaccine")
colnames(table1p2) <- c("follow_c","%","follow_m","%","follow_h","%")
table1p2

# End of script
