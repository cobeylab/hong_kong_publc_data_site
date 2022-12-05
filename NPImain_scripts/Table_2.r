#
# Script to reproduce information in Table 2 from:
#
# Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, et al.
# Facemasks and hand hygiene to prevent influenza transmission 
# in households, a randomized trial.
# Annals of Internal Medicine, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Sep 08, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV3/"
source("http://www.hku.hk/bcowling/influenza/NPImain_scripts/Analyzed_hh.r")

q1data <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
intervention <- read.csv(paste(dir, "randomarm_407.csv", sep=""))
av <- read.csv(paste(dir, "antiviral_m.csv", sep=""))

table2p1 <- matrix(rep(NA,336),ncol=12,byrow=FALSE)

####################################### Part 1-1: For randomized index subjects (n=407) #####################################################

q1 <- q1data
intervention$hhID <- as.numeric(substr(intervention$hhID,5,7))

tab2 <- merge(intervention[2:3],q1,by="scrID",all.x=TRUE)
tab2 <- merge(tab2,av[av$member==0,c(1,3)],by="hhID",all.x=TRUE)
tab2$av <- as.character(tab2$av)
tab2$av[is.na(tab2$av)] <- 0

tab2_c <- tab2[tab2$intervention==1,]      # control
tab2_h <- tab2[tab2$intervention==3,]      # hand hygiene
tab2_m <- tab2[tab2$intervention==4,]      # hand hygiene + mask

c(dim(tab2_c)[1],dim(tab2_h)[1],dim(tab2_m)[1]) # Randomly assigned

####### Age group #######
table2p1[1:5,1] <- c(dim(tab2_c[tab2_c$age<=5,])[1], dim(tab2_c[tab2_c$age<=15&tab2_c$age>5,])[1],dim(tab2_c[tab2_c$age<=30&tab2_c$age>15,])[1],
	dim(tab2_c[tab2_c$age<=50&tab2_c$age>30,])[1],dim(tab2_c[tab2_c$age>50,])[1])

table2p1[1:5,5] <- c(dim(tab2_h[tab2_h$age<=5,])[1], dim(tab2_h[tab2_h$age<=15&tab2_h$age>5,])[1],dim(tab2_h[tab2_h$age<=30&tab2_h$age>15,])[1],
	dim(tab2_h[tab2_h$age<=50&tab2_h$age>30,])[1],dim(tab2_h[tab2_h$age>50,])[1])

table2p1[1:5,9] <- c(dim(tab2_m[tab2_m$age<=5,])[1], dim(tab2_m[tab2_m$age<=15&tab2_m$age>5,])[1],dim(tab2_m[tab2_m$age<=30&tab2_m$age>15,])[1],
	dim(tab2_m[tab2_m$age<=50&tab2_m$age>30,])[1],dim(tab2_m[tab2_m$age>50,])[1])
	
# Median age & IQR
table2p1[6,c(1:2,5:6,9:10)] <- c( round(quantile(tab2_c$age,0.5)), paste(round(quantile(tab2_c$age,0.25)),"-",round(quantile(tab2_c$age,0.75)),sep=""),
                                  round(quantile(tab2_h$age,0.5)), paste(round(quantile(tab2_h$age,0.25)),"-",round(quantile(tab2_h$age,0.75)),sep=""),
                                  round(quantile(tab2_m$age,0.5)), paste(round(quantile(tab2_m$age,0.25)),"-",round(quantile(tab2_m$age,0.75)),sep=""))

####### Sex #######
table2p1[7,c(1,5,9)] <- c(dim(tab2_c[tab2_c$male==1,])[1], dim(tab2_h[tab2_h$male==1,])[1], dim(tab2_m[tab2_m$male==1,])[1])

####### Symptoms #######
table2p1[8:14,1] <- colSums(tab2_c[,c(6:12)],na.rm=TRUE)
table2p1[8:14,5] <- colSums(tab2_h[,c(6:12)],na.rm=TRUE)
table2p1[8:14,9] <- colSums(tab2_m[,c(6:12)],na.rm=TRUE)

####### Delay from symptom onset to randomization #######
table2p1[15:19,1] <- table(tab2_c$onsettime)
table2p1[15:19,5] <- table(tab2_h$onsettime)
table2p1[15:19,9] <- table(tab2_m$onsettime)

table2p1[c(1:5,7:19),2] <- round(as.numeric(table2p1[c(1:5,7:19),1])/nrow(tab2_c),2)
table2p1[c(1:5,7:19),6] <- round(as.numeric(table2p1[c(1:5,7:19),5])/nrow(tab2_h),2)
table2p1[c(1:5,7:19),10] <- round(as.numeric(table2p1[c(1:5,7:19),9])/nrow(tab2_m),2)


####################################### Part 1-2: For 259 analyzed index subjects #######################################

hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))

hchar$analyzed <- hculture$analyzed[hculture$member==0]
hchar <- merge(hchar,q1data[c(2,12)],by="hhID",all.x=TRUE)
hchar$onsettime <- hchar$onsettime*12
hchar$delay <- hchar$onset_v1_delay-hchar$onsettime
hchar <- hchar[hchar$analyzed==1,c(1:3,12)] # hhID, intervention, familysize, delay
tab2.2 <- merge(tab2[-3],hchar,by="hhID")
tab2.2$fever <- 1*(tab2.2$bodytemp>=37.8)

tab2_c <- tab2.2[tab2.2$intervention==1,]
tab2_h <- tab2.2[tab2.2$intervention==3,]
tab2_m <- tab2.2[tab2.2$intervention==4,]

c(dim(tab2_c)[1],dim(tab2_h)[1],dim(tab2_m)[1])   # Analyzed

####### Age group #######
table2p1[1:5,3] <- c(dim(tab2_c[tab2_c$age<=5,])[1], dim(tab2_c[tab2_c$age<=15&tab2_c$age>5,])[1],dim(tab2_c[tab2_c$age<=30&tab2_c$age>15,])[1],
	dim(tab2_c[tab2_c$age<=50&tab2_c$age>30,])[1],dim(tab2_c[tab2_c$age>50,])[1])

table2p1[1:5,7] <- c(dim(tab2_h[tab2_h$age<=5,])[1], dim(tab2_h[tab2_h$age<=15&tab2_h$age>5,])[1],dim(tab2_h[tab2_h$age<=30&tab2_h$age>15,])[1],
	dim(tab2_h[tab2_h$age<=50&tab2_h$age>30,])[1],dim(tab2_h[tab2_h$age>50,])[1])

table2p1[1:5,11] <- c(dim(tab2_m[tab2_m$age<=5,])[1], dim(tab2_m[tab2_m$age<=15&tab2_m$age>5,])[1],dim(tab2_m[tab2_m$age<=30&tab2_m$age>15,])[1],
	dim(tab2_m[tab2_m$age<=50&tab2_m$age>30,])[1],dim(tab2_m[tab2_m$age>50,])[1])
	
# Median age & IQR
table2p1[6,c(3:4,7:8,11:12)] <- c( round(quantile(tab2_c$age,0.5)), paste(round(quantile(tab2_c$age,0.25)),"-",round(quantile(tab2_c$age,0.75)),sep=""),
                                  round(quantile(tab2_h$age,0.5)), paste(round(quantile(tab2_h$age,0.25)),"-",round(quantile(tab2_h$age,0.75)),sep=""),
                                  round(quantile(tab2_m$age,0.5)), paste(round(quantile(tab2_m$age,0.25)),"-",round(quantile(tab2_m$age,0.75)),sep=""))

####### Sex #######
table2p1[7,c(3,7,11)] <- c(dim(tab2_c[tab2_c$male==1,])[1], dim(tab2_h[tab2_h$male==1,])[1], dim(tab2_m[tab2_m$male==1,])[1])

####### Symptoms #######
table2p1[8:14,3] <- colSums(tab2_c[,c(20,6:11)],na.rm=TRUE)
table2p1[8:14,7] <- colSums(tab2_h[,c(20,6:11)],na.rm=TRUE)
table2p1[8:14,11] <- colSums(tab2_m[,c(20,6:11)],na.rm=TRUE)

####### Delay from symptom onset to randomization #######
table2p1[15:19,3] <- table(tab2_c$onsettime,exclude=NULL)
table2p1[15:19,7] <- table(tab2_h$onsettime)
table2p1[15:19,11] <- table(tab2_m$onsettime)

####### delay from randomization to intervention #######
table2p1[20:23,3] <- c(dim(tab2_c[tab2_c$delay<=12,])[1],dim(tab2_c[tab2_c$delay>12&tab2_c$delay<=24,])[1],dim(tab2_c[tab2_c$delay>24&tab2_c$delay<=36,])[1],dim(tab2_c[tab2_c$delay>36,])[1])
table2p1[20:23,7] <- c(dim(tab2_h[tab2_h$delay<=12,])[1],dim(tab2_h[tab2_h$delay>12&tab2_h$delay<=24,])[1],dim(tab2_h[tab2_h$delay>24&tab2_h$delay<=36,])[1],dim(tab2_h[tab2_h$delay>36,])[1])
table2p1[20:23,11] <- c(dim(tab2_m[tab2_m$delay<=12,])[1],dim(tab2_m[tab2_m$delay>12&tab2_m$delay<=24,])[1],dim(tab2_m[tab2_m$delay>24&tab2_m$delay<=36,])[1],dim(tab2_m[tab2_m$delay>36,])[1])

# prescribed antiviral #
table2p1[24:27,3] <- c(dim(tab2_c[tab2_c$av=="tamiflu",])[1],dim(tab2_c[tab2_c$av=="amantadine",])[1],dim(tab2_c[tab2_c$av=="relenza",])[1],dim(tab2_c[tab2_c$av=="ribavirin",])[1])
table2p1[24:27,7] <- c(dim(tab2_h[tab2_h$av=="tamiflu",])[1],dim(tab2_h[tab2_h$av=="amantadine",])[1],dim(tab2_h[tab2_h$av=="relenza",])[1],dim(tab2_h[tab2_h$av=="ribavirin",])[1])
table2p1[24:27,11] <- c(dim(tab2_m[tab2_m$av=="tamiflu",])[1],dim(tab2_m[tab2_m$av=="amantadine",])[1],dim(tab2_m[tab2_m$av=="relenza",])[1],dim(tab2_m[tab2_m$av=="ribavirin",])[1])

# Median household size
table2p1[28,c(3,7,11)] <- c(median(tab2_c$familysize), median(tab2_h$familysize), median(tab2_m$familysize))

table2p1[c(1:5,7:27),4] <- round(as.numeric(table2p1[c(1:5,7:27),3])/nrow(tab2_c),2)
table2p1[c(1:5,7:27),8] <- round(as.numeric(table2p1[c(1:5,7:27),7])/nrow(tab2_h),2)
table2p1[c(1:5,7:27),12] <- round(as.numeric(table2p1[c(1:5,7:27),11])/nrow(tab2_m),2)

rownames(table2p1) <- c("<=5yr","6-15yr","16-30yr","31-50yr","50+yr","median age (IQR)","men","fever","headache","sthroat","cough","myalgia","rnose","phlegm",
                        "0-12h","12-24h","24-36h","36-48h","48-60h","0-12h","12-24h","24-36h","36-48h","oseltamivir","amantadine","zanamivir","ribavirin","median household size")
colnames(table2p1) <- c("ramdom_c","%","follow_c","%","ramdom_h","%","follow_h","%","ramdom_m","%","follow_m","%")
table2p1


####################################### Part 2: For contacts in 259 analyzed households ################################################

baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""))

demog <- baseflu[1:4] # hhID, member, age, sex
demog$vaccine <- baseflu$vaccine08
tab2.3 <- merge(demog,hchar,by="hhID")
contact <- tab2.3[tab2.3$member>0,]

ct_c <- contact[contact$intervention==1,]
ct_h <- contact[contact$intervention==3,]
ct_m <- contact[contact$intervention==4,]

c(dim(ct_c)[1],dim(ct_h)[1],dim(ct_m)[1])   # contact in analyzed households

table2p2 <- matrix(rep(NA,54),ncol=6,byrow=FALSE)

#### Age group ####
table2p2[1:6,1] <- c(length(na.exclude(ct_c$age[ct_c$age<=5])),length(na.exclude(ct_c$age[ct_c$age<=15&ct_c$age>5])),length(na.exclude(ct_c$age[ct_c$age<=30&ct_c$age>15])),
	                   length(na.exclude(ct_c$age[ct_c$age<=50&ct_c$age>30])),length(na.exclude(ct_c$age[ct_c$age>50])),length(ct_c$age[is.na(ct_c$age)]))
table2p2[1:6,3] <- c(length(na.exclude(ct_h$age[ct_h$age<=5])),length(na.exclude(ct_h$age[ct_h$age<=15&ct_h$age>5])),length(na.exclude(ct_h$age[ct_h$age<=30&ct_h$age>15])),
	                   length(na.exclude(ct_h$age[ct_h$age<=50&ct_h$age>30])),length(na.exclude(ct_h$age[ct_h$age>50])),length(ct_h$age[is.na(ct_h$age)]))
table2p2[1:6,5] <- c(length(na.exclude(ct_m$age[ct_m$age<=5])),length(na.exclude(ct_m$age[ct_m$age<=15&ct_m$age>5])),length(na.exclude(ct_m$age[ct_m$age<=30&ct_m$age>15])),
	                   length(na.exclude(ct_m$age[ct_m$age<=50&ct_m$age>30])),length(na.exclude(ct_m$age[ct_m$age>50])),length(ct_m$age[is.na(ct_m$age)]))

# Median age & IQR
table2p2[7,] <- c( round(quantile(ct_c$age,0.5,na.rm=TRUE)), paste(round(quantile(ct_c$age,0.25,na.rm=TRUE)),"-",round(quantile(ct_c$age,0.75,na.rm=TRUE)),sep=""),
                   round(quantile(ct_h$age,0.5,na.rm=TRUE)), paste(round(quantile(ct_h$age,0.25,na.rm=TRUE)),"-",round(quantile(ct_h$age,0.75,na.rm=TRUE)),sep=""),
                   round(quantile(ct_m$age,0.5,na.rm=TRUE)), paste(round(quantile(ct_m$age,0.25,na.rm=TRUE)),"-",round(quantile(ct_m$age,0.75,na.rm=TRUE)),sep=""))
                   
#### Sex ####
table2p2[8,c(1,3,5)] <- c(sum(ct_c$male, na.rm = TRUE),sum(ct_h$male, na.rm = TRUE),sum(ct_m$male, na.rm = TRUE))

#### Vaccine history ####
table2p2[9,c(1,3,5)] <- c(sum(ct_c$vaccine, na.rm = TRUE),sum(ct_h$vaccine, na.rm = TRUE),sum(ct_m$vaccine, na.rm = TRUE))

table2p2[c(1:6,8:9),2] <- round(as.numeric(table2p2[c(1:6,8:9),1])/nrow(ct_c),2)
table2p2[c(1:6,8:9),4] <- round(as.numeric(table2p2[c(1:6,8:9),3])/nrow(ct_h),2)
table2p2[c(1:6,8:9),6] <- round(as.numeric(table2p2[c(1:6,8:9),5])/nrow(ct_m),2)

rownames(table2p2) <- c("<=5yr","6-15yr","16-30yr","31-50yr","50+yr","unknown","median age (IQR)","men","vaccine")
colnames(table2p2) <- c("follow_c","%","follow_m","%","follow_h","%")
table2p2

# End of script

