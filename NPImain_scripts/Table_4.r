#
# Script to reproduce information in Table 4 from:
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

###-----------------------------------------------Lab-confirmed secondary cases --------------------------------------------------------###

table4 <- matrix(rep(NA,170),ncol=10,byrow=FALSE)

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i]!= 0 & hculture$analyzed[i]==1 & hculture$exd_contact[i]==0 &
             ( (hculture$V2[i]!=0 & !is.na(hculture$V2[i])) | (hculture$V3[i]!=0 & !is.na(hculture$V3[i]))  ) )
              {hculture$labsedcase[i] <- 1}
	 else {hculture$labsedcase[i] <- 0}
}

###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################

baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""))
av <- read.csv(paste(dir, "antiviral_m.csv", sep=""))

### Add variable 'vaccine' <- 1 if the subject received vaccination in past one year
hculture$vaccine <- baseflu$vaccine08

### Add variable 'indexsex' & 'indexage'
indexage <- baseflu[baseflu$member==0,c(1,4)]
names(indexage)[2] <- "indexage"
indexsex <- baseflu[baseflu$member==0,c(1,3)]
names(indexsex)[2] <- "indexsex"
hculture <- merge(hculture,indexage)
hculture <- merge(hculture,indexsex)
hculture$indexagegp <- 1*(hculture$indexage<=15)+1*(hculture$indexage<=5)

### Add variable 'intervention' and 'indexdelay' (time from symptom onset to intervention)
hculture <- merge(hculture,housechar[,c(1,2,9)])
hculture$arm <- factor(hculture$intervention, levels=c(1,3,4), labels=c("control", "hand", "handmask"))
hculture$delay[hculture$onset_v1_delay<=36] <- "d1"
hculture$delay[hculture$onset_v1_delay>36] <- "d2"

### Add variable 'age' and 'sex' (for household contact)
hculture$age <- baseflu$age
hculture$sex <- baseflu$male

# set age (10,30) according to the relationship to index - in order to get agegp
hculture$age[hculture$hhID==148&hculture$member==3] <- 10
hculture$age[hculture$hhID==317&hculture$member==4] <- 10
hculture$age[is.na(hculture$age)] <- 30
hculture$agegp <- 1*(hculture$age<=15)+1*(hculture$age<=5)

# add index antiviral use
av <- av[av$member==0,]
av$indexav <- 1
hculture <- merge(hculture,av[c(1,6)],by="hhID",all.x=TRUE)
hculture$indexav[is.na(hculture$indexav)] <- 0
hculture <- hculture[order(hculture$hhID,hculture$member),]

#####################################  Compute SAR based on lab-confirmed results  ########################################################

library(gee)
hculture.sar <- hculture[hculture$member!=0 & hculture$analyzed == 1,]

inf.gee <- gee(labsedcase~arm+factor(agegp)+sex+vaccine+factor(indexagegp)+indexsex+indexav, id=factor(hhID),
  data=hculture.sar, corstr = "exchangeable", family="binomial")

results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
results <- round(results, 2)                 # Odds ratio and 95% CI

table4[c(1,4,7,9,11,14,16),c(2,5,8)] <- 1.00
table4[c(2,3,5,6,8,10,12,13,15,17),2:4] <- cbind(results$OR[-1],results$lower.CI[-1],results$upper.CI[-1])

#### n 

new <- hculture.sar[!is.na(hculture.sar$agegp),]

table4[,1] <- c(dim(new[new$arm=="control",])[1],dim(new[new$arm=="hand",])[1],dim(new[new$arm=="handmask",])[1],
  dim(new[new$age>15,])[1], dim(new[new$age<=15&new$age>5,])[1], dim(new[new$age<=5,])[1],
  dim(new[new$sex==0,])[1], dim(new[new$sex==1,])[1],
  dim(new[new$vaccine==0,])[1],dim(new[new$vaccine==1,])[1],
  length(unique(new$hhID[new$indexage>15])),length(unique(new$hhID[new$indexage<=15&new$indexage>5])),length(unique(new$hhID[new$indexage<=5])),
  length(unique(new$hhID[new$indexsex==0])),length(unique(new$hhID[new$indexsex==1])),
  length(unique(new$hhID[new$indexav==0])),length(unique(new$hhID[new$indexav==1])))

###----------------------------------------------------- Clinic secondary cases --------------------------------------------------------###

symptom <- read.csv(paste(dir, "symptomday_d.csv", sep=""))

for (k in 1:2){

if(k==1){
  # Case 1 (Monto 2002)
  symptom$fever <- 1*(!is.na(symptom$bodytemp)&symptom$bodytemp>=37.8)
  symptom$indicate <- symptom$fever+symptom$cough+symptom$headache+symptom$sthroat+symptom$pmuscle
  symptom$flu <- 1*(symptom$indicate>=2)
}

if(k==2){
  # Case 2 (WHO)
  symptom$fever <- 1*(!is.na(symptom$bodytemp)&symptom$bodytemp>=37.8)
  symptom$indicate <- symptom$cough+symptom$sthroat
  symptom$flu <- 1*(symptom$fever==1 & symptom$indicate>=1)
}

###########################################################################################################################################

hclinic <- data.frame(hhID = hculture$hhID)
hclinic$member <- hculture$member
symptom.temp <- reshape(symptom[c(1:3,13)], timevar="day", idvar=c("hhID","member"), direction="wide", v.names="flu")
hclinic <- merge(hclinic,symptom.temp, by=c("hhID","member"), all.x=TRUE)
hclinic <- hclinic[order(hclinic$hhID,hclinic$member),]

names(hclinic) <- c("hhID","member","day0","day1","day2","day3","day4","day5","day6","day7","day8","day9")

hclinic$exd_contact <- hculture$exd_contact
hclinic$analyzed <- hculture$analyzed

## Define secondary cases
for (i in 1:nrow(hclinic)){
    if ( hclinic$member[i] != 0 & hclinic$analyzed[i] == 1 & hclinic$exd_contact[i]==0 &
         ( (!is.na(hclinic$day1[i]) & hclinic$day1[i]==1) | (!is.na(hclinic$day2[i]) & hclinic$day2[i]==1) |
	   (!is.na(hclinic$day3[i]) & hclinic$day3[i]==1) | (!is.na(hclinic$day4[i]) & hclinic$day4[i]==1) |
	   (!is.na(hclinic$day5[i]) & hclinic$day5[i]==1) | (!is.na(hclinic$day6[i]) & hclinic$day6[i]==1) |
	   (!is.na(hclinic$day7[i]) & hclinic$day7[i]==1) | (!is.na(hclinic$day8[i]) & hclinic$day8[i]==1) |
	   (!is.na(hclinic$day9[i]) & hclinic$day9[i]==1) ))
              {hclinic$clinicsedcase[i] <- 1}
	 else {hclinic$clinicsedcase[i] <- 0}
}

###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################

hclinic$arm <- hculture$arm
hclinic$indexagegp <- hculture$indexagegp
hclinic$indexsex <- hculture$indexsex
hclinic$agegp <- hculture$agegp
hclinic$sex <- hculture$sex
hclinic$vaccine <- hculture$vaccine
hclinic$indexav <- hculture$indexav
hclinic$delay <- hculture$delay

#####################################  Compute SAR based on clinic-confirmed results  ######################################################

library(gee)
hclinic.sar <- hclinic[hclinic$member!=0 & hclinic$analyzed == 1,]

inf.gee <- gee(clinicsedcase~arm+factor(agegp)+sex+vaccine+factor(indexagegp)+indexsex+indexav, id=factor(hhID),
  data=hclinic.sar, corstr = "exchangeable", family="binomial")

results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results$p.value <- 2*pnorm(-1*abs(results$beta/results$se))
results <- round(results, 2)

table4[c(2,3,5,6,8,10,12,13,15,17),2:4+3*k] <- cbind(results$OR[-1],results$lower.CI[-1],results$upper.CI[-1])
}

rownames(table4) <- c("control","mask","hand","adult","child6-15","child<=5","female","male","not vac","vaccine",
              "adult index","child6-15 index","child<=5 index","female index","male index","no antiviral","antiviral")
colnames(table4) <- c("n","lab-OR","CI-low","CI-up","cd1-OR","CI-low","CI-up","cd2-OR","CI-low","CI-up")
table4

# End of script

