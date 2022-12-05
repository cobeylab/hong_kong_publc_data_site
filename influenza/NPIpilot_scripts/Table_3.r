#
# Script to reproduce information in Table 3 (Factors affecting SARs) from:
#
# Cowling BJ, Fung ROP, Cheng KY, Fang VJ, Chan KH, Seto WH, et al.
# Preliminary findings of a randomized trial of non-pharmaceutical
# interventions to prevent influenza transmission in households.
# PLoS ONE, 2008; 3(5):e2101.
# http://dx.doi.org/10.1371/journal.pone.0002101
#
# Last updated by Vicky Fang and Ben Cowling
# January 5, 2009

# note - results may differ slightly from those in the PLoS ONE article for three reasons, in order of likelihood:
#  1 .. inclusion of new data (particularly new laboratory test results) subsequent to publication.
#  2 .. updated data (e.g. transcription errors detected by ongoing data cleaning processes).
#  3 .. bugs in earlier scripts.

# make sure you have downloaded the package 'gee'
require(gee)

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

hc <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""), header=TRUE)

###---------------------------------------------------- Lab-confirmed secondary cases ------------------------------------------------###

table3 <- matrix(rep(NA,169),ncol=13,byrow=FALSE)

hc <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
for (i in 1:nrow(hc)){
     if(!is.na(hc$PCR[i]) & (is.na(hc$culture[i])|(!is.na(hc$culture[i])&hc$culture[i]==0)) ) hc$culture[i] <- hc$PCR[i]
}
hc <- hc[-5]

hculture <- data.frame(hhID = baseflu$hhID, member = baseflu$member)

hc.temp <- reshape(hc, timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="culture")
hculture <- merge(hculture,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hculture) <- c("hhID","member","V0","V1","V2","V3","V4")

## exd_index: none of V0/V1 culture is A/B; d_index: both V0 and V1 culture is 0; Contact exclusion: V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
                  | (is.na(hculture$V1[i]) & hculture$V0[i]==0) | (hculture$V0[i]==0 & hculture$V1[i]==0) )) 
              {hculture$exd_index[i]<-1}    else {hculture$exd_index[i]<-0}
     if(hculture$member[i]==0 &  !is.na(hculture$V0[i]) & !is.na(hculture$V1[i]) & hculture$V0[i]==0 & hculture$V1[i]==0)
              {hculture$d_index[i]<-1}      else {hculture$d_index[i]<-0}
     if(hculture$member[i]!=0 & ( !is.na(hculture$V1[i]) & (hculture$V1[i]=="A" | hculture$V1[i]=="B") ))
              {hculture$exd_contact[i]=1}   else{hculture$exd_contact[i]=0}
}

## Define the household which should be excluded as long as index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]   # for calculating secondary cases
d_index <- hculture[hculture$member==0,c(1,9)]     # for calculating SAR

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
hculture <- merge(hculture[,-8], d_index)
dim(hculture)

for (i in 1:nrow(hculture)){
    if ( hculture$exd_index[i]==1 | hculture$exd_contact[i] ==1)
       {hculture$exclude[i] <-1}
       else  {hculture$exclude[i] <-0}
}

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$exclude[i] == 0 & !is.na(hculture$V1[i]) & ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) |
               (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) | (hculture$V4[i] !=0 & !is.na(hculture$V4[i])) ) )
              {hculture$labsedcase[i] <- 1}
     else {hculture$labsedcase[i] <- 0}
}

###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################

qv <- read.csv(paste(dir, "clinicdat_h.csv", sep=""), header=TRUE)
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""), header=TRUE)

### Add variable 'vaccine' <- 1 if the subject received vaccination in any of the past 1 year
for ( i in 1:nrow(hculture)){
    if ( (!is.na(baseflu$vaccine07[i]) & baseflu$vaccine07[i] == 1) )
          {hculture$vaccine[i] <-1}    else {hculture$vaccine[i] <- 0}
}

### Add variable 'intervention' and 'indexdelay' (time from symptom onset to intervention)
hchar <- merge(hchar,qv[,c(2,25)],by="hhID",all.x=TRUE)
hchar$onsettime[hchar$onsettime!=5] <- floor(hchar$onsettime[hchar$onsettime!=5]/2)
hchar$onsettime[hchar$onsettime==5] <- 3
hchar$v1_day <- hchar$v1_day+hchar$onsettime # from symptom onset to intervention
hculture <- merge(hculture,hchar[,c(1,2,5)],by="hhID",all.x=TRUE)
hculture$intervention <- factor(hculture$intervention, levels=1:3, labels=c("control", "mask", "hand"))
hculture$indexdelay <- 1*(hculture$v1_day<=1)

### Add variable 'indexsex' & 'indexage'
hculture <- merge(hculture,baseflu[baseflu$member==0,c(1,3,4)],by="hhID",all.x=TRUE)
names(hculture)[17:18] <- c("indexsex","indexage")
hculture$indexagegp <- cut(hculture$indexage, c(-0.1,15,100))

### Add variable 'age' and 'sex' (for household contact)
hculture$age <- baseflu$age
hculture$sex <- baseflu$sex
hculture$agegp <- cut(hculture$age,c(-0.1,15,100))

#####################################  Compute SAR based on lab-confirmed results  ########################################################

hculture.sar <- hculture[hculture$member!=0 & hculture$d_index == 0,]

inf.gee <- gee(labsedcase~intervention+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hculture.sar, corstr = "exchangeable", family="binomial")

results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results <- round(results, 2)

table3[c(1,4,6,8,10,12),c(2,5,8,11)] <- 1.00
table3[c(2,3,5,7,9,11,13),2:4] <- cbind(results$OR[-1],results$lower.CI[-1],results$upper.CI[-1])

###################################### n ##################################################################################################

new <- hculture.sar[!is.na(hculture.sar$agegp),]

table3[,1] <- c(dim(new[new$intervention=="control",])[1], dim(new[new$intervention=="mask",])[1], dim(new[new$intervention=="hand",])[1],
        dim(new[new$age<=15,])[1], dim(new[new$age>15,])[1], dim(new[new$sex=="f",])[1], dim(new[new$sex=="m",])[1],
        dim(new[new$vaccine==0,])[1], dim(new[new$vaccine==1,])[1], 
        length(unique(new$hhID[new$indexage<=15])), length(unique(new$hhID[new$indexage>15])),
        length(unique(new$hhID[new$indexsex=="f"])), length(unique(new$hhID[new$indexsex=="m"])) )


###----------------------------------------------------- Clinic secondary cases --------------------------------------------------------###

symptom <- read.csv(paste(dir, "symptomday_d.csv", sep=""), header=TRUE)
symptom$day <- rep(0:9,nrow(symptom)/10)

for (k in 1:3){

if(k==1){
    # Clinical definition 1
    symptom$fever <- 1*(as.numeric(symptom$bodytemp)>=38)
    symptom$indicate <- symptom$cough+symptom$sthroat+symptom$rnose+symptom$tired+symptom$headache+symptom$sjoint+symptom$pmuscle
    symptom$flu <- 1*( symptom$fever==1 | symptom$indicate>=2 )
}

if(k==2){
    # Clinical definition 2
    symptom$fever <- 1*(as.numeric(symptom$bodytemp)>=37.8)
    symptom$indicate <- symptom$fever+symptom$cough+symptom$headache+symptom$sthroat+symptom$pmuscle
    symptom$flu <- 1*(symptom$indicate>=2)
}

if(k==3){
    # Clinical definition 3
    symptom$fever <- 1*(as.numeric(symptom$bodytemp)>=37.8)
    symptom$indicate <- symptom$cough+symptom$sthroat
    symptom$flu <- 1*(symptom$fever==1 & symptom$indicate>=1)
}

##### Followed by either one clinical flu definition...

hclinic <- data.frame(hhID = baseflu$hhID)
hclinic$member <- 0
for(i in 2:nrow(hclinic)){
  if(hclinic$hhID[i]==hclinic$hhID[i-1]) hclinic$member[i] <- hclinic$member[i-1]+1
}

symptom.temp <- reshape(symptom[c(1:3,25)], timevar="day", idvar=c("hhID","member"), direction="wide", v.names="flu")
hclinic <- merge(hclinic,symptom.temp, by=c("hhID","member"), all.x=TRUE)
names(hclinic) <- c("hhID","member","day0","day1","day2","day3","day4","day5","day6","day7","day8","day9")

hclinic$exclude <- hculture$exclude

## Define secondary cases
for (i in 1:nrow(hclinic)){
    if ( hclinic$member[i] != 0 & hclinic$exclude[i] == 0 & !is.na(hclinic$day0[i]) &
         ( (!is.na(hclinic$day1[i]) & hclinic$day1[i]==1) | (!is.na(hclinic$day2[i]) & hclinic$day2[i]==1) | 
       (!is.na(hclinic$day3[i]) & hclinic$day3[i]==1) | (!is.na(hclinic$day4[i]) & hclinic$day4[i]==1) | 
       (!is.na(hclinic$day5[i]) & hclinic$day5[i]==1) | (!is.na(hclinic$day6[i]) & hclinic$day6[i]==1) | 
       (!is.na(hclinic$day7[i]) & hclinic$day7[i]==1) | (!is.na(hclinic$day8[i]) & hclinic$day8[i]==1) | 
       (!is.na(hclinic$day9[i]) & hclinic$day9[i]==1) ))
              {hclinic$clinicsedcase[i] <- 1}
     else {hclinic$clinicsedcase[i] <- 0}
}

###################### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp ##################################

hclinic$d_index <- hculture$d_index
hclinic$intervention <- hculture$intervention
hclinic$indexagegp <- hculture$indexagegp
hclinic$indexsex <- hculture$indexsex
hclinic$agegp <- hculture$agegp
hclinic$sex <- hculture$sex
hclinic$vaccine <- hculture$vaccine

#####################################  Compute SAR based on lab-confirmed results  ######################################################

library(gee)
hclinic.sar <- hclinic[hclinic$member!=0 & hclinic$d_index == 0,]

inf.gee <- gee(clinicsedcase~intervention+agegp+sex+vaccine+indexagegp+indexsex, id=factor(hhID),
  data=hclinic.sar, corstr = "exchangeable", family="binomial")

results <- data.frame(beta=inf.gee$coef, se=sqrt(diag(inf.gee[[20]])),
  row.names=names(inf.gee$coef))
results$OR <- exp(results$beta)
results$lower.CI <- exp(results$beta - 1.96*results$se)
results$upper.CI <- exp(results$beta + 1.96*results$se)
results <- round(results, 2) 

table3[c(2,3,5,7,9,11,13),2:4+3*k] <- cbind(results$OR[-1],results$lower.CI[-1],results$upper.CI[-1])
}

rownames(table3) <- c("control","mask","hand","child","adult","female","male","not vac","vaccine",
              "child index","adult index","female index","male index")
colnames(table3) <- c("n","lab-OR","CI-low","CI-up","cd1-OR","CI-low","CI-up","cd2-OR","CI-low","CI-up","cd3-OR","CI-low","CI-up")
table3

# End of script
