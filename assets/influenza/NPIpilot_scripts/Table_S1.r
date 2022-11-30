#
# Script to reproduce information in Table S1 (clinical flu definition vs. gold standard) from:
#
# Cowling BJ, Fung ROP, Cheng KY, Fang VJ, Chan KH, Seto WH, et al.
# Preliminary findings of a randomized trial of non-pharmaceutical
# interventions to prevent influenza transmission in households.
# PLoS ONE, 2008; 3(5):e2101.
# http://dx.doi.org/10.1371/journal.pone.0002101
#
# Last updated by Vicky Fang and Ben Cowling
# January 5, 2009

# make sure you have downloaded the package 'ROCR'
require(ROCR)

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

hc <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""), header=TRUE)
tableS1 <- matrix(rep(NA,15),ncol=5,byrow=FALSE)

####################################################### Lab-confirmed secondary cases ####################################################

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


####################################################### Clinic secondary cases ###########################################################

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

############################################ Calculate sensitivity / specificity / Area under ROC ########################################

x <- hculture[,c(1,2,11,12)] # hhID,member, exclude, labsedcase 
y <- hclinic[,c(1,2,14)]
z <- merge(x,y)
zz <- z[z$member!=0&z$exclude==0,]

clinic <- zz$clinicsedcase
lab <- zz$labsedcase          # gold standard
p1 <- prediction(clinic, lab)
ss <- performance(p1, "sens", "spec")
p2 <- performance(p1, "auc")

################# Bootstrap for area under ROC ########################

boot.auc <- function(reps,fludata){
    n <- length(fludata[,1])
    temp <- rep(NA,reps)
    output <- data.frame(auc=temp)
    for(i in 1:reps){
        a <- sample(1:n,replace=TRUE)
    temp.fludata <- fludata[a,]
        temp.clinic <- temp.fludata$clinicsedcase
    temp.lab <- temp.fludata$labsedcase
    temp.p1 <- prediction(temp.clinic,temp.lab)
    temp.p2 <- performance(temp.p1,"auc")
    output$auc[i] <- temp.p2@y.values[[1]]
    }
    output
}

flu.boot <- boot.auc(1000,zz)
# sensitivity, specificity, area under roc, and 95% CI
tableS1[k,] <- round(c(ss@y.values[[1]][2], ss@x.values[[1]][2], p2@y.values[[1]], quantile(flu.boot$auc,c(0.025,0.975))),2)
}

rownames(tableS1) <- c("cd1","cd2","cd3")
colnames(tableS1) <- c("sensitivity","specificity","AUR","CI-low","CI-up")
tableS1


# End of script
