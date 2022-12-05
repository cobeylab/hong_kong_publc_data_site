#
# Script to reproduce the data used in:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

#
# Construct the data frame for analysis (n=384 index case in total) (n=331 analyzed hhs)
#

dir.pilot <- "http://sph.hku.hk/data/HongKongNPIpilotV4/"
dir.main <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"

clinic.pilot <- read.csv(paste(dir.pilot,"clinicdat_h.csv", sep=""))
clinic.main <- read.csv(paste(dir.main,"clinicdat_h.csv", sep=""))
clinic.pilot$hhID <- as.numeric(substr(clinic.pilot$hhID,3,6))
clinic.main$hhID <- clinic.main$hhID+8000
clinic.pilot$male <- 1*(clinic.pilot$sex=="m")
clinic <- rbind(clinic.pilot[c("hhID","male","age","self_fever","headache","sthroat","cough","pmuscle","rnose","phlegm","onsettime",
                               "QVres","bodytemp","antiviral","antibiotics","antihistamine","antipyretic","steroid")],
                clinic.main[c("hhID","male","age","self_fever","headache","sthroat","cough","pmuscle","rnose","phlegm","onsettime",
                              "QVres","bodytemp","antiviral","antibiotics","antihistamine","antipyretic","steroid")])

hc.pilot <- read.csv(paste(dir.pilot, "home_culture_qPCR.csv", sep=""),as.is=T)
hc.main <- read.csv(paste(dir.main, "home_pcr.csv", sep=""))
hc.pilot$hhID <-  as.numeric(substr(hc.pilot$hhID,3,6))
hc.main$hhID <- hc.main$hhID+8000
hc.pilot$culture <- 1*(hc.pilot$culture>0)
hc.pilot$PCR <- 1*(hc.pilot$qPCR>0)
hc.pilot$PCR[is.na(hc.pilot$PCR)] <- hc.pilot$culture[is.na(hc.pilot$PCR)]
hc.pilot$PCR[!is.na(hc.pilot$culture)&hc.pilot$culture>0] <- hc.pilot$culture[!is.na(hc.pilot$culture)&hc.pilot$culture>0]
hc.main$PCR <- 1*(hc.main$qPCR>0)
hc <- rbind(hc.pilot[c("hhID","member","visit","PCR")],hc.main[c("hhID","member","visit","PCR")])

demog.pilot <- read.csv(paste(dir.pilot, "adherence_m.csv", sep=""))
demog.main <- read.csv(paste(dir.main, "adherence_m.csv", sep=""))
demog.pilot$hhID <-  as.numeric(substr(demog.pilot$hhID,3,6))
demog.main$hhID <- demog.main$hhID+8000
demog.pilot$male <- 1*(demog.pilot$sex=="m")
names(demog.pilot)[5:6] <- c("vaccine","mask")
names(demog.main)[5] <- c("vaccine")
demog.main$mask_usage <- demog.main$given_mask-demog.main$remain_mask
demog.main$smallgel_usage <- demog.main$smallgel_given-demog.main$smallgel_remain
demog <- rbind(demog.pilot[c("hhID","member","male","age","vaccine","mask","mask_usage","soap","handrub","washhand","smallgel_usage","chronic_disease")],
               demog.main[c("hhID","member","male","age","vaccine","mask","mask_usage","soap","handrub","washhand","smallgel_usage","chronic_disease")])

av.pilot <- read.csv(paste(dir.pilot, "antiviral_m.csv", sep=""))
av.main <- read.csv(paste(dir.main, "antiviral_m.csv", sep=""))
av.pilot$hhID <- as.numeric(substr(av.pilot$hhID,3,6))
av.main$hhID <- av.main$hhID+8000
av <- rbind(av.pilot,av.main)

hchar.pilot <- read.csv(paste(dir.pilot, "hchar_h.csv", sep=""))
hchar.main <- read.csv(paste(dir.main, "hchar_h.csv", sep=""))
hchar.pilot$hhID <- as.numeric(substr(hchar.pilot$hhID,3,6))
hchar.main$hhID <- hchar.main$hhID+8000
hchar.main$v4_day <- NA
hchar <- rbind(hchar.pilot[c("hhID","intervention","familysize","house_size","clinic_date","clinic_day","v1_day","v2_day","v3_day","v4_day")],
               hchar.main[c("hhID","intervention","familysize","house_size","clinic_date","clinic_day","v1_day","v2_day","v3_day","v4_day")])

symptom.pilot <- read.csv(paste(dir.pilot, "symptomday_d.csv", sep=""))
symptom.main <- read.csv(paste(dir.main, "symptomday_d.csv", sep=""))
symptom.pilot$hhID <- as.numeric(substr(symptom.pilot$hhID,3,6))
symptom.main$hhID <- symptom.main$hhID+8000
symptom <- rbind(symptom.pilot[c("hhID","member","day","bodytemp","headache","sthroat","cough","pmuscle","rnose","phlegm")],
               symptom.main[c("hhID","member","day","bodytemp","headache","sthroat","cough","pmuscle","rnose","phlegm")])

### Index Symptom Scores at clinic
clinic$fever <- 1*(clinic$bodytemp>=37.8)
clinic$fever[is.na(clinic$fever)] <- clinic$self_fever[is.na(clinic$fever)]
clinic$b_m0_score<-clinic$fever+clinic$headache+clinic$sthroat+clinic$cough+clinic$pmuscle+clinic$rnose+clinic$phlegm
clinic$b_m0_rscore <-clinic$sthroat+clinic$cough+clinic$rnose+clinic$phlegm

#only include those with Homevisit
clinic2 <- clinic[!is.na(clinic$hhID),]

# PCR results
hc1 <- hc[hc$visit>0,]
hc2 <- reshape(hc1,timevar="visit",idvar=c("hhID","member"),direction="wide")
hc3 <- merge(demog[1:2],hc2[1:5],by=c("hhID","member"),all.x=TRUE)

index <-hc3[hc3$member==0,]
index$als.index <- 1*(index$PCR.1==1&!is.na(index$PCR.1))
clinic2 <- merge(clinic2,index[c("hhID","als.index")],all.x=TRUE)
clinic3 <- clinic2[clinic2$als.index==1,1:(ncol(clinic2)-1)]

# add antiviral usage
clinic4 <- merge(clinic3,av[av$member==0,-2],all.x=TRUE)
clinic4$av <- as.character(clinic4$av)
clinic4$av[is.na(clinic4$av)] <- "none"

# add interval from onset to tamiflu usage
clinic5 <- clinic4[clinic4$av=="none"|clinic4$av=="tamiflu",]  # exclude index who took antiviral but not tamiflu
clinic5$tami_onset <- NA
for (i in 1:nrow(clinic5)) {
    if (clinic5$av[i]=="none") clinic5$tami_onset[i] <- "1none"
    else if (clinic5$onsettime[i]==5) clinic5$tami_onset[i] <- "2tamiflu>48hrs"
    else if (clinic5$onsettime[i]==4|clinic5$onsettime[i]==3) clinic5$tami_onset[i] <- "3tamiflu25-48hrs"
    else clinic5$tami_onset[i] <- "4tamiflu<24hrs"
}

# add PCR results
clinic6 <- merge(clinic5,hc2[hc2$member==0,-2],all.x=TRUE)

# add age group
clinic6$agegp <- cut(clinic6$age,c(0,5,12,17,100))
clinic6$agegp <- factor(clinic6$agegp,labels=c("2youngerchildren0-5yo","3olderchildren6-12yo","4adolescents13-17yo","1adults18+yo"))

# add vaccine status
clinic7 <- merge(clinic6,demog[demog$member==0,c(1,5,12)],all.x=TRUE)
clinic7$vaccine[is.na(clinic7$vaccine)] <- 0
clinic7$chronic_disease[is.na(clinic7$chronic_disease)] <- 0

clinic7$flu.type <- "A"
clinic7$flu.type[clinic7$QVres==2] <- "B"
clinic7$flu.type[clinic7$hhID==7159] <- "A"  # culture result A
index.tami <- clinic7

# remove hhs with co-index cases
co.index <- unique(hc3$hhID[hc3$PCR.1==1&!is.na(hc3$PCR.1)&hc3$member>0])
clinic8 <- clinic7[!(clinic7$hhID%in%co.index),]

index.hh <- clinic8

#
# alleviation of all symptoms ---------------------------------------------------------------------------------------------
#

symptom.index <-symptom[symptom$member==0,]
symptom.index$fever <- 1*(symptom.index$bodytemp>=37.8)

#Replacing NA symptoms scores by 0 if not all individual symptom score on that day are NA
for (i in 1:nrow(symptom.index)){
    for (j in 5:11){
      if ((!is.na(symptom.index$headache[i]) | !is.na(symptom.index$sthroat[i]) | !is.na(symptom.index$fever[i])
          | !is.na(symptom.index$cough[i]) | !is.na(symptom.index$pmuscle[i]) | !is.na(symptom.index$rnose[i])
          | !is.na(symptom.index$phlegm[i])) & is.na(symptom.index[i,j])) symptom.index[i,j] <- 0
    }
}

#summing up to total symptom score (/7)
symptom.index$hv_score <- symptom.index$fever+symptom.index$headache+symptom.index$sthroat+symptom.index$cough+symptom.index$pmuscle+symptom.index$rnose+symptom.index$phlegm
symptom.index$hv_rscore <-symptom.index$sthroat+symptom.index$cough+symptom.index$rnose+symptom.index$phlegm

symptom.index2 <-symptom.index[,c("hhID","day","hv_score","fever","hv_rscore")]
symptom.index3 <- reshape(symptom.index2,idvar="hhID",v.names=c("hv_score","fever","hv_rscore"),timevar="day",direction="wide")
symptom.index4 <- symptom.index3[symptom.index3$hhID %in% clinic7$hhID,]
m0_hv <- symptom.index4[c("hhID","hv_score.0","hv_score.1","hv_score.2","hv_score.3","hv_score.4","hv_score.5","hv_score.6","hv_score.7","hv_score.8","hv_score.9",
                        "fever.0","fever.1","fever.2","fever.3","fever.4","fever.5","fever.6","fever.7","fever.8","fever.9",
                        "hv_rscore.0","hv_rscore.1","hv_rscore.2","hv_rscore.3","hv_rscore.4","hv_rscore.5","hv_rscore.6","hv_rscore.7","hv_rscore.8","hv_rscore.9")]

# Create dataframe for daily symptom score (clinic day, day 0-9)

m0_score1 <- data.frame(hhID=symptom.index4$hhID)

m0_score1$clinic <- clinic7$b_m0_score

m0_score1$day14 <- m0_score1$day13 <- m0_score1$day12 <- m0_score1$day11 <- m0_score1$day10 <- m0_score1$day9 <-
  m0_score1$day8 <- m0_score1$day7 <- m0_score1$day6 <- m0_score1$day5 <- m0_score1$day4 <-
  m0_score1$day3 <- m0_score1$day2 <- m0_score1$day1 <- m0_score1$day0 <- NA

## adding onset to clinic delay
m0_score2 <-merge(m0_score1,clinic7[,c("hhID","onsettime")], by="hhID", all.x=T)
hvdelay <- hchar[c(1,6,7)]
hvdelay$hvdelay <- hvdelay$v1_day-hvdelay$clinic_day
m0_score2 <- merge(m0_score2,hvdelay[c(1,4)],all.x=TRUE)
m0_score <- m0_score2


### Imputing Total Symptoms Scores for Day 0 to Day 13

for (i in 1:nrow(m0_score)) {
    for (j in 1:9) {
      if (m0_score$onsettime[i]==1 & m0_score$hvdelay[i]==0) {
         m0_score$day0[i] <- m0_score$clinic[i]
         m0_score[i,j+3] <- m0_hv[i,j+2]
      }
      else if (m0_score$onsettime[i]==1 & m0_score$hvdelay[i]==1) {
         m0_score$day0[i] <- m0_score$clinic[i]
         m0_score[i,j+3] <- m0_hv[i,j+1]
      }
      else if (m0_score$onsettime[i]==1 & m0_score$hvdelay[i]==2) {
         m0_score$day0[i] <- m0_score$clinic[i]
         m0_score[i,j+4] <- m0_hv[i,j+1]
      }
      else if ((m0_score$onsettime[i]==2 | m0_score$onsettime[i]==3) & m0_score$hvdelay[i]==0) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$clinic[i]
         m0_score[i,j+4] <- m0_hv[i,j+2]
      }
      else if ((m0_score$onsettime[i]==2 | m0_score$onsettime[i]==3) & m0_score$hvdelay[i]==1) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$clinic[i]
         m0_score[i,j+4] <- m0_hv[i,j+1]
      }
      else if ((m0_score$onsettime[i]==2 | m0_score$onsettime[i]==3) & m0_score$hvdelay[i]==2) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$clinic[i]
         m0_score[i,j+5] <- m0_hv[i,j+1]
         m0_score[i,5] <- m0_hv[i,2]
      }
      else if (m0_score$onsettime[i]==4 & m0_score$hvdelay[i]==0) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$day2[i]<- m0_score$clinic[i]
         m0_score[i,j+5] <- m0_hv[i,j+2]
      }
      else if (m0_score$onsettime[i]==4 & m0_score$hvdelay[i]==1) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$day2[i]<- m0_score$clinic[i]
         m0_score[i,j+5] <- m0_hv[i,j+1]
      }
      else if (m0_score$onsettime[i]==4 & m0_score$hvdelay[i]==2) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$day2[i]<- m0_score$clinic[i]
         m0_score[i,j+6] <- m0_hv[i,j+1]
      }
      else if (m0_score$onsettime[i]==5 & m0_score$hvdelay[i]==0) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$day2[i]<- m0_score$day3[i]<-  m0_score$clinic[i]
         m0_score[i,j+6] <- m0_hv[i,j+2]
      }
      else if (m0_score$onsettime[i]==5 & m0_score$hvdelay[i]==1) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$day2[i]<- m0_score$day3[i]<-  m0_score$clinic[i]
         m0_score[i,j+6] <- m0_hv[i,j+1]
      }
      else if (m0_score$onsettime[i]==5 & m0_score$hvdelay[i]==2) {
         m0_score$day0[i]<- m0_score$day1[i]<- m0_score$day2[i]<- m0_score$day3[i]<-  m0_score$clinic[i]
         m0_score[i,j+7] <- m0_hv[i,j+1]
      }
    }
}

m0_score$all_event <- m0_score$all_time <- 0

for(i in 1:nrow(m0_score)){  #cycle through the rows 1 at a time
  for(j in 3:16){  #cycle through the columns 1 at a time
     if(!is.na(m0_score[i,j])){
      if(m0_score[i,j]>0){
        m0_score$all_time[i] <- j-3
      }
      if(m0_score[i,j]==0){
        m0_score$all_time[i] <- j-3
        m0_score$all_event[i] <- 1
        break
      }
    }
  }
}

output.all <- merge(clinic7, m0_score[c("hhID", "all_time", "all_event")])


#
# alleviation of fever -----------------------------------------------------------------------------------------------------
#

m0_fscore1 <- data.frame(hhID=symptom.index4$hhID)

m0_fscore1$f_clinic <- clinic7$fever

m0_fscore1$f_day14 <-m0_fscore1$f_day13 <-m0_fscore1$f_day12 <-m0_fscore1$f_day11 <- m0_fscore1$f_day10 <-
m0_fscore1$f_day9 <- m0_fscore1$f_day8 <- m0_fscore1$f_day7 <- m0_fscore1$f_day6 <- m0_fscore1$f_day5 <-
m0_fscore1$f_day4 <-  m0_fscore1$f_day3 <- m0_fscore1$f_day2 <- m0_fscore1$f_day1 <- m0_fscore1$f_day0 <- NA

## adding onset to clinic delay
m0_fscore2 <-merge(m0_fscore1,clinic7[,c("hhID","onsettime")], by="hhID", all.x=T)
m0_fscore2 <- merge(m0_fscore2,hvdelay[c(1,4)],all.x=TRUE)
m0_fscore <- m0_fscore2

###Imputing Day 0 to Day 14 fever score with adjustment of hvdelay and clinic delay from onset

for (i in 1:nrow(m0_fscore)) {
    for (j in 1:9) {
      if (m0_fscore$onsettime[i]==1 & m0_fscore$hvdelay[i]==0) {
         m0_fscore$f_day0[i] <- m0_fscore$f_clinic[i]
         m0_fscore[i,j+3] <- m0_hv[i,j+12]
      }
      else if (m0_fscore$onsettime[i]==1 & m0_fscore$hvdelay[i]==1) {
         m0_fscore$f_day0[i] <- m0_fscore$f_clinic[i]
         m0_fscore[i,j+3] <- m0_hv[i,j+11]
      }
      else if (m0_fscore$onsettime[i]==1 & m0_fscore$hvdelay[i]==2) {
         m0_fscore$f_day0[i] <- m0_fscore$f_clinic[i]
         m0_fscore[i,j+4] <- m0_hv[i,j+11]
      }
      else if ((m0_fscore$onsettime[i]==2 | m0_fscore$onsettime[i]==3) & m0_fscore$hvdelay[i]==0) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_clinic[i]
         m0_fscore[i,j+4] <- m0_hv[i,j+12]
      }
      else if ((m0_fscore$onsettime[i]==2 | m0_fscore$onsettime[i]==3) & m0_fscore$hvdelay[i]==1) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_clinic[i]
         m0_fscore[i,j+4] <- m0_hv[i,j+11]
      }
      else if ((m0_fscore$onsettime[i]==2 | m0_fscore$onsettime[i]==3) & m0_fscore$hvdelay[i]==2) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_clinic[i]
         m0_fscore[i,j+5] <- m0_hv[i,j+11]
         m0_fscore[i,5] <-m0_hv[i,6]
      }
      else if (m0_fscore$onsettime[i]==4 & m0_fscore$hvdelay[i]==0) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_day2[i]<- m0_fscore$f_clinic[i]
         m0_fscore[i,j+5] <- m0_hv[i,j+12]
      }
      else if (m0_fscore$onsettime[i]==4 & m0_fscore$hvdelay[i]==1) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_day2[i]<- m0_fscore$f_clinic[i]
         m0_fscore[i,j+5] <- m0_hv[i,j+11]
      }
      else if (m0_fscore$onsettime[i]==4 & m0_fscore$hvdelay[i]==2) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_day2[i]<- m0_fscore$f_clinic[i]
         m0_fscore[i,j+6] <- m0_hv[i,j+11]
      }
      else if (m0_fscore$onsettime[i]==5 & m0_fscore$hvdelay[i]==0) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_day2[i]<- m0_fscore$f_day3[i]<-  m0_fscore$f_clinic[i]
         m0_fscore[i,j+6] <- m0_hv[i,j+12]
      }
      else if (m0_fscore$onsettime[i]==5 & m0_fscore$hvdelay[i]==1) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_day2[i]<- m0_fscore$f_day3[i]<-  m0_fscore$f_clinic[i]
         m0_fscore[i,j+6] <- m0_hv[i,j+11]
      }
      else if (m0_fscore$onsettime[i]==5 & m0_fscore$hvdelay[i]==2) {
         m0_fscore$f_day0[i]<- m0_fscore$f_day1[i]<- m0_fscore$f_day2[i]<- m0_fscore$f_day3[i]<-  m0_fscore$f_clinic[i]
         m0_fscore[i,j+7] <- m0_hv[i,j+11]
      }
    }
}

## Loop determining alleviation

m0_fscore$f_event <- m0_fscore$f_time <- 0

for(i in 1:nrow(m0_fscore)){  #cycle through the rows 1 at a time
  for(j in 3:16){  #cycle through the columns 1 at a time
    if(!is.na(m0_fscore[i,j])){
      if(m0_fscore[i,j]>0){
        m0_fscore$f_time[i] <- j-3
      }
      if(m0_fscore[i,j]==0){
        m0_fscore$f_time[i] <- j-3
        m0_fscore$f_event[i] <- 1
        break
      }
    }
  }
}

output.fever <- merge(clinic7, m0_fscore[m0_fscore$f_day0>0,c("hhID", "f_time", "f_event")])

#
# alleviation of respiratory symptoms --------------------------------------------------------------------------------------
#

m0_rscore1 <- data.frame(hhID=symptom.index4$hhID)

m0_rscore1$r_clinic <- clinic7$b_m0_rscore

m0_rscore1$r_day14 <-m0_rscore1$r_day13 <-m0_rscore1$r_day12 <-m0_rscore1$r_day11 <- m0_rscore1$r_day10 <-
m0_rscore1$r_day9 <- m0_rscore1$r_day8 <- m0_rscore1$r_day7 <- m0_rscore1$r_day6 <- m0_rscore1$r_day5 <-
m0_rscore1$r_day4 <-  m0_rscore1$r_day3 <- m0_rscore1$r_day2 <- m0_rscore1$r_day1 <- m0_rscore1$r_day0 <- NA

## adding onset to clinic delay
m0_rscore2 <-merge(m0_rscore1,clinic7[,c("hhID","onsettime")], by="hhID", all.x=T)
m0_rscore2 <- merge(m0_rscore2,hvdelay[c(1,4)],all.x=TRUE)
m0_rscore <- m0_rscore2

for (i in 1:nrow(m0_rscore)) {
    for (j in 1:9) {
      if (m0_rscore$onsettime[i]==1 & m0_rscore$hvdelay[i]==0) {
         m0_rscore$r_day0[i] <- m0_rscore$r_clinic[i]
         m0_rscore[i,j+3] <- m0_hv[i,j+22]
      }
      else if (m0_rscore$onsettime[i]==1 & m0_rscore$hvdelay[i]==1) {
         m0_rscore$r_day0[i] <- m0_rscore$r_clinic[i]
         m0_rscore[i,j+3] <- m0_hv[i,j+21]
      }
      else if (m0_rscore$onsettime[i]==1 & m0_rscore$hvdelay[i]==2) {
         m0_rscore$r_day0[i] <- m0_rscore$r_clinic[i]
         m0_rscore[i,j+4] <- m0_hv[i,j+21]
      }
      else if ((m0_rscore$onsettime[i]==2 | m0_rscore$onsettime[i]==3) & m0_rscore$hvdelay[i]==0) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_clinic[i]
         m0_rscore[i,j+4] <- m0_hv[i,j+22]
      }
      else if ((m0_rscore$onsettime[i]==2 | m0_rscore$onsettime[i]==3) & m0_rscore$hvdelay[i]==1) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_clinic[i]
         m0_rscore[i,j+4] <- m0_hv[i,j+21]
      }
      else if ((m0_rscore$onsettime[i]==2 | m0_rscore$onsettime[i]==3) & m0_rscore$hvdelay[i]==2) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_clinic[i]
         m0_rscore[i,j+5] <- m0_hv[i,j+21]
         m0_rscore[i,5] <-m0_hv[i,6]
      }
      else if (m0_rscore$onsettime[i]==4 & m0_rscore$hvdelay[i]==0) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_day2[i]<- m0_rscore$r_clinic[i]
         m0_rscore[i,j+5] <- m0_hv[i,j+22]
      }
      else if (m0_rscore$onsettime[i]==4 & m0_rscore$hvdelay[i]==1) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_day2[i]<- m0_rscore$r_clinic[i]
         m0_rscore[i,j+5] <- m0_hv[i,j+21]
      }
      else if (m0_rscore$onsettime[i]==4 & m0_rscore$hvdelay[i]==2) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_day2[i]<- m0_rscore$r_clinic[i]
         m0_rscore[i,j+6] <- m0_hv[i,j+21]
      }
      else if (m0_rscore$onsettime[i]==5 & m0_rscore$hvdelay[i]==0) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_day2[i]<- m0_rscore$r_day3[i]<-  m0_rscore$r_clinic[i]
         m0_rscore[i,j+6] <- m0_hv[i,j+22]
      }
      else if (m0_rscore$onsettime[i]==5 & m0_rscore$hvdelay[i]==1) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_day2[i]<- m0_rscore$r_day3[i]<-  m0_rscore$r_clinic[i]
         m0_rscore[i,j+6] <- m0_hv[i,j+21]
      }
      else if (m0_rscore$onsettime[i]==5 & m0_rscore$hvdelay[i]==2) {
         m0_rscore$r_day0[i]<- m0_rscore$r_day1[i]<- m0_rscore$r_day2[i]<- m0_rscore$r_day3[i]<-  m0_rscore$r_clinic[i]
         m0_rscore[i,j+7] <- m0_hv[i,j+21]
      }
    }
}

## Loop picking up alleviation

m0_rscore$r_event <- m0_rscore$r_time <- 0

for(i in 1:nrow(m0_rscore)){  #cycle through the rows 1 at a time
  for(j in 3:17){  #cycle through the columns 1 at a time
    if(!is.na(m0_rscore[i,j])){
      if(m0_rscore[i,j]>0){
        m0_rscore$r_time[i] <- j-3
      }
      if(m0_rscore[i,j]==0){
        m0_rscore$r_time[i] <- j-3
        m0_rscore$r_event[i] <- 1
        break
      }
    }
  }
}

output.rs <- merge(clinic7, m0_rscore[m0_rscore$r_day0>0,c("hhID", "r_time", "r_event")])


# Duration of viral shedding (interval censored) -----------------------------------------------------------------------

hc.pcr <- reshape(hc[hc$member==0,],timevar="visit",idvar=c("hhID","member"),direction="wide")
hc.pcr <- merge(hc.pcr,hchar[c(1,6:10)],all.x=TRUE)
vshed <- hc.pcr[hc.pcr$hhID%in%clinic7$hhID,]

vshed$mark <- vshed$timeR <- vshed$timeL <-  rep(NA,nrow(vshed))

for (i in 1:nrow(vshed)){
    if (!is.na(vshed$PCR.4[i])&vshed$PCR.4[i]==1){
       vshed$timeL[i] <- vshed$v4_day[i]
       vshed$mark[i] <- 4
    }
    else if(!is.na(vshed$PCR.3[i])&vshed$PCR.3[i]==1){
       vshed$timeL[i] <- vshed$v3_day[i]
       vshed$mark[i] <- 3
    }
    else if(!is.na(vshed$PCR.2[i])&vshed$PCR.2[i]==1){
       vshed$timeL[i] <- vshed$v2_day[i]
       vshed$mark[i] <- 2
    }
    else if(!is.na(vshed$PCR.1[i])&vshed$PCR.1[i]==1){
       vshed$timeL[i] <- vshed$v1_day[i]
       vshed$mark[i] <- 1
    }
    else if(!is.na(vshed$PCR.0[i])&vshed$PCR.0[i]==1){
       vshed$timeL[i] <- vshed$clinic_day[i]
       vshed$mark[i] <- 0
    }
}

for (i in 1:nrow(vshed)){
     if(!is.na(vshed$PCR.1[i])&vshed$mark[i]<1) vshed$timeR[i] <- vshed$v1_day[i]
     else if(!is.na(vshed$PCR.2[i])&vshed$mark[i]<2) vshed$timeR[i] <- vshed$v2_day[i]
     else if(!is.na(vshed$PCR.3[i])&vshed$mark[i]<3) vshed$timeR[i] <- vshed$v3_day[i]
     else if(!is.na(vshed$PCR.4[i])&vshed$mark[i]<4) vshed$timeR[i] <- vshed$v4_day[i]
     else vshed$timeR[i] <- 100
}

vshed <- vshed[c("hhID","member","timeL","timeR")]

# Define secondary cases --------------------------------------------------------------------------------------------------

pcr.all <- reshape(hc,timevar="visit",idvar=c("hhID","member"),direction="wide")
trans <- merge(demog[c(1:5,12)],pcr.all,by=c("hhID","member"),all.x=TRUE)
trans <- trans[trans$hhID%in%index.hh$hhID,]
trans.c <- trans[trans$member>0,]

for (i in 1:nrow(trans.c)){
  if( !is.na(trans.c$PCR.1[i])&trans.c$PCR.1[i]==0 & ( (!is.na(trans.c$PCR.2[i])&trans.c$PCR.2[i]==1) 
       | (!is.na(trans.c$PCR.3[i])&trans.c$PCR.3[i]==1) | (!is.na(trans.c$PCR.4[i])&trans.c$PCR.4[i]==1) ) )  trans.c$labsedcase[i] <- 1 
  else  trans.c$labsedcase[i] <- 0 
}

##
symptom$fever <- 1*(symptom$bodytemp>=37.8)
symptom$indicate <- symptom$fever+symptom$cough+symptom$headache+symptom$sthroat+symptom$pmuscle
symptom$flu <- 1*(symptom$indicate>=2)
symptom.temp <- reshape(symptom[c(1:3,13)], timevar="day", idvar=c("hhID","member"), direction="wide", v.names="flu")
trans.c2 <- merge(trans.c[c(1,2,7)],symptom.temp, by=c("hhID","member"), all.x=TRUE)
names(trans.c2) <- c("hhID","member","PCR.0","day0","day1","day2","day3","day4","day5","day6","day7","day8","day9")

## Define secondary cases
for (i in 1:nrow(trans.c2)){
    if ( !is.na(trans.c$PCR.1[i]) & trans.c$PCR.1[i]==0 &
         ( (!is.na(trans.c2$day1[i]) & trans.c2$day1[i]==1) | (!is.na(trans.c2$day2[i]) & trans.c2$day2[i]==1) | 
       (!is.na(trans.c2$day3[i]) & trans.c2$day3[i]==1) | (!is.na(trans.c2$day4[i]) & trans.c2$day4[i]==1) | 
       (!is.na(trans.c2$day5[i]) & trans.c2$day5[i]==1) | (!is.na(trans.c2$day6[i]) & trans.c2$day6[i]==1) | 
       (!is.na(trans.c2$day7[i]) & trans.c2$day7[i]==1) | (!is.na(trans.c2$day8[i]) & trans.c2$day8[i]==1) | 
       (!is.na(trans.c2$day9[i]) & trans.c2$day9[i]==1) ))     trans.c2$clinicsedcase[i] <- 1
     else trans.c2$clinicsedcase[i] <- 0
}
trans.c$clinicsedcase <- trans.c2$clinicsedcase

#
# End of script
#

