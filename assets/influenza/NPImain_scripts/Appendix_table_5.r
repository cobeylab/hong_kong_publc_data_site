#
# Script to reproduce information in Appendix Table 5 from:
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
symptom <- read.csv(paste(dir, "symptomday_d.csv", sep=""))

appt5 <- matrix(rep(NA,80),ncol=10,byrow=FALSE)

## Define lab-confirmed secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$analyzed[i]==1 & hculture$exd_contact[i]==0 &
             ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) | (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) ) )
         hculture$labsedcase[i] <- 1
	  else hculture$labsedcase[i] <- 0
}

## Define clinic-confirmed secondary cases

for (k in 1:2){

if(k==1){
    # Clinical definition 1
    symptom$fever <- 1*(!is.na(symptom$bodytemp)&symptom$bodytemp>=37.8)
    symptom$indicate <- symptom$fever+symptom$cough+symptom$headache+symptom$sthroat+symptom$pmuscle
    symptom$flu <- 1*(symptom$indicate>=2)
}

if(k==2){
    # Clinical definition 2
   symptom$fever <- 1*(!is.na(symptom$bodytemp)&symptom$bodytemp>=37.8)
   symptom$indicate <- symptom$cough+symptom$sthroat
   symptom$flu <- 1*(symptom$fever==1 & symptom$indicate>=1)
}

##### Followed by either one clinical flu definition...

hclinic <- data.frame(hhID = hculture$hhID)
hclinic$member <- 0
for(i in 2:nrow(hclinic)){
  if(hclinic$hhID[i]==hclinic$hhID[i-1]) hclinic$member[i] <- hclinic$member[i-1]+1
}

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
if(k==1) hculture$c2nd1 <- hclinic$clinicsedcase
else if (k==2) hculture$c2nd2 <- hclinic$clinicsedcase
}

### Add variable 'intervention' and 'indexdelay' (time from symptom onset to intervention)
hculture <- merge(hculture[c(1,2,9:12)],housechar[,c(1,2,9)])
hculture$arm <- factor(hculture$intervention, levels=c(1,3,4), labels=c("control", "hand", "handmask"))
hculture$delay[hculture$onset_v1_delay<=36] <- "d1"
hculture$delay[hculture$onset_v1_delay>36] <- "d2"

# Four definitions for influenza
hculture$def1 <- 1*(hculture$labsedcase==1|hculture$c2nd1==1)
hculture$def2 <- 1*(hculture$labsedcase==1&hculture$c2nd1==1)
hculture$def3 <- 1*(hculture$labsedcase==1|hculture$c2nd2==1)
hculture$def4 <- 1*(hculture$labsedcase==1&hculture$c2nd2==1)

hculture2 <- hculture[-c(4:6)]

#-------------------------------------- Calculate lab-confirmed SAR and CI ---------------------------------------------#

hculture.cc <- merge(hculture2,housechar[,c(1,3)])
names(hculture.cc)[dim(hculture.cc)[2]]<- "hhsize"
hculture.cc <- hculture.cc[hculture.cc$member!=0&hculture.cc$analyzed==1,]
hculture.cc$hhsize <- hculture.cc$hhsize -1
hculture.cc$hhsize[hculture.cc$member>1] <- 0
hc_contact <- hculture.cc[c(1,2,6:12)]

## Begin function 'icc' for intracluster correlation coefficient
icc <- function(data_id,data_arm,data_clustersize,data_success){
 M <- length(data_id)
 Y <- sum(data_success)
 P <- Y/M
 K <- length(unique(data_id))

 arm <- c("control","hand","handmask")   # interventions
 Yij <- Mij <- Pij <- hhID <- matrix(nrow=3,ncol=max(length(unique(data_id[data_arm==arm[1]])),
                                                     length(unique(data_id[data_arm==arm[2]])),
                                                     length(unique(data_id[data_arm==arm[3]])) ))
 for (i in 1:length(arm)){
      hhID[i,1:length(unique(data_id[data_arm==arm[i]]))] <- unique(data_id[data_arm==arm[i]])
 }

 for (i in 1:length(arm)){
      for (j in 1:length(unique(data_id[data_arm==arm[i]]))){
        Yij[i,j] <- sum(data_success[data_id==hhID[i,j]])
        Mij[i,j] <- sum(data_clustersize[data_id==hhID[i,j]])
        Pij[i,j] <- Yij[i,j]/Mij[i,j]
      }
 }

 Y_i <- M_i <- P_i <- m_Ai <- rho_i <- rep(NA,length(arm))
 output <- data.frame(Y_i,M_i,m_Ai,P_i,rho_i)
 for (i in 1:length(arm)){
     output$Y_i[i] <- sum(na.exclude(Yij[i,]))
     output$M_i[i] <- sum(na.exclude(Mij[i,]))
     output$P_i[i] <- output$Y_i[i]/output$M_i[i]
     output$m_Ai[i] <- sum(na.exclude(Mij[i,]^2))/output$M_i[i]
 }

 m_0 <- (M-sum(output$m_Ai))/(K-2)
 MSW <- sum(rowSums(Mij*Pij*(1-Pij),na.rm = TRUE)) / (M-K)
 MSC <- sum(rowSums(Mij*(Pij-output$P_i)*(Pij-output$P_i),na.rm = TRUE)) / (K-2)
 rho <- (MSC-MSW)/(MSC+(m_0-1)*MSW)
 output$rho_i <- rho
 output
}
## End of function 'icc'

## Function cluster bootstrap CI
cbCI <- function(data,def,nboot){
   idlist <- unique(data$hhID)
   
   if(def==1) sed <- sum(data$def1)
   else if (def==2) sed <- sum(data$def2)
   else if (def==3) sed <- sum(data$def3)
   else sed <- sum(data$def4)
   
   est <- sed/dim(data)[1]  # estimated SAR
   sar <- rep(NA,nboot)
   for (i in 1:nboot){
	 idsample <- sort(sample(idlist,replace=TRUE))
	 newdata <- data[data$hhID==idsample[1],]
	 for(j in 2:length(idsample)){
            temp <- data[data$hhID==idsample[j],]
	    newdata <- rbind(newdata,temp)
	 }
	 if(def==1) sar[i] <- sum(newdata$def1)/dim(newdata)[1]
	 else if(def==2) sar[i] <- sum(newdata$def2)/dim(newdata)[1]
	 else if(def==3) sar[i] <- sum(newdata$def3)/dim(newdata)[1]
	 else sar[i] <- sum(newdata$def4)/dim(newdata)[1]
   }
   round(c(est,quantile(sar,c(0.025,0.975))),2)
}
## End of function

################################################## Full data #############################################################

c(dim(hc_contact[hc_contact$arm=="control",])[1],dim(hc_contact[hc_contact$arm=="hand",])[1],dim(hc_contact[hc_contact$arm=="handmask",])[1])

# SAR & cluster bootstrap CIs
c2nd_all <- hc_contact[hc_contact$arm=="control",]
h2nd_all <- hc_contact[hc_contact$arm=="hand",]
m2nd_all <- hc_contact[hc_contact$arm=="handmask",]
set.seed(12345)
for(i in 1:4){
   appt5[i,1:9] <- c(cbCI(c2nd_all,i,1000), cbCI(h2nd_all,i,1000), cbCI(m2nd_all,i,1000))
}

## Chi-square test
for(j in 1:4){
if(j==1) chi_stat <- icc(hc_contact$hhID,hc_contact$arm,hc_contact$hhsize,hc_contact$def1)
if(j==2) chi_stat <- icc(hc_contact$hhID,hc_contact$arm,hc_contact$hhsize,hc_contact$def2)
if(j==3) chi_stat <- icc(hc_contact$hhID,hc_contact$arm,hc_contact$hhsize,hc_contact$def3)
if(j==4) chi_stat <- icc(hc_contact$hhID,hc_contact$arm,hc_contact$hhsize,hc_contact$def4)

C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
appt5[j,10] <- round(pchisq(achi_lab,df=2,lower.tail=FALSE),3)
}
## End of chi-square test

######################## Repeat above stats stratified by symptom onset-intervention <= 36 hrs ###########################

c(dim(hc_contact[hc_contact$arm=="control" & hc_contact$delay=="d1",])[1],
  dim(hc_contact[hc_contact$arm=="hand" & hc_contact$delay=="d1",])[1],
  dim(hc_contact[hc_contact$arm=="handmask" & hc_contact$delay=="d1",])[1])

# SAR & cluster bootstrap CIs
c2nd_s <- hc_contact[hc_contact$arm=="control" & hc_contact$delay=="d1",]
h2nd_s <- hc_contact[hc_contact$arm=="hand" & hc_contact$delay=="d1",]
m2nd_s <- hc_contact[hc_contact$arm=="handmask" & hc_contact$delay=="d1",]
set.seed(12345)
for(i in 1:4){
   appt5[i+4,1:9] <- c(cbCI(c2nd_s,i,1000), cbCI(h2nd_s,i,1000), cbCI(m2nd_s,i,1000))
}

## Chi-square test
hc_lt <- hc_contact[hc_contact$delay=="d1",]
for(j in 1:4){
if(j==1) chi_stat <- icc(hc_lt$hhID,hc_lt$arm,hc_lt$hhsize,hc_lt$def1)
if(j==2) chi_stat <- icc(hc_lt$hhID,hc_lt$arm,hc_lt$hhsize,hc_lt$def2)
if(j==3) chi_stat <- icc(hc_lt$hhID,hc_lt$arm,hc_lt$hhsize,hc_lt$def3)
if(j==4) chi_stat <- icc(hc_lt$hhID,hc_lt$arm,hc_lt$hhsize,hc_lt$def4)

C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
appt5[j+4,10] <- round(pchisq(achi_lab,df=2,lower.tail=FALSE),3)
}
## End of chi-square test

rownames(appt5) <- c("lab or cdef1","lab and cdef1","lab or cdef2","lab and cdef2",
                      "<36 - lab or cdef1","<36 - lab and cdef1","<36 - lab or cdef2","<36 - lab and cdef2")
colnames(appt5) <- c("Control","CI-low","CI-up","Hany hygiene","CI-low","CI-up","Mask+HH","CI-low","CI-up","p-value")
appt5

# End of script


