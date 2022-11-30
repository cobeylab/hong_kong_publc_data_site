#
# Script to reproduce information in Table 2 (secondary attack ratios) from:
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

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

#source the Figure_1 file to avoid repeating some of the data cleaning procedures
source("http://www.hku.hk/bcowling/influenza/NPI_scripts/Figure_1.r")  

table2 <- matrix(rep(NA,120),ncol=10,byrow=FALSE)
baseflu <- baseflu[,1:4]         # hhID, member, age, sex
analysis <- analysis[,c(1,2,9)]  # hhID, intervention, baseline
tab2 <- merge(baseflu,analysis)
tab2_als <- tab2[tab2$baseline==0,]

tab2_c <- tab2_als[tab2_als$intervention==1,]
tab2_m <- tab2_als[tab2_als$intervention==2,]
tab2_h <- tab2_als[tab2_als$intervention==3,]

############################ For secondary cases calculation ##############################################################################

## Begin function 'icc' for intracluster correlation coefficient
icc <- function(data_id,data_arm,data_clustersize,data_success){
 M <- length(data_id)
 Y <- sum(data_success)
 P <- Y/M
 K <- length(unique(data_id))

 arm <- c("control","mask","hand")   # interventions
 Yij <- Mij <- Pij <- hhID <- matrix(nrow=3,ncol=length(unique(data_id[data_arm==arm[1]])))
 for (i in 1:length(arm)){
      hhID[i,1:length(unique(data_id[data_arm==arm[i]]))] <- as.character(unique(data_id[data_arm==arm[i]]))
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

###---------------------------------------------- Lab confirmed ------------------------------------------------------------------------###

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
###

hchar <- merge(hchar,qv[,c(2,25)],by="hhID",all.x=TRUE)
hchar$onsettime[hchar$onsettime!=5] <- floor(hchar$onsettime[hchar$onsettime!=5]/2)
hchar$onsettime[hchar$onsettime==5] <- 3
hchar$v1_day <- hchar$v1_day+hchar$onsettime # from symptom onset to intervention
hculture <- merge(hculture,hchar[,c(1:3,5)],by="hhID",all.x=TRUE)
hculture$intervention <- factor(hculture$intervention, levels=1:3, labels=c("control", "mask", "hand"))
hculture$indexdelay <- 1*(hculture$v1_day<=1)
names(hculture)[14]<- "hhsize"
hculture$hhsize <- hculture$hhsize -1
hculture$hhsize[hculture$member!=1] <- 0

hc_contact <- hculture[hculture$member>0&hculture$d_index==0,c(1,2,12:14,16)]

########################################################### Full data #####################################################################

c2nd_all <- dim(hc_contact[hc_contact$intervention=="control",])[1]
m2nd_all <- dim(hc_contact[hc_contact$intervention=="mask",])[1]
h2nd_all <- dim(hc_contact[hc_contact$intervention=="hand",])[1]
c2nd_all 
m2nd_all 
h2nd_all 

c2nd <- dim(hc_contact[hc_contact$intervention=="control"&hc_contact$labsedcase==1,])[1]
m2nd <- dim(hc_contact[hc_contact$intervention=="mask"&hc_contact$labsedcase==1,])[1]
h2nd <- dim(hc_contact[hc_contact$intervention=="hand"&hc_contact$labsedcase==1,])[1]

## Adjusted chi-square test 
chi_stat <- icc(hc_contact$hhID,hc_contact$intervention,hc_contact$hhsize,hc_contact$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
chi <- round(pchisq(achi_lab,df=2,lower.tail=FALSE),2)
## End of chi-square test

table2[1,] <- c(round(c2nd/c2nd_all,2), round(binom.test(c2nd,c2nd_all)$conf[1:2],2),
        round(m2nd/m2nd_all,2), round(binom.test(m2nd,m2nd_all)$conf[1:2],2),
        round(h2nd/h2nd_all,2), round(binom.test(h2nd,h2nd_all)$conf[1:2],2), chi)

############################## Repeat above stats stratified by symptom onset-intervention <= 36 hrs #####################################

c2nd_all_s <- dim(hc_contact[hc_contact$intervention=="control" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
m2nd_all_s <- dim(hc_contact[hc_contact$intervention=="mask" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
h2nd_all_s <- dim(hc_contact[hc_contact$intervention=="hand" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
c2nd_all_s 
m2nd_all_s 
h2nd_all_s 

c2nd <- dim(hc_contact[hc_contact$intervention=="control"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
m2nd <- dim(hc_contact[hc_contact$intervention=="mask"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]
h2nd <- dim(hc_contact[hc_contact$intervention=="hand"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,])[1]

## Adjusted chi-square test 
hc.lt <- hc_contact[!is.na(hc_contact$indexdelay) & hc_contact$indexdelay==1,] 
chi_stat <- icc(hc.lt$hhID,hc.lt$intervention,hc.lt$hhsize,hc.lt$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
chi <- round(pchisq(achi_lab,df=2,lower.tail=FALSE),2)
## End of chi-square test

table2[5,] <- c(round(c2nd/c2nd_all_s,2), round(binom.test(c2nd,c2nd_all_s)$conf[1:2],2),
        round(m2nd/m2nd_all_s,2), round(binom.test(m2nd,m2nd_all_s)$conf[1:2],2),
        round(h2nd/h2nd_all_s,2), round(binom.test(h2nd,h2nd_all_s)$conf[1:2],2), chi)

############################### Repeat above stats stratified by symptom onset-intervention > 36 hrs ######################################

c2nd_all_g <- dim(hc_contact[hc_contact$intervention=="control" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
m2nd_all_g <- dim(hc_contact[hc_contact$intervention=="mask" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
h2nd_all_g <- dim(hc_contact[hc_contact$intervention=="hand" & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
c2nd_all_g 
m2nd_all_g 
h2nd_all_g 

c2nd <- dim(hc_contact[hc_contact$intervention=="control"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
round(c2nd/c2nd_all_g,2)
round(binom.test(c2nd,c2nd_all_g)$conf[1:2],2)

m2nd <- dim(hc_contact[hc_contact$intervention=="mask"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
round(m2nd/m2nd_all_g,2)
round(binom.test(m2nd,m2nd_all_g)$conf[1:2],2)

h2nd <- dim(hc_contact[hc_contact$intervention=="hand"&hc_contact$labsedcase==1 & !is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,])[1]
round(h2nd/h2nd_all,2)
round(binom.test(h2nd,h2nd_all_g)$conf[1:2],2)

## Adjusted chi-square test 
hc.gt <- hc_contact[!is.na(hc_contact$indexdelay) & hc_contact$indexdelay==0,] 
chi_stat <- icc(hc.gt$hhID,hc.gt$intervention,hc.gt$hhsize,hc.gt$labsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_lab <- sum(chi)
chi <- round(pchisq(achi_lab,df=2,lower.tail=FALSE),2)
## End of chi-square test

table2[9,] <- c(round(c2nd/c2nd_all_g,2), round(binom.test(c2nd,c2nd_all_g)$conf[1:2],2),
        round(m2nd/m2nd_all_g,2), round(binom.test(m2nd,m2nd_all_g)$conf[1:2],2),
        round(h2nd/h2nd_all_g,2), round(binom.test(h2nd,h2nd_all_g)$conf[1:2],2), chi)


###---------------------------------------------- 3 clinic definitions -----------------------------------------------------------------###

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

hclinic$intervention <- hculture$intervention
hclinic$d_index <- hculture$d_index
hclinic$indexdelay <- hculture$indexdelay
hclinic$hhsize <- hculture$hhsize
clicon <- hclinic[hclinic$member>0&hclinic$d_index==0,]

########################################################### Full data ####################################################################

c2nd <- dim(clicon[clicon$intervention=="control"&clicon$clinicsedcase==1,])[1]
m2nd <- dim(clicon[clicon$intervention=="mask"&clicon$clinicsedcase==1,])[1]
h2nd <- dim(clicon[clicon$intervention=="hand"&clicon$clinicsedcase==1,])[1]

## Adjusted chi-square test 
chi_stat <- icc(clicon$hhID,clicon$intervention,clicon$hhsize,clicon$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_clinic <- sum(chi)
chi <- round(pchisq(achi_clinic,df=2,lower.tail=FALSE),2)
## End of chi-square test

table2[1+k,] <- c(round(c2nd/c2nd_all,2), round(binom.test(c2nd,c2nd_all)$conf[1:2],2),
        round(m2nd/m2nd_all,2), round(binom.test(m2nd,m2nd_all)$conf[1:2],2),
        round(h2nd/h2nd_all,2), round(binom.test(h2nd,h2nd_all)$conf[1:2],2), chi)

#################################### Repeat above stats stratified by symptom onset-intervention <= 36 hrs ################################

c2nd <- dim(clicon[clicon$intervention=="control"&clicon$clinicsedcase==1 & clicon$indexdelay==1,])[1]
m2nd <- dim(clicon[clicon$intervention=="mask"&clicon$clinicsedcase==1 & clicon$indexdelay==1,])[1]
h2nd <- dim(clicon[clicon$intervention=="hand"&clicon$clinicsedcase==1 & clicon$indexdelay==1,])[1]

## Adjusted chi-square test 
clicon.lt <- clicon[!is.na(clicon$indexdelay) & clicon$indexdelay==1,] 
chi_stat <- icc(clicon.lt$hhID,clicon.lt$intervention,clicon.lt$hhsize,clicon.lt$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_clinic <- sum(chi)
chi <- round(pchisq(achi_clinic,df=2,lower.tail=FALSE),2)
## End of chi-square test

table2[5+k,] <- c(round(c2nd/c2nd_all_s,2), round(binom.test(c2nd,c2nd_all_s)$conf[1:2],2),
        round(m2nd/m2nd_all_s,2), round(binom.test(m2nd,m2nd_all_s)$conf[1:2],2),
        round(h2nd/h2nd_all_s,2), round(binom.test(h2nd,h2nd_all_s)$conf[1:2],2), chi)

###################################### Repeat above stats stratified by symptom onset-intervention > 36 hrs ###############################

c2nd <- dim(clicon[clicon$intervention=="control"&clicon$clinicsedcase==1 & clicon$indexdelay==0,])[1]
m2nd <- dim(clicon[clicon$intervention=="mask"&clicon$clinicsedcase==1 & clicon$indexdelay==0,])[1]
h2nd <- dim(clicon[clicon$intervention=="hand"&clicon$clinicsedcase==1 & clicon$indexdelay==0,])[1]

## Adjusted chi-square test 
clicon.gt <- clicon[!is.na(clicon$indexdelay) & clicon$indexdelay==0,] 
chi_stat <- icc(clicon.gt$hhID,clicon.gt$intervention,clicon.gt$hhsize,clicon.gt$clinicsedcase)
C_i <- chi <- rep(NA,3)
P <- sum(chi_stat$Y_i)/sum(chi_stat$M_i)
for(i in 1:3){
    C_i[i] <- 1+(chi_stat$m_Ai[i])*chi_stat$rho_i[1]
    chi[i] <- chi_stat$M_i[i]*(chi_stat$P_i[i]-P)^2/(C_i[i]*P*(1-P))
}
achi_clinic <- sum(chi)
chi <- round(pchisq(achi_clinic,df=2,lower.tail=FALSE),2)
## End of chi-square test

table2[9+k,] <- c(round(c2nd/c2nd_all_g,2), round(binom.test(c2nd,c2nd_all_g)$conf[1:2],2),
        round(m2nd/m2nd_all_g,2), round(binom.test(m2nd,m2nd_all_g)$conf[1:2],2),
        round(h2nd/h2nd_all_g,2), round(binom.test(h2nd,h2nd_all_g)$conf[1:2],2), chi)
}

rownames(table2) <- c("Any-lab","Any-cd1","Any-cd2","Any-cd3","<36-lab","<36-cd1","<36-cd2","<36-cd3",">36-lab",">36-cd1",">36-cd2",">36-cd3")
colnames(table2) <- c("control","CI-low","CI-up","mask","CI-low","CI-up","hand","CI-low","CI-up","p-value")
table2


# End of script
