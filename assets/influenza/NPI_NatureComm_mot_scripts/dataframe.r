#
# R syntax to reformat raw data for:
#
# Cowling BJ, Ip DKM, Fang VJ, et al.
# Aerosol transmission is an important mode of influenza A virus spread
# Nature Communications, 2013 (in press).
#
# Last updated by ang VJ and Cowling BJ.
# March 22, 2013
#


#
# NPI 2008 data
#


source("http://www.hku.hk/bcowling/influenza/NPImain_scripts/Analyzed_hh.r")
dir <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
demog <- read.csv(paste(dir, "adherence_m.csv", sep=""))
symp <- read.csv(paste(dir, "symptomday_d.csv", sep=""))

hc <- merge(hc,clinic[c("hhID","QVres")],all.x=T); hc <- hc[order(hc$hhID,hc$member,hc$visit),]
fluA <- unique(hc$hhID[hc$member==0&hc$QVres==1&!is.na(hc$QVres)])
hculture <- hculture[hculture$hhID%in%fluA,]

hculture$V1[is.na(hculture$V1)] <- hculture$V2[is.na(hculture$V2)] <- hculture$V3[is.na(hculture$V3)] <- 0
hculture <- hculture[hculture$analyzed==1,]
hculture$infect <- 1*((hculture$member>0&(hculture$V2==1|hculture$V3==1))|hculture$member==0)
symp$fever <- 1*(symp$bodytemp>=37.8)
symp$ILI <- 1*(symp$bodytemp>=37.8&symp$cough==1)
symp$symp <- 1*(rowSums(symp[,c(5:11)],na.rm=T)>=1)
symp2a <- reshape(symp[c("hhID","member","day","fever")],timevar="day",idvar=c("hhID","member"),v.names="fever",direction="wide")
symp2b <- reshape(symp[c("hhID","member","day","cough")],timevar="day",idvar=c("hhID","member"),v.names="cough",direction="wide")
symp2a$fever <- 1*(rowSums(symp2a[,3:12],na.rm=T)>=1)
symp2b$cough <- 1*(rowSums(symp2b[,3:12],na.rm=T)>=1)
symp3 <- merge(symp2a[c("hhID","member","fever")],symp2b[c("hhID","member","cough")],by=c("hhID","member"))
symp3$ILI <- 1*(symp3$fever==1&symp3$cough==1)

tmp <- merge(hculture,symp3[c("hhID","member","ILI")],by=c("hhID","member"),all.x=T)

# find symptom onset day for secondary cases - T defined as min(onset,PCR+ day)-1
symp4 <- unique(symp[1:2]); symp4$onset <- NA
for(i in 1:nrow(symp4)){
  ind <- symp$day[symp$hhID==symp4$hhID[i]&symp$member==symp4$member[i]&!is.na(symp$symp)&symp$symp==1]
  if(length(ind)!=0) symp4$onset[i] <- min(ind)
}

# infect=0&ILI=1: infection not detected by PCR? or ILI from other non-flu infections?
tmp$ILI[tmp$infect==0&tmp$ILI==1] <- 0

mdata <- merge(tmp[tmp$analyzed==1&tmp$member>0,],hchar[c(1:2,7:9)],all.x=T)
mdata$delta <- mdata$infect
mdata <- merge(mdata,symp4,by=c("hhID","member"),all.x=T); mdata$onset <- mdata$onset+mdata$v1_day
incubation <- 1
for(i in 1:nrow(mdata)){
  if(mdata$V2[i]==0&mdata$V3[i]==0) mdata$T0[i] <- mdata$v3_day[i]
  if(mdata$V2[i]==1) mdata$T0[i] <- min(mdata$v2_day[i],mdata$onset[i],na.rm=T)-incubation
  if(mdata$V2[i]==0&mdata$V3[i]==1) mdata$T0[i] <- min(mdata$v3_day[i],mdata$onset[i],na.rm=T)-incubation
}
mdata$T0 <- pmax(mdata$v1_day+1,mdata$T0)

# construct the model data frame
cr <- mdata[c("hhID","member","T0","delta","v1_day")]  #v1_day can be treated as truncation day
names(cr)[5] <- "trunc_day";
cr$X_hh <- 1*(mdata$intervention==3|mdata$intervention==4)
cr$X_fm <- 1*(mdata$intervention==4)
cr <- merge(cr,tmp[c("hhID","member","ILI")],by=c("hhID","member"),all.x=T)
cr <- merge(cr,demog[c("hhID","member","age")],by=c("hhID","member"),all.x=T)

#
# NPI 2009 winter data
#

dir <- "http://sph.hku.hk/data/HongKongNPIstudy2009V1/"
hc <- read.csv(paste(dir, "qPCR.csv", sep=""))
demog <- read.csv(paste(dir, "demog_m.csv", sep=""))
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
symp <- read.csv(paste(dir, "symp_d.csv", sep=""))
clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))

hc$flu <- 1*(hc$qPCR>0)
hc <- merge(hc,clinic[c("hhID","QVres")],all.x=T); hc <- hc[order(hc$hhID,hc$member,hc$visit),]
fluA <- unique(hc$hhID[hc$member==0&((hc$swine==1&!is.na(hc$swine))|(hc$AH1==1&!is.na(hc$AH1))|(hc$AH3==1&!is.na(hc$AH3))|(hc$QVres==1&!is.na(hc$QVres)))])
hculture <- reshape(hc[hc$hhID%in%fluA,c(c("hhID","member","visit","flu"))], timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="flu")
names(hculture) <- c("hhID","member","V1","V2","V3")
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
hculture <- merge(hculture[,-6], exd_index)
hculture <- merge(hculture[,-6], d_contact,all.x=TRUE)
hculture$d_contact[is.na(hculture$d_contact)] <- 0
hculture <- hculture[order(hculture$hhID,hculture$member),]
hculture$analyzed <- 1*(hculture$exd_index==0&hculture$d_contact==0)

hculture$V1[is.na(hculture$V1)] <- hculture$V2[is.na(hculture$V2)] <- hculture$V3[is.na(hculture$V3)] <- 0
hculture <- hculture[hculture$analyzed==1,]
hculture$infect <- 1*((hculture$member>0&(hculture$V2==1|hculture$V3==1))|hculture$member==0)

symp$fever <- 1*(symp$bodytemp>=37.8)
symp$symp <- 1*(rowSums(symp[,c(5:11)],na.rm=T)>=1)
symp2a <- reshape(symp[c("hhID","member","day","fever")],timevar="day",idvar=c("hhID","member"),v.names="fever",direction="wide")
symp2b <- reshape(symp[c("hhID","member","day","cough")],timevar="day",idvar=c("hhID","member"),v.names="cough",direction="wide")
symp2a$fever <- 1*(rowSums(symp2a[,3:12],na.rm=T)>=1)
symp2b$cough <- 1*(rowSums(symp2b[,3:12],na.rm=T)>=1)
symp3 <- merge(symp2a[c("hhID","member","fever")],symp2b[c("hhID","member","cough")],by=c("hhID","member"))
symp3$ILI <- 1*(symp3$fever==1&symp3$cough==1)
tmp <- merge(hculture,symp3[c("hhID","member","ILI")],by=c("hhID","member"),all.x=T)

# find symptom onset day for secondary cases - T defined as min(onset,PCR+ day)-1
symp4 <- unique(symp[1:2]); symp4$onset <- NA
for(i in 1:nrow(symp4)){
  ind <- symp$day[symp$hhID==symp4$hhID[i]&symp$member==symp4$member[i]&!is.na(symp$symp)&symp$symp==1]
  if(length(ind)!=0) symp4$onset[i] <- min(ind)
}

# infect=0&ILI=1: infection not detected by PCR? or ILI from other non-flu infections?
tmp$ILI[tmp$infect==0&tmp$ILI==1&!is.na(tmp$ILI)] <- 0

mdata <- merge(tmp[tmp$analyzed==1&tmp$member>0,],hchar[c(1,2,7:9)],all.x=T); mdata <- mdata[!is.na(mdata$v1_day)&!is.na(mdata$v2_day),]
mdata$delta <- mdata$infect
mdata <- merge(mdata,symp4,by=c("hhID","member"),all.x=T); mdata$onset <- mdata$onset+mdata$v1_day
incubation <- 1
for(i in 1:nrow(mdata)){
  if(mdata$V2[i]==0&mdata$V3[i]==0) mdata$T0[i] <- mdata$v3_day[i]
  if(mdata$V2[i]==1) mdata$T0[i] <- min(mdata$v2_day[i],mdata$onset[i],na.rm=T)-incubation
  if(mdata$V2[i]==0&mdata$V3[i]==1) mdata$T0[i] <- min(mdata$v3_day[i],mdata$onset[i],na.rm=T)-incubation
}
mdata$T0 <- pmax(mdata$v1_day+1,mdata$T0)

# construct the model data frame
cr2 <- mdata[c("hhID","member","T0","delta","v1_day")]  #v1_day can be treated as truncation day
names(cr2)[5] <- "trunc_day";
cr2$X_hh <- 1*(mdata$intervention==3|mdata$intervention==4)
cr2$X_fm <- 1*(mdata$intervention==4)
cr2 <- merge(cr2,tmp[c("hhID","member","ILI")],by=c("hhID","member"),all.x=T)
cr2 <- cr2[!is.na(cr2$ILI)&cr2$hhID<665,]
cr2 <- merge(cr2,demog[c("hhID","member","age")],by=c("hhID","member"),all.x=T)

crdata <- rbind(cr,cr2)

# impute adult/child according to relationship
crdata$age[crdata$hhID==148&crdata$member<=2] <- 40; crdata$age[crdata$hhID==148&crdata$member==3] <- 10
crdata$age[crdata$hhID==273] <- 40;  crdata$age[crdata$hhID==283&crdata$member==3] <- 40;
crdata$age[crdata$hhID==301&crdata$member==1] <- 40;
crdata$age[crdata$hhID==317&crdata$member==4] <- 10; crdata$age[crdata$hhID==317&crdata$member!=4] <- 40
crdata$age[crdata$hhID==398&crdata$member<=2] <- 40
crdata$age[crdata$hhID==419&crdata$member==1] <- 40; crdata$age[crdata$hhID==435&crdata$member==1] <- 40
crdata$age[crdata$hhID==444&crdata$member==1] <- 40; crdata$age[crdata$hhID==466&crdata$member==1] <- 40
crdata$age[crdata$hhID==549&crdata$member==1] <- 40; crdata$age[crdata$hhID==603&crdata$member==1] <- 40
crdata$age[crdata$hhID==626&crdata$member==1] <- 40; crdata$age[crdata$hhID==645&(crdata$member==2|crdata$member==5)] <- 40
crdata$age[crdata$hhID==647&crdata$member==1] <- 40

crdata.adult <- crdata[crdata$age>=16,]
crdata.child <- crdata[crdata$age<16,]

crdata2 <- crdata
crdata2$child <- 1*(crdata2$age<16)

## consider cluster effect
ninfect <- table(crdata$hhID[crdata$delta==1])
ninfect2 <- data.frame(hhID=as.numeric(names(ninfect)),ninfect=as.numeric(ninfect))
crdata3 <- merge(crdata2,ninfect2,all.x=T); crdata3$ninfect[is.na(crdata3$ninfect)] <- 0
crdata3 <- crdata3[order(crdata3$hhID,crdata3$T0,-crdata3$delta,crdata3$member),]

hkdata <- crdata3
# mark the possible tertiary cases
hkdata$mark2 <- hkdata$mark1 <- 0; hkdata$T2 <- hkdata$T1 <- hkdata$T0 # time to infection from index/first case/second case
id <- unique(hkdata[hkdata$ninfect>=2,c("hhID","ninfect")])
for(i in 1:nrow(id)){
  rownum <- which(hkdata$hhID==id$hhID[i])[1]
  if(id$ninfect[i]==2){
    if(hkdata$T0[rownum+1]-hkdata$T0[rownum]>0){
      hkdata$mark1[rownum+1] <- 1
      hkdata$T1[rownum+1] <- hkdata$T0[rownum+1]-hkdata$T0[rownum]
    }
  }
  else if(id$ninfect[i]==3){
    # for the second infected case
    if(hkdata$T0[rownum+1]-hkdata$T0[rownum]>0){
      hkdata$mark1[rownum+1] <- 1
      hkdata$T1[rownum+1] <- hkdata$T0[rownum+1]-hkdata$T0[rownum]
    }
    # for the third infected case
    if(hkdata$T0[rownum+2]-hkdata$T0[rownum]>0){
      hkdata$mark1[rownum+2] <- 1
      hkdata$T1[rownum+2] <- hkdata$T0[rownum+2]-hkdata$T0[rownum]
    }
    if(hkdata$T0[rownum+2]-hkdata$T0[rownum+1]>0){
      hkdata$mark2[rownum+2] <- 1
      hkdata$T2[rownum+2] <- hkdata$T0[rownum+2]-hkdata$T0[rownum+1]
    }
  }
}

#
# End of script.
#
