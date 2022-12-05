#
# R syntax to reproduce information for Table 3 from:
#
# Cowling BJ, Fang VJ, Nishiura H, et al.
# Increased risk of noninfluenza respiratory virus infections
# sssociated with receipt of inactivated influenza vaccine
# CID, 2012.
#
# Last updated by Fang VJ and Cowling BJ.
# January 28, 2013
#

require(chron)
dir <- "http://sph.hku.hk/data/KiddivaxPilotV2/"
sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep="")); symp <- symp[symp$member==0,]
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
hchar <- read.csv(paste(dir, "housechar_h.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep="")); swab <- swab[swab$member==0,]
resp <- read.csv(paste(dir, "resplex.csv", sep=""))


sero <- merge(sero[sero$member==0&sero$hhID!=9200&sero$hhID!=9219,],random,all.x=T)
sero$day0 <- as.numeric(dates(as.character(sero$start),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$day1 <- as.numeric(dates(as.character(sero$end),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$len <- sero$day1-sero$day0

#
symp$day <- dates(as.character(symp$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")
symp$ARI <- 1*((1*(symp$bodytemp>=37.8)+symp$headache+symp$cough+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm)>=2)
symp$ILI <- 1*(symp$bodytemp>=37.8&(symp$cough==1|symp$sthroat==1))
symp <- merge(symp,sero[c("hhID","day0","day1")],all.x=T); symp <- symp[symp$day>=symp$day0&symp$day<=symp$day1,1:18]
symp <- symp[order(symp$hhID,symp$member,symp$day),]
missSD <- c(9134,9135,9138,9156,9166,9170,9172,9180,9200,9205,9215,9219)

# ARI
epi <- symp[symp$ARI==1,]
temp <- unique(epi[1:2])
temp$subject <- 1:nrow(temp)
epi <- merge(epi,temp,by=c("hhID","member"),all.x=T)

for (i in 1:nrow(temp)){
   init.row <- nrow(epi[epi$subject<i,])+1
   epi$epi[init.row] <- 1
   j <- 1
   epi.row <- init.row
   while (j <= nrow(epi[epi$subject==i,])-1){
        if(epi$day[init.row+j]-epi$day[init.row+j-1]<7){
           epi$epi[init.row+j] <- epi$epi[init.row+j-1]
        }
        else {
           epi$epi[init.row+j] <- epi$epi[init.row+j-1]+1
           epi.row <- init.row+j
        }
        j <- j+1
   }
   temp$nepi[i] <- epi$epi[init.row+j-1]
}

alld <- merge(sero[c("hhID","member","intervention","len")],temp[c("hhID","nepi")],all.x=T); alld$nepi[is.na(alld$nepi)] <- 0

# construct the data frame - hhID, member, episode count, resplex result (maybe NA)
# look at any swab result within that episode (between start and end period)
epi2 <- unique(epi[c("hhID","member","epi")])      # nrow(epi2): 135
for(i in 1:nrow(epi2)){
  temp.epi <- epi$day[epi$hhID==epi2$hhID[i]&epi$member==epi2$member[i]&epi$epi==epi2$epi[i]]
  epi2$epi.start[i] <- min(temp.epi)
  epi2$epi.end[i] <- max(temp.epi)
}

# swabs collected and labresult
resp$date <- as.Date(resp$date,format="%d/%m/%Y")
names(swab)[2] <- "date"; swab$date <- as.Date(swab$date,format="%d/%m/%Y") 
labresult <- merge(resp[1:4],swab[c(1,2,7:10)],all.x=T)
labresult$labresult <- as.character(labresult$Resplex.II.plus)
labresult$labresult[is.na(labresult$labresult)|(!is.na(labresult$labresult)&labresult$labresult=="")] <- "Neg"
labresult$labresult[labresult$labresult=="Swine H1N1"] <- "pH1"
labresult$labresult[labresult$sub.sh1==1&!is.na(labresult$sub.sh1)] <- "sH1"
labresult$labresult[labresult$sub.sh3==1&!is.na(labresult$sub.sh3)] <- "sH3"
labresult$labresult[labresult$sub.ph1==1&!is.na(labresult$sub.ph1)] <- "pH1"
labresult$day <- dates(as.character(labresult$date),format="y-m-d")-dates("1/1/2009",format="d/m/y")

#
cuminc <- function(epi,lab){
    epi$detect <- 0
    for(i in 1:nrow(epi)){
      temp.lab <- lab[lab$hhID==epi$hhID[i]&lab$day>=epi$epi.start[i]&lab$day<=epi$epi.end[i]+5,]
      if(nrow(temp.lab)==0) next
      else epi$detect[i] <- 1
    }
    alld2 <- alld; alld2$ndetect <- 0
    for(i in 1:nrow(alld2)){
      temp.epi <- epi[epi$hhID==alld2$hhID[i]&epi$detect==1,]
      if(nrow(temp.epi)==0) next
      alld2$ndetect[i] <- nrow(temp.epi)
    }
    #
   alld2$group <- 1*(alld2$intervention=="placebo")
   m11 <- glm(ndetect ~ intervention+offset(log(len)),family="poisson",data=alld2); m10 <- glm(ndetect ~ group+offset(log(len)),family="poisson",data=alld2)
   m12 <- glm(ndetect ~ offset(log(len)),family="poisson",data=alld2)
   dev.dff <- anova(m12, m11); p.value <- 1-pchisq(dev.dff[2,4], dev.dff[2,3]) # dev.dff[2,8]
   p1 <- exp(m10$coef[1]+log(365))*1000;
   p1.CI <- exp(c(m10$coef[1]-1.96*sqrt(diag(vcov(m10)))[1]+log(365),m10$coef[1]+1.96*sqrt(diag(vcov(m10)))[1]+log(365)))*1000
   p2 <- exp(m11$coef[1]+log(365))*1000;
   p2.CI <- exp(c(m11$coef[1]-1.96*sqrt(diag(vcov(m11)))[1]+log(365),m11$coef[1]+1.96*sqrt(diag(vcov(m11)))[1]+log(365)))*1000
   RR <- exp(m11$coef[2]); RR.CI <- exp(c(m11$coef[2]-1.96*sqrt(diag(vcov(m11)))[2],m11$coef[2]+1.96*sqrt(diag(vcov(m11)))[2]))

   if(p1.CI[2]==Inf) p1.CI[2] <- 3/(sum(alld2$len[alld2$intervention=="TIV"])/365)*1000
   if(p2.CI[2]==Inf) p2.CI[2] <- 3/(sum(alld2$len[alld2$intervention=="placebo"])/365)*1000
   output <- as.numeric(c(sum(alld2$ndetect[alld2$intervention=="TIV"]),round(c(p1,p1.CI),0),
                          sum(alld2$ndetect[alld2$intervention=="placebo"]),round(c(p2,p2.CI),0),round(c(p.value),2)))
   output
}

tab <- matrix(NA,ncol=9,nrow=11)
colnames(tab) <- c("TIV","rate","CI_low","CI_up","Placebo","rate","CI_low","CI_up","p-value"); 
rownames(tab) <- c("Any flu","sH1","sH3","B","pH1","Any non-flu","Rhino","Coxsackie","other","ARInovirus","ARInospecimen")

# any seasonal influenza
tab[1,] <- cuminc(epi2,labresult[labresult$labresult=="sH1"|labresult$labresult=="sH3"|labresult$labresult=="Flu B",])
tab[2,] <- cuminc(epi2,labresult[labresult$labresult=="sH1",])
tab[3,] <- cuminc(epi2,labresult[labresult$labresult=="sH3",])
tab[4,] <- cuminc(epi2,labresult[labresult$labresult=="Flu B",])
# pandemic influenza
tab[5,] <- cuminc(epi2,labresult[labresult$labresult=="pH1",])

# any non-influenza virus
tab[6,] <- cuminc(epi2,labresult[labresult$labresult=="Rhino"|labresult$labresult=="Coxsackie/Echo"|
                                 labresult$labresult=="Para 3"|labresult$labresult=="NL63"|labresult$labresult=="OC43"|
                                 labresult$labresult=="HMPV"|labresult$labresult=="229E",])
tab[7,] <- cuminc(epi2,labresult[labresult$labresult=="Rhino",])
tab[8,] <- cuminc(epi2,labresult[labresult$labresult=="Coxsackie/Echo",])
tab[9,] <- cuminc(epi2,labresult[labresult$labresult=="Para 3"|labresult$labresult=="NL63"|labresult$labresult=="OC43"|
                      labresult$labresult=="HMPV"|labresult$labresult=="229E",])

# ARI episode with specimen but no virus detected (1) / ARI episode with no specimen collected (2)
for(j in 1:2){
epi3 <- epi2
epi3$labresult <- NA;
for(i in 1:nrow(epi3)){
  temp.swab <- labresult[labresult$hhID==epi3$hhID[i]&labresult$day>=epi3$epi.start[i]&labresult$day<=epi3$epi.end[i]+5,]
  if(nrow(temp.swab)==0) next
  lab.temp <- unique(temp.swab$labresult)
  if(length(lab.temp)>1) epi3$labresult[i] <- 1
  else epi3$labresult[i] <- 1*(lab.temp!="Neg")
}
#table(epi3$labresult,exclude=NULL)     # 65 episodes with specimens (32 with virus and 33 without virus)
   alld2 <- alld
   alld2$group <- 1*(alld2$intervention=="placebo")
   alld2$epicount <- 0
   for(i in 1:nrow(alld2)){
     if(j==1) alld2$epicount[i] <- nrow(epi3[epi3$hhID==alld2$hhID[i]&epi3$labresult==0&!is.na(epi3$labresult),])  # (1)
     else alld2$epicount[i] <- nrow(epi3[epi3$hhID==alld2$hhID[i]&is.na(epi3$labresult),])  # (2)
   }
   m11 <- glm(epicount ~ intervention+offset(log(len)),family="poisson",data=alld2); m10 <- glm(epicount ~ group+offset(log(len)),family="poisson",data=alld2)
   m12 <- glm(epicount ~ offset(log(len)),family="poisson",data=alld2)
   dev.dff <- anova(m12, m11); p.value <- 1-pchisq(dev.dff[2,4], dev.dff[2,3]) # dev.dff[2,8]
   p1 <- exp(m10$coef[1]+log(365))*1000;
   p1.CI <- exp(c(m10$coef[1]-1.96*sqrt(diag(vcov(m10)))[1]+log(365),m10$coef[1]+1.96*sqrt(diag(vcov(m10)))[1]+log(365)))*1000
   p2 <- exp(m11$coef[1]+log(365))*1000;
   p2.CI <- exp(c(m11$coef[1]-1.96*sqrt(diag(vcov(m11)))[1]+log(365),m11$coef[1]+1.96*sqrt(diag(vcov(m11)))[1]+log(365)))*1000
   tab[j+9,] <- c(sum(alld2$epicount[alld2$intervention=="TIV"]),round(c(p1,p1.CI),0),
                  sum(alld2$epicount[alld2$intervention=="placebo"]),round(c(p2,p2.CI),0),round(p.value,2))
}

tab

# End of script
