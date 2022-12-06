#
# R syntax to reproduce information for Table 2 from:
#
# Cowling BJ, Fang VJ, Nishiura H, et al.
# Increased risk of noninfluenza respiratory virus infections
# sssociated with receipt of inactivated influenza vaccine
# CID, 2012.
#
# Last updated by Fang VJ and Cowling BJ.
# January 24, 2013
#

require(chron)
dir <- "http://sph.hku.hk/data/KiddivaxPilotV2/"
sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep="")); symp <- symp[symp$member==0,]
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))

sero <- merge(sero[sero$member==0&sero$hhID!=9200&sero$hhID!=9219,],random,all.x=T)
sero$miday <- as.numeric(dates(as.character(sero$date.mids),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$miday[is.na(sero$miday)] <- as.numeric(dates("15/4/2009",format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$day0 <- as.numeric(dates(as.character(sero$start),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$day1 <- as.numeric(dates(as.character(sero$end),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$wlen <- sero$miday-sero$day0; sero$slen <- sero$day1-sero$miday  # length of study period
sero$len <- sero$day1-sero$day0

#
symp$day <- dates(as.character(symp$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")
symp <- merge(symp,sero[c("hhID","miday")],all.x=T)
symp$ARI <- 1*((1*(symp$bodytemp>=37.8)+symp$headache+symp$cough+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm)>=2)
symp$ILI <- 1*(symp$bodytemp>=37.8&(symp$cough==1|symp$sthroat==1))
symp <- merge(symp,sero[c("hhID","day0","day1")],all.x=T); symp <- symp[symp$day>=symp$day0&symp$day<=symp$day1,1:19]
symp <- symp[order(symp$hhID,symp$member,symp$day),]
missSD <- c(9134,9135,9138,9156,9166,9170,9172,9180,9200,9205,9215,9219)

cuminc <- function(epi,len){
  temp <- unique(epi[c("hhID","member","intervention")])
  temp$subject <- 1:nrow(temp)
  epi <- merge(epi,temp[-3],by=c("hhID","member"),all.x=T)
  
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
             epi$epi[init.row+j] <- epi$epi[init.row+j-1]+1             # round(sum(da$len[da$intervention=="TIV"]/365))
             epi.row <- init.row+j
          }
          j <- j+1
     }
     temp$nepi[i] <- epi$epi[init.row+j-1]
  }

   alld <- merge(len,temp[c("hhID","nepi")],all.x=T); alld$nepi[is.na(alld$nepi)] <- 0
   n1 <- sum(alld$nepi[alld$intervention=="TIV"]); n2 <- sum(alld$nepi[alld$intervention=="placebo"])
   N1 <- round(sum(len$len[len$intervention=="TIV"]/365)); N2 <- round(sum(len$len[len$intervention=="placebo"]/365)) # per person-year
   
   alld$group <- 1*(alld$intervention=="placebo")
   m11 <- glm(nepi ~ intervention+offset(log(len)),family="poisson",data=alld); m10 <- glm(nepi ~ group+offset(log(len)),family="poisson",data=alld)
   m12 <- glm(nepi ~ offset(log(len)),family="poisson",data=alld)
   dev.dff <- anova(m12, m11); p.value <- 1-pchisq(dev.dff[2,4], dev.dff[2,3]) # dev.dff[2,8]
   p1 <- exp(m10$coef[1]+log(365))*1000; 
   p1.CI <- exp(c(m10$coef[1]-1.96*sqrt(diag(vcov(m10)))[1]+log(365),m10$coef[1]+1.96*sqrt(diag(vcov(m10)))[1]+log(365)))*1000
   p2 <- exp(m11$coef[1]+log(365))*1000; 
   p2.CI <- exp(c(m11$coef[1]-1.96*sqrt(diag(vcov(m11)))[1]+log(365),m11$coef[1]+1.96*sqrt(diag(vcov(m11)))[1]+log(365)))*1000
   RR <- exp(m11$coef[2]); RR.CI <- exp(c(m11$coef[2]-1.96*sqrt(diag(vcov(m11)))[2],m11$coef[2]+1.96*sqrt(diag(vcov(m11)))[2]))
   
   output <- as.numeric(c(round(c(p1,p1.CI,p2,p2.CI),0),round(c(RR,RR.CI,p.value),2)))
   output  
}

tab <- matrix(NA,ncol=10,nrow=4)
colnames(tab) <- c("TIV","CI_low","CI_up","Placebo","CI_low","CI_up","RR","CI_low","CI_up","p-value"); rownames(tab) <- c("Womter ARI","Winter FARI","Summer ARI","Sumemr FARI")

# winter ARI
wARI <- merge(symp[symp$ARI==1&symp$day<=symp$miday,c("hhID","member","day")],sero[c("hhID","intervention")],all.x=T)
wlen <- sero[!(sero$hhID%in%missSD),c("hhID","intervention","wlen")]; names(wlen)[3] <- "len"
tab[1,] <- cuminc(wARI,wlen)

# winter FARI
wFARI <- merge(symp[symp$ILI==1&symp$day<=symp$miday,c("hhID","member","day")],sero[c("hhID","intervention")],all.x=T)
wlen <- sero[!(sero$hhID%in%missSD),c("hhID","intervention","wlen")]; names(wlen)[3] <- "len"
tab[2,] <- cuminc(wFARI,wlen)

# sumemr ARI
sARI <- merge(symp[symp$ARI==1&symp$day>symp$miday,c("hhID","member","day")],sero[c("hhID","intervention")],all.x=T)
slen <- sero[!(sero$hhID%in%missSD),c("hhID","intervention","slen")]; names(slen)[3] <- "len"
tab[3,] <- cuminc(sARI,slen)

# summer FARI
sFARI <- merge(symp[symp$ILI==1&symp$day>symp$miday,c("hhID","member","day")],sero[c("hhID","intervention")],all.x=T)
slen <- sero[!(sero$hhID%in%missSD),c("hhID","intervention","slen")]; names(slen)[3] <- "len"
tab[4,] <- cuminc(sFARI,slen)

tab

# End of script


