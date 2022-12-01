#
# R syntax to reproduce information for Figure 1 from:
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
symp <- read.csv(paste(dir, "symptom_d.csv", sep="")); symp <- symp[symp$member==0,]
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep="")); swab <- swab[swab$member==0,]
resp <- read.csv(paste(dir, "resplex.csv", sep=""))
qmh <- read.csv(paste(dir, "QMHisolate.csv", sep=""))

symp$day <- dates(as.character(symp$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")
symp$ARI <- 1*((1*(symp$bodytemp>=37.8)+symp$headache+symp$cough+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm)>=2)
symp <- symp[order(symp$hhID,symp$member,symp$day),]

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

# construct the data frame - hhID, member, episode count, resplex result (maybe NA)
# look at any swab result within that episode (between start and end period)
epi2 <- unique(epi[c("hhID","member","epi")])      # nrow(epi2): 135
for(i in 1:nrow(epi2)){
  temp.epi <- epi$day[epi$hhID==epi2$hhID[i]&epi$member==epi2$member[i]&epi$epi==epi2$epi[i]]
  epi2$epi.start[i] <- min(temp.epi)
  epi2$epi.end[i] <- max(temp.epi)
}

# virus detected by RT-PCR/multiplex
resp$day <- dates(as.character(resp$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")
names(swab)[2] <- "date"; swab$day <- dates(as.character(swab$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")

swab$stype <- NA; swab$stype[swab$sub.sh1==1] <- "sH1"; swab$stype[swab$sub.sh3==1] <- "sH3";
swab$stype[swab$sub.ph1==1] <- "pH1"; swab$stype[swab$sub.b==1] <- "B"
rs <- merge(resp[c(1,2,4,7)],swab[c(1,3,11,12)],by=c("hhID","member","day"),all.x=T)
rs$labresult <- as.character(rs$Resplex.II.plus); rs$labresult[(is.na(rs$Resplex.II.plus))|(!is.na(rs$Resplex.II.plus)&rs$Resplex.II.plus=="")] <- "Neg"
rs$labresult[!is.na(rs$stype)] <- rs$stype[!is.na(rs$stype)]; rs$labresult[rs$labresult=="Swine H1N1"] <- "pH1"

a <- resp
a$ARI <- NA
for(i in 1:nrow(a)){
  temp.epi <- epi2[epi2$hhID==a$hhID[i]&epi2$member==a$member[i],]
  if(nrow(temp.epi)==0) next
  for(j in 1:nrow(temp.epi)){
    if(a$day[i]>=temp.epi$epi.start[j]-5&a$day[i]<=temp.epi$epi.end[j]+5) a$ARI[i] <- 1
  }
}
a <- a[a$ARI==1&!is.na(a$ARI),]
resp2 <- merge(rs,a[c("hhID","day")],all.y=T)   # specimens with ARI

### for the plot
respc <- merge(resp2,random); respc <- respc[respc$labresult!="Neg",]
flu.record <- list(NA,NA,NA,NA)
for (k in 1:4){
if(k==1) {flu <- respc[respc$member==0&(respc$labresult=="sH1"|respc$labresult=="sH3"|respc$labresult=="pH1"|respc$labresult=="B"),];
          id <- demog[demog$member==0,1:2]}
if(k==2) {flu <- respc[respc$member==0&(respc$labresult=="Rhino"),]; id <- demog[demog$member==0,1:2]}
if(k==3) {flu <- respc[respc$member==0&(respc$labresult=="Coxsackie/Echo"),]; id <- demog[demog$member==0,1:2]}
if(k==4) {flu <- respc[respc$member==0&!(respc$labresult=="sH1"|respc$labresult=="sH3"|respc$labresult=="pH1"|respc$labresult=="B"|
               respc$labresult=="Rhino"|respc$labresult=="Coxsackie/Echo"),]; id <- demog[demog$member==0,1:2]}

flu$epi <- NA; flu$epi[1] <- 1
for(i in 2:nrow(flu)){
  if(flu$hhID[i]==flu$hhID[i-1]&flu$member[i]==flu$member[i-1]&flu$day[i]-flu$day[i-1]<7) flu$epi[i] <- 0
  else flu$epi[i] <- 1
}

flu.record[[k]] <- flu[flu$epi==1,]
}

#
qmh <- qmh[54:104,]
qmh$day <- NA; qmh$day[1] <- as.numeric(dates("10/1/2009",format="d/m/y")-dates("1/1/2009",format="d/m/y"))
for(i in 2:nrow(qmh)){qmh$day[i] <- qmh$day[i-1]+7}
qmh$ph1n1[qmh$day<=180] <- NA # plot start from July

# barplot of positive samples along with timeline

windows(width=4,height=9)
layout(matrix(1:5,ncol=1,byrow=T))
par(mar=c(3,4,2,1))
maintitle <- c("Influenza virus","Rhinovirus","Coxsackie/Echovirus","Other virus")
for(k in 1:4){
nday <- table(flu.record[[k]]$day); bday <- sort(unique(flu.record[[k]]$day)); gp <- unique(flu.record[[k]][order(flu.record[[k]]$day),c("day","intervention")])[,2]

plot(NA,xlim=c(0,350),ylim=c(0,2),axes=F,main=maintitle[k],xlab="",ylab="Number of detections")
for(i in 1:length(bday)){
  lines(rep(bday[i],2),c(0,nday[i]),lty=2-1*(gp[i]=="TIV"),col=1*(gp[i]=="TIV")+1)
}
axis(1,at=cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31)),labels=c("Jan","","","Apr","","","Jul","","","Oct","","","Jan"))
axis(2,at=0:2,las=1)
}

# QMH surveillance plot
plot(NA,xlim=c(0,350),ylim=c(0,0.4),main="Influenza surveillance data",xlab="",ylab="Proportion of positive detections",axes=F)
lines(qmh$day-3.5,qmh$ph1n1/qmh$Total.NPA,lty=1)
lines(qmh$day-3.5,qmh$Seasonal.H1/qmh$Total.NPA,lty=2)
lines(qmh$day-3.5,qmh$Seasonal.H3/qmh$Total.NPA,lty=3)
axis(1,at=cumsum(c(0,31,28,31,30,31,30,31,31,30,31,30,31)),labels=c("Jan","","","Apr","","","Jul","","","Oct","","","Jan"))
axis(2,at=0:4/10,labels=c("0%","10%","20%","30%","40%"),las=1)
legend(0,0.4,legend=c("Seasonal influenza A(H1N1)","Seasonal influenza A(H3N2)","Pandemic influenza A(H1N1)"),
       lty=c(2,3,1),bty="n")

# End of script





