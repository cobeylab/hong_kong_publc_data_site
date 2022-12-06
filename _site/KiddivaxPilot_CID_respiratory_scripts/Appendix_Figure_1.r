#
# R syntax to reproduce information for Appendix Figure 1 from:
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

sero <- merge(sero[sero$member==0&sero$hhID!=9200&sero$hhID!=9219,],random,all.x=T)
sero$day0 <- as.numeric(dates(as.character(sero$start),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$day1 <- as.numeric(dates(as.character(sero$end),format="d/m/y")-dates("1/1/2009",format="d/m/y"))

symp$day <- dates(as.character(symp$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")
symp$ARI <- 1*((1*(symp$bodytemp>=37.8)+symp$headache+symp$cough+symp$sthroat+symp$pmuscle+symp$rnose+symp$phlegm)>=2)
symp$ILI <- 1*(symp$bodytemp>=37.8&(symp$cough==1|symp$sthroat==1))
symp <- merge(symp,sero[c("hhID","day0","day1")],all.x=T); symp <- symp[symp$day>=symp$day0&symp$day<=symp$day1,1:19]
symp <- symp[order(symp$hhID,symp$member,symp$day),]

windows(width=7,height=8)
layout(matrix(1:2,ncol=1))
par(mar=c(3,4,2,1))

for(k in 1:2){

# ARI
if(k==1) epi <- symp[symp$ARI==1,]
if(k==2) epi <- symp[symp$ILI==1,]
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

epi2 <- unique(epi[c("hhID","member","epi")])      # nrow(epi2): 135
for(i in 1:nrow(epi2)){
  temp.epi <- epi$day[epi$hhID==epi2$hhID[i]&epi$member==epi2$member[i]&epi$epi==epi2$epi[i]]
  epi2$epi.start[i] <- min(temp.epi)
  epi2$epi.end[i] <- max(temp.epi)
}
epi2 <- merge(epi2,sero[c("hhID","intervention")],all.x=T)

# plot
plotdata <- data.frame(wkcount=1:12,wkstart=0:11*28-18,wkend=0:11*28+9)
for(i in 1:nrow(plotdata)){
  plotdata$nepi.onset.t[i] <- nrow(epi2[epi2$epi.start>=plotdata$wkstart[i]&epi2$epi.start<=plotdata$wkend[i]&epi2$intervention=="TIV",])
  plotdata$nepi.onset.p[i] <- nrow(epi2[epi2$epi.start>=plotdata$wkstart[i]&epi2$epi.start<=plotdata$wkend[i]&epi2$intervention=="placebo",])
  tmp <- sero[sero$day0<=plotdata$wkend[i]&sero$day1>=plotdata$wkstart[i],c("intervention","day0","day1")]
  tmp$len <- pmin(tmp$day1,plotdata$wkend[i])-pmax(tmp$day0,plotdata$wkstart[i])+1
  plotdata$ratet[i] <- plotdata$nepi.onset.t[i]/(sum(tmp$len[tmp$intervention=="TIV"])/365)
  plotdata$ratep[i] <- plotdata$nepi.onset.p[i]/(sum(tmp$len[tmp$intervention=="placebo"])/365)
}

#
#

if(k==1){
plot(NA,xlim=c(0,273),ylim=c(0,4),main="",axes=F,xlab="",ylab="Incidence rate per 1,000 p-y")
lines(plotdata$wkstart,plotdata$ratet);lines(plotdata$wkstart,plotdata$ratep,lty=2)
axis(1,at=c(-31,cumsum(c(0,31,28,31,30,31,30,31,31,30,31))),labels=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
axis(2,at=0:4,labels=0:4*1000,las=1,cex.axis=0.8)
mtext("A",side=3,at=-44,line=0.5,font=2)
}
##
if(k==2){
plot(NA,xlim=c(0,273),ylim=c(0,4),main="",axes=F,xlab="",ylab="Incidence rate per 1,000 p-y")
lines(plotdata$wkstart,plotdata$ratet);lines(plotdata$wkstart,plotdata$ratep,lty=2)
axis(1,at=c(-31,cumsum(c(0,31,28,31,30,31,30,31,31,30,31))),labels=c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))
axis(2,at=0:4,labels=0:4*1000,las=1,cex.axis=0.8)
mtext("B",side=3,at=-44,line=0.5,font=2)
}

}

# End of script


