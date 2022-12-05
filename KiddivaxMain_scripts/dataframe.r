#
# R syntax to reformat raw data for:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

require(chron)
dir <- "http://sph.hku.hk/data/KiddivaxMainV1/"

random <- read.csv(paste(dir, "randomcode.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
arr <- read.csv(paste(dir, "ARR.csv", sep=""))
demog <- read.csv(paste(dir, "demographic.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep=""))
sero <- read.csv(paste(dir, "serology.csv", sep=""))

# Define RT-PCR (swab) confirmed infections
swab <- merge(swab,sero[c(1,2,5,6)],by=c("hhID","member"),all.x=T)
swab$mark <- 1*(as.Date(as.character(swab$date),format="%Y%m%d")-as.Date(as.character(swab$start.date),format="%d/%m/%Y")>=0&
                as.Date(as.character(swab$date),format="%Y%m%d")-as.Date(as.character(swab$end.date),format="%d/%m/%Y")<=0)
swab <- swab[swab$mark==1,1:10]

swab.pH1 <- unique(swab[swab$Swine.H1=="P" ,c("hhID","member")]);
swab.sH3 <- unique(swab[swab$H3=="P",c("hhID","member")]); swab.B <- unique(swab[swab$FluB=="P",c("hhID","member")])
swab.pH1$swab.pH1 <-1; swab.sH3$swab.sH3 <-1; swab.B$swab.B <-1

# Define serology-confirmed infections (>=4 fold rise in antibody titers)
sero$pH1 <- 1*(sero$post.season.pH1/sero$postvax.pH1>=4)
sero$sH3 <- 1*(sero$post.season.sH3/sero$postvax.sH3>=4)
sero$B <- 1*(sero$post.season.B.Brisbane/sero$postvax.B.Brisbane>=4)

# construct the data frame
index.pre <- merge(demog,random,by="hhID",all.x=T) 
index.pre$age[index.pre$age<6] <- 6

# merge in serology results
index.pre <- merge(index.pre,sero,by=c("hhID","member"),all.x=T)

# vac0910.pre=1 if intervention=TIV
index.pre$vac0910.pre[index.pre$intervention=="TIV"]<-1
index.pre$vac0910.pre[index.pre$intervention=="placebo"]<-0

# merge in ARI/ILI
symp <- symp[symp$hhID!=2654,]  # hhID=2654 should be dropout, input error??
symp <- merge(symp,sero[c(1,2,5,6)],by=c("hhID","member"),all.x=T)
symp$mark <- 1*(as.Date(as.character(symp$date),format="%d/%m/%Y")-as.Date(as.character(symp$start.date),format="%d/%m/%Y")>=0&
                as.Date(as.character(symp$date),format="%d/%m/%Y")-as.Date(as.character(symp$end.date),format="%d/%m/%Y")<=0)
symp <- symp[symp$mark==1,1:14]

ARI <- unique(symp[symp$ARI==1,1:2]); ARI$ARI <- 1
ILI <- unique(symp[symp$ILI==1,1:2]); ILI$ILI <- 1
index.pre <- merge(index.pre,ARI, by=c("hhID","member"),all.x=T); index.pre <- merge(index.pre,ILI, by=c("hhID","member"),all.x=T)
index.pre$ARI[is.na(index.pre$ARI)] <- index.pre$ILI[is.na(index.pre$ILI)] <- 0

# merge swab A and B +ve
index.pre <- merge(index.pre,swab.pH1,by=c("hhID","member"),all.x=T)
index.pre <- merge(index.pre,swab.sH3,by=c("hhID","member"),all.x=T)
index.pre <- merge(index.pre,swab.B,by=c("hhID","member"),all.x=T)
index.pre$swab.pH1[is.na(index.pre$swab.pH1)] <- 0
index.pre$swab.sH3[is.na(index.pre$swab.sH3)] <- 0
index.pre$swab.B[is.na(index.pre$swab.B)] <- 0

sero <- sero[-c(5,6)]

# count ARI/ILI episodes
cnames <- c("nARI","nILI")
for(k in 1:2){
  if(k==1) epi <- symp[symp$ARI==1,]
  if(k==2) epi <- symp[symp$ILI==1,]
  
  epi$day <- dates(as.character(epi$date),format="d/m/Y")-dates("1/1/2010",format="d/m/Y")
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
  index.pre <- merge(index.pre,temp[c(1,4)],all.x=T); index.pre$nepi[is.na(index.pre$nepi)] <- 0
  names(index.pre)[ncol(index.pre)] <- cnames[k]
}

#
# End of script.
#