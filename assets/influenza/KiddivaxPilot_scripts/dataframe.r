#
# R syntax to reformat raw data for:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 11, 2010

dir <- "http://sph.hku.hk/data/KiddivaxPilotV1/"
require(chron)

###Reading raw data files
sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))
symp <- read.csv(paste(dir, "symptom_d.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
arr <- read.csv(paste(dir, "ARR_h.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))

##Flu infection by swab---------------------------------------------------------------------------------------------

swab$swab.summer<- 1*(dates(as.character(swab$sub.date),format="d/m/y")-dates("05/01/2009",format="m/d/y")>0)
swab.sh1w<-unique(swab[swab$sub.sh1==1 & swab$swab.summer==0,c("hhID","member")])
swab.sh1s<-unique(swab[swab$sub.sh1==1 & swab$swab.summer==1,c("hhID","member")])
swab.sh3w<-unique(swab[swab$sub.sh3==1 & swab$swab.summer==0,c("hhID","member")])
swab.sh3s<-unique(swab[swab$sub.sh3==1 & swab$swab.summer==1,c("hhID","member")])
swab.bw<-unique(swab[swab$sub.b==1 & swab$swab.summer==0,c("hhID","member")])
swab.bs<-unique(swab[swab$sub.b==1 & swab$swab.summer==1,c("hhID","member")])
swab.ph1<-unique(swab[swab$sub.ph1==1 ,c("hhID","member")])

swab.sh1w$swab.sh1w <-1
swab.sh3w$swab.sh3w <-1
swab.ph1$swab.ph1 <-1
swab.bw$swab.bw <-1

##Flu infection by serology------------------------------------------------------------------------------------------
#convert ph1 titer to numeric & NA for Positive PositiveOne & 5 for Negative
sero$ph1.mid.titer[sero$ph1.mid.titer=="Negative"] <-5
sero$ph1.post.titer[sero$ph1.post.titer=="Negative" | sero$ph1.post.titer=="PositiveOne"] <-5

#pandemic A/H1N1 infection
sero$ph1s <- 1*(as.numeric(as.character(sero$ph1.post.titer))/as.numeric(as.character(sero$ph1.mid.titer))>=4)
sero$ph1s[is.na(sero$ph1s) & sero$ph1.post.titer=="Positive"] <- 1

sero$ph1.mid.titer <-as.numeric(as.character(sero$ph1.mid.titer))
sero$ph1.post.titer <-as.numeric(as.character(sero$ph1.post.titer))

#Flu B infection by serology -----------------------------------------------------------------------------------
# B Florida
sero$flo.w <-NA
sero$flo.w[sero$member==0] <-1*(sero$FluB.Florida.mids[sero$member==0]/sero$FluB.Florida.postv[sero$member==0]>=4)
sero$flo.w[sero$member!=0] <-1*(sero$FluB.Florida.mids[sero$member!=0]/sero$FluB.Florida.pre[sero$member!=0]>=4)
sero$flo.s <-NA
sero$flo.s <-1*(sero$FluB.Florida.posts/sero$FluB.Florida.mids>=4)
sero$flo <-1* (sero$flo.w==1 | sero$flo.s==1)

#Seasonal H1N1 winter & summer -------------------------------------------------------------------------------------
sero$sh1w[sero$member==0] <- 1*(sero$sh1.mids[sero$member==0]/sero$sh1.postv[sero$member==0] >=4)
sero$sh1w[sero$member>0] <- 1*(sero$sh1.mids[sero$member>0]/sero$sh1.pre[sero$member>0] >=4)
sero$sh1s <- 1*(sero$sh1.posts/sero$sh1.mids>=4)

#Seasonal H3N2 winter & summer
sero$sh3w[sero$member==0] <- 1*(sero$sh3.mids[sero$member==0]/sero$sh3.postv[sero$member==0] >=4)
sero$sh3w[sero$member>0] <- 1*(sero$sh3.mids[sero$member>0]/sero$sh3.pre[sero$member>0] >=4)
sero$sh3s <- 1*(sero$sh3.posts/sero$sh3.mids>=4)

#Edit ILI & ARI data---------------------------------------------------------------------------------------------------

#convert btemp to fever or not
symp$fever <- 1*(symp$bodytemp>=37.8)
symp[is.na(symp)]<-0
symp$total <- symp$sthroat + symp$cough + symp$phlegm + symp$rnose + symp$pmuscle + symp$fever + symp$headache
symp$total.cdc <- symp$fever + 1*(symp$cough + symp$sthroat>=1)

##Dividing symp into before and after may
symp$summer.d <-dates(as.character(symp$date),format="d/m/y")-dates("01/05/2009",format="d/m/y")
symp$summer <- 1*(symp$summer.d>=0)
symp.summer<-symp[symp$summer==1,]
symp.winter<-symp[symp$summer==0,]

###ILI (2/7 symptoms) List
iliw.list<-unique(symp.winter[symp.winter$total.cdc>1,c("hhID","member")])
ilis.list<-unique(symp.summer[symp.summer$total.cdc>1,c("hhID","member")])
iliw.list$iliw<-1
ilis.list$ilis<-1

##ili in august sept
pday <-"1/9/2009"  #change for diff cutoff
ilip.list<-unique(symp[symp$total.cdc>1 & dates(as.character(symp$date),format="d/m/y")-dates(pday,format="d/m/y")>0, c("hhID","member")])
ilip.list$ilip<-1

####ARI fever + cough/sore throat List
ariw.list<-unique(symp.winter[symp.winter$total>1,c("hhID","member")])
aris.list<-unique(symp.summer[symp.summer$total>1,c("hhID","member")])
ariw.list$ariw<-1
aris.list$aris<-1

### Index.pre---------------------------------------------------------------------------------------------------------
index.pre.raw <- demog[demog$member==0,]
index.pre <-merge(index.pre.raw,random,by="hhID",all.x=T) #Take the 1st entry data after correction
index.pre$age[index.pre$age<6] <- 6

##merge in serology results
index.pre <-merge(index.pre,sero[sero$member==0,],by=c("hhID","member"),all.x=T)
index.pre$flo.s[is.na(index.pre$flo.s)] <-0

## impute chronic illnesses 9167,9195=asthma, 9163=asthma, chronic respiratory, lung or breathing problem (not asthma), 9217=cardiac disease
index.pre$chron <-0
index.pre$chron[index.pre$hhID==9167|index.pre$hhID==9163|index.pre$hhID==9195|index.pre$hhID==9217] <-1

##convert vac08 NA into 0
index.pre$vac08[is.na(index.pre$vac08)] <- 0

##create vac 09==intervention
index.pre$vac09[index.pre$intervention=="TIV"]<-1
index.pre$vac09[index.pre$intervention=="placebo"]<-0

##merge in ili
index.pre<-merge(index.pre,iliw.list, by=c("hhID","member"),all.x=T)
index.pre<-merge(index.pre,ilis.list, by=c("hhID","member"),all.x=T)
index.pre<-merge(index.pre,ilip.list,by=c("hhID","member"),all.x=T)

index.pre$iliw[is.na(index.pre$iliw)]<-0
index.pre$ilis[is.na(index.pre$ilis)]<-0
index.pre$ilip[is.na(index.pre$ilip)]<-0

##merge in ARI
index.pre<-merge(index.pre,ariw.list, by=c("hhID","member"),all.x=T)
index.pre<-merge(index.pre,aris.list, by=c("hhID","member"),all.x=T)
index.pre$ariw[is.na(index.pre$ariw)]<-0
index.pre$aris[is.na(index.pre$aris)]<-0

#make ilip NA if posts date is after ilip cutline
index.pre$ilip[dates(as.character(index.pre$date.posts),format="d/m/y")<dates(pday,format="d/m/y")] <- NA

#merge swab A and B +ve
index.pre<-merge(index.pre,swab.sh1w,by=c("hhID","member"),all.x=T)
index.pre<-merge(index.pre,swab.sh3w,by=c("hhID","member"),all.x=T)
index.pre<-merge(index.pre,swab.ph1,by=c("hhID","member"),all.x=T)
index.pre<-merge(index.pre,swab.bw,by=c("hhID","member"),all.x=T)

index.pre$swab.sh1w[is.na(index.pre$swab.sh1w)] <-0
index.pre$swab.sh3w[is.na(index.pre$swab.sh3w)] <-0
index.pre$swab.ph1[is.na(index.pre$swab.ph1)] <-0
index.pre$swab.bw[is.na(index.pre$swab.bw)] <-0
index.pre$swab.bs <-0
index.pre$swab.sh1s <-0
index.pre$swab.sh3s <-0

#number of days from 31/8 posts
index.pre$day.p <- dates(as.character(index.pre$date.posts),format="d/m/y")-dates("31/08/09",format="d/m/y")

# make ILI/ARI/miss-sch NA for those hhs without symptom diary returned
missSD <- c(9134,9135,9138,9156,9166,9170,9172,9180,9200,9205,9215,9219)
index.pre$iliw[index.pre$hhID%in%missSD] <- index.pre$ariw[index.pre$hhID%in%missSD] <- 
index.pre$ilis[index.pre$hhID%in%missSD] <- index.pre$aris[index.pre$hhID%in%missSD] <- 
index.pre$ilip[index.pre$hhID%in%missSD] <-  NA


### member.pre---------------------------------------------------------------------------------------------
member.pre <- demog[demog$member>0,]
member.pre <-merge(member.pre,random,by="hhID",all.x=T) #Take the 1st entry data after correction
member.pre$age[member.pre$hhID==9158 & member.pre$member==4] <- 60
member.pre$age[member.pre$hhID==9193 & member.pre$member==4] <- 60
member.pre$age[member.pre$hhID==9111 & member.pre$member==4] <- 1

##merge in serology results and exclude cases without complete blood
member.pre <-merge(member.pre,sero[sero$member!=0,],by=c("hhID","member"),all.x=T)
member.pre$vac09 <- member.pre$vac09.pre
member.pre$vac09[is.na(member.pre$vac09)] <- 0

#Revert sero flu for vaccine reaction
###############################################

member.pre$sh1w[member.pre$vac09.mid==1 & member.pre$vac09!=1 & member.pre$swab.sh1w!=1] <-0
member.pre$sh3w[member.pre$vac09.mid==1 & member.pre$vac09!=1 & member.pre$swab.sh3w!=1] <-0
member.pre$flo.w[member.pre$vac09.mid==1 & member.pre$vac09!=1] <-0

member.pre$sh1s[member.pre$vac09.post==1 & member.pre$vac09.mid!=1 & member.pre$vac09!=1] <-0
member.pre$sh3s[member.pre$vac09.post==1 & member.pre$vac09.mid!=1 & member.pre$vac09!=1] <-0
member.pre$flo.s[member.pre$vac09.post==1 & member.pre$vac09.mid!=1 & member.pre$vac09!=1] <-0

#vac08
member.pre$vac08[is.na(member.pre$vac08)] <- 0

##merge in ili
member.pre<-merge(member.pre,iliw.list, by=c("hhID","member"),all.x=T)
member.pre<-merge(member.pre,ilis.list, by=c("hhID","member"),all.x=T)
member.pre<-merge(member.pre,ilip.list,by=c("hhID","member"),all.x=T)

member.pre$iliw[is.na(member.pre$iliw)]<-0
member.pre$ilis[is.na(member.pre$ilis)]<-0
member.pre$ilip[is.na(member.pre$ilip)]<-0

##merge in ARI
member.pre<-merge(member.pre,ariw.list, by=c("hhID","member"),all.x=T)
member.pre<-merge(member.pre,aris.list, by=c("hhID","member"),all.x=T)
member.pre$ariw[is.na(member.pre$ariw)]<-0
member.pre$aris[is.na(member.pre$aris)]<-0

#merge swab A and B +ve
member.pre<-merge(member.pre,swab.sh1w,by=c("hhID","member"),all.x=T)
member.pre<-merge(member.pre,swab.sh3w,by=c("hhID","member"),all.x=T)
member.pre<-merge(member.pre,swab.ph1,by=c("hhID","member"),all.x=T)

member.pre$swab.sh1w[is.na(member.pre$swab.sh1w)] <-0
member.pre$swab.sh3w[is.na(member.pre$swab.sh3w)] <-0
member.pre$swab.ph1[is.na(member.pre$swab.ph1)] <-0
member.pre$swab.sh1s <-0
member.pre$swab.sh3s <-0
member.pre$swab.bw <-0
member.pre$swab.bs <-0

#make ilip NA if posts date is after ilip cutline
member.pre$ilip[dates(as.character(member.pre$date.posts),format="d/m/y")<dates(pday,format="d/m/y")] <- NA

#no of days from 31/8 for posts
member.pre$day.p <- dates(as.character(member.pre$date.posts),format="d/m/y")-dates("31/08/09",format="d/m/y")

# make ILI/ARI/miss-sch NA for those hhs without symptom diary returned
missSD <- c(9134,9135,9138,9156,9166,9170,9172,9180,9200,9205,9215,9219)
member.pre$iliw[member.pre$hhID%in%missSD] <- member.pre$ariw[member.pre$hhID%in%missSD] <-
member.pre$ilis[member.pre$hhID%in%missSD] <- member.pre$aris[member.pre$hhID%in%missSD] <-
member.pre$ilip[member.pre$hhID%in%missSD] <- NA

member.pre$flo <-1* (member.pre$flo.w==1 | member.pre$flo.s==1)

#merge index.pre and member.pre to all for plot
all <-rbind(index.pre[,names(index.pre) %in% names(member.pre)], member.pre[,names(member.pre) %in% names(index.pre)])
all$agegp <-cut(all$age,c(0,5,12,17,30,50,65,100))


## AR dataframe--------------------------------------------------------------------------------------------------------
arr <- data.frame(lapply(arr,function(x,...){x[is.na(x)] <- 0 ; x}))

arr$pain <- arr$bruise <- arr$red <- arr$swell <-
            arr$mpain <- arr$cough <- arr$headache <- arr$tired <- arr$chill <-arr$fever <-  NA
for (i in 1:nrow(arr)){
  for (j in 46:55){
   arr[i,j] <- max(arr[i,(j-46)*4+6:9],na.rm=TRUE)
  }
}

#merge intervention group to arr
arr2 <- merge(arr, index.pre[,c("hhID","intervention")] , by="hhID",all.y=T)
arr2[is.na(arr2)]<-0


#
# End of script.
#

