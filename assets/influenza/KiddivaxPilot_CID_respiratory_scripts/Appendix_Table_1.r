#
# R syntax to reproduce information for Appendix Table 1 from:
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
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep="")); swab <- swab[swab$member==0,]
resp <- read.csv(paste(dir, "resplex.csv", sep=""))

#
sero <- merge(sero[sero$member==0&sero$hhID!=9200&sero$hhID!=9219,],random,all.x=T)
sero$day0 <- as.numeric(dates(as.character(sero$start),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$day1 <- as.numeric(dates(as.character(sero$end),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$len <- sero$day1-sero$day0

sero$ph1.mid.titer <- as.numeric(as.character(sero$ph1.mid.titer)); sero$ph1.post.titer <- as.numeric(as.character(sero$ph1.post.titer))
sero$ph1s <- 1*(sero$ph1.post.titer/sero$ph1.mid.titer>=4)
sero$sh1w <- 1*(sero$sh1.mids/sero$sh1.postv>=4); sero$sh1s <- 1*(sero$sh1.posts/sero$sh1.mids>=4)
sero$sh3w <- 1*(sero$sh3.mids/sero$sh3.postv>=4); sero$sh3s <- 1*(sero$sh3.posts/sero$sh3.mids>=4)
sero$b.w <- 1*(sero$FluB.Florida.mids/sero$FluB.Florida.postv>=4); sero$b.s <- 1*(sero$FluB.Florida.posts/sero$FluB.Florida.mids>=4)
sero$b.w[sero$hhID==9120|sero$hhID==9128|sero$hhID==9137|sero$hhID==9160] <- 1  # titer against B.Brisbane (not include in online raw data)
sero$swab.sh1 <- 1*(sero$hhID%in%swab$hhID[swab$sub.sh1==1]); sero$swab.sh3 <- 1*(sero$hhID%in%swab$hhID[swab$sub.sh3==1])
sero$swab.ph1 <- 1*(sero$hhID%in%swab$hhID[swab$sub.ph1==1]); sero$swab.b <- 1*(sero$hhID%in%swab$hhID[swab$sub.b==1])

index <- sero[c("hhID","member","sh1w","sh1s","sh3w","sh3s","b.w","b.s","ph1s","swab.sh1","swab.sh3","swab.b","swab.ph1")]
index <- data.frame(lapply(index,function(x,...){x[is.na(x)] <- 0 ; x}))
#
index$sh1s[index$hhID==9101] <- 0; index$b.w[index$hhID==9137] <- 0; index$sh3s[index$hhID==9193] <- 0;
index$sh1w[index$hhID==9199] <- 0; index$sh1s[index$hhID==9211] <- 0; index$ph1s[index$hhID==9212] <- 0
#
index$sH1N1 <- pmax(1*(index$sh1w==1)+1*(index$sh1s==1),1*(index$swab.sh1==1))
index$sH3N2 <- pmax(1*(index$sh3w==1)+1*(index$sh3s==1),1*(index$swab.sh3==1))
index$sB <-  pmax(1*(index$b.w==1)+1*(index$b.s==1),1*(index$swab.b==1))
index$pH1N1 <- 1*(index$ph1s==1|index$swab.ph1==1)

labresult <- index[c("hhID","member","sH1N1","sH3N2","sB","pH1N1")]
resp <- resp[resp$Resplex.II.plus!="Neg"&resp$Resplex.II.plus!="Flu A"&resp$Resplex.II.plus!="Flu B"&resp$Resplex.II.plus!="Swine H1N1"&resp$Resplex.II.plus!="",]
resp$day <- dates(as.character(resp$date),format="d/m/y")-dates("1/1/2009",format="d/m/y")
resp <- resp[order(resp$hhID,resp$day),]
for(i in 1:4){
  if(i==1) temp <- resp
  if(i==2) temp <- resp[resp$Resplex.II.plus=="Rhino",]
  if(i==3) temp <- resp[resp$Resplex.II.plus=="Coxsackie/Echo",]
  if(i==4) temp <- resp[resp$Resplex.II.plus!="Rhino"&resp$Resplex.II.plus!="Coxsackie/Echo",]
  temp$mark <- 0
  for(j in 2:nrow(temp)){
    if(temp$hhID[j]==temp$hhID[j-1]&temp$day[j]-temp$day[j-1]<=7) temp$mark[j] <- 1
  }
  temp <- temp[temp$mark==0,]
  id <- data.frame(hhID=names(table(temp$hhID)), resp=as.numeric(table(temp$hhID))); names(id)[2] <- paste("resp",i,sep="")
  labresult <- merge(labresult,id,by="hhID",all.x=T)
}
labresult$resp1[is.na(labresult$resp1)] <- labresult$resp2[is.na(labresult$resp2)] <-
labresult$resp3[is.na(labresult$resp3)] <- labresult$resp4[is.na(labresult$resp4)] <- 0

labresult <- merge(labresult,random,all.x=T)
labresult <- merge(labresult,sero[c("hhID","len")],all.x=T)

#
cuminc2 <- function(ldata){
   n <- c(sum(ldata$event[ldata$intervention=="TIV"]),sum(ldata$event[ldata$intervention=="placebo"]))
   N <- rev(table(ldata$intervention))

   ldata$group <- 1*(ldata$intervention=="placebo")
   m11 <- glm(event ~ intervention+offset(log(len)),family="poisson",data=ldata); m10 <- glm(event ~ group+offset(log(len)),family="poisson",data=ldata)
   m12 <- glm(event ~ offset(log(len)),family="poisson",data=ldata)
   dev.dff <- anova(m12, m11); p.value <- 1-pchisq(dev.dff[2,4], dev.dff[2,3]) # dev.dff[2,8]
   p1 <- exp(m10$coef[1]+log(365))*1000;
   p1.CI <- exp(c(m10$coef[1]-1.96*sqrt(diag(vcov(m10)))[1]+log(365),m10$coef[1]+1.96*sqrt(diag(vcov(m10)))[1]+log(365)))*1000
   p2 <- exp(m11$coef[1]+log(365))*1000;
   p2.CI <- exp(c(m11$coef[1]-1.96*sqrt(diag(vcov(m11)))[1]+log(365),m11$coef[1]+1.96*sqrt(diag(vcov(m11)))[1]+log(365)))*1000
   output <- c(n[1],round(c(p1,p1.CI)),n[2],round(c(p2,p2.CI)),round(p.value,2))
   output
}

tab <- matrix(NA,ncol=9,nrow=9)
colnames(tab) <- c("TIV","rate","CI_low","CI_up","Placebo","rate","CI_low","CI_up","p-value"); 
rownames(tab) <- c("Any flu","sH1","sH3","B","pH1","Any non-flu","Rhino","Coxsackie","other")

# any influenza
temp <- labresult; temp$event <- temp$sH1N1+temp$sH3N2+temp$sB; tab[1,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$sH1N1; tab[2,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$sH3N2; tab[3,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$sB; tab[4,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$pH1N1; tab[5,] <- cuminc2(temp)

# any non-influenza virus
temp <- labresult; temp$event <- temp$resp1; tab[6,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$resp2; tab[7,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$resp3; tab[8,] <- cuminc2(temp)
temp <- labresult; temp$event <- temp$resp4; tab[9,] <- cuminc2(temp)

tab

# End of script

