#
# Script to reproduce information in Appendix Table 8 from:
#
# Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, et al.
# Facemasks and hand hygiene to prevent influenza transmission 
# in households, a randomized trial.
# Annals of Internal Medicine, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Sep 08, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV3/"

# Analyzed households: all 331 households

appt8 <- matrix(rep(NA,190),ncol=10,byrow=FALSE)

hc <- read.csv(paste(dir, "home_pcr.csv", sep=""))
housechar <- read.csv(paste(dir, "hchar_h.csv", sep=""))       
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""))      
symptom <- read.csv(paste(dir, "symptomday_d.csv", sep=""))

### Lab-confirmed secondary cases

mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
hculture <- data.frame(hhID = baseflu$hhID, member = baseflu$member)

hc.temp <- reshape(hc, timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="PCR")
hculture <- merge(hculture,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hculture) <- c("hhID","member","V1","V2","V3")

for (i in 1:nrow(hculture)){
     if (hculture$member[i]!=0 &  !is.na(hculture$V1[i]) & hculture$V1[i]!=0)   hculture$d_contact[i]<-1
     else hculture$d_contact[i]<-0
}

co_index <- data.frame(hhID=unique(hculture$hhID[hculture$d_contact==1]))
co_index$co_index <- 1
hculture <- merge(hculture,co_index,all.x=TRUE)
hculture$co_index[is.na(hculture$co_index)] <- 0
hculture <- hculture[order(hculture$hhID,hculture$member),]

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$d_contact[i]==0 & ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) |
               (hculture$V3[i] !=0 & !is.na(hculture$V3[i]))  ) )
              {hculture$labsedcase[i] <- 1}
	 else {hculture$labsedcase[i] <- 0}
}

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

hclinic$d_contact <- hculture$d_contact

## Define secondary cases
for (i in 1:nrow(hclinic)){
    if ( hclinic$member[i] != 0 & hclinic$d_contact[i]==0 &
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
hculture <- hculture[-c(3:5)]

#### Get arm, index agegp, sex index, vaccine contact, contact sex, contact agegp 

av <- read.csv(paste(dir, "antiviral_m.csv", sep=""))
      
### Add variable 'vaccine' <- 1 if the subject received vaccination in past one year
hculture$vaccine <- baseflu$vaccine08

### Add variable 'indexsex' & 'indexage'
indexage <- baseflu[baseflu$member==0,c(1,4)]
names(indexage)[2] <- "indexage"
indexsex <- baseflu[baseflu$member==0,c(1,3)]
names(indexsex)[2] <- "indexsex"
hculture <- merge(hculture,indexage)
hculture <- merge(hculture,indexsex)

### Add variable 'intervention' and 'indexdelay' (time from symptom onset to intervention)
hculture <- merge(hculture,housechar[,c(1,2,9)])
hculture$arm <- factor(hculture$intervention, levels=c(1,3,4), labels=c("control", "hand", "handmask"))

### Add variable 'age' and 'sex' (for household contact)
hculture$age <- baseflu$age
hculture$sex <- baseflu$male

# add index antiviral use
av <- av[av$member==0,]
av$indexav <- 1
hculture <- merge(hculture,av[c(1,6)],by="hhID",all.x=TRUE)
hculture$indexav[is.na(hculture$indexav)] <- 0
hculture <- hculture[order(hculture$hhID,hculture$member),]

datai <- read.csv(paste(dir, "incomplete_m.csv", sep=""))
datai$co_index <- datai$d_contact <- 0
datai$c2nd2 <- datai$c2nd1 <- datai$labsedcase <- NA
names(datai)[4:5] <- c("sex","vaccine")

datac <- rbind(hculture[c(1,2,14,15,8,11,16,12,3:7)],datai)

### Add variable 'indexsex' & 'indexage'
indexage <- datac[datac$member==0,c(1,3)]
names(indexage)[2] <- "indexage"
indexsex <- datac[datac$member==0,c(1,4)]
names(indexsex)[2] <- "indexsex"
datac <- merge(datac,indexage)
datac <- merge(datac,indexsex)
datac$indexagegp <- 1*(datac$indexage<=15)+1*(datac$indexage<=5)

datac$delay[datac$onset_v1_delay<=36] <- "d1"
datac$delay[datac$onset_v1_delay>36] <- "d2"


# MI

library(Hmisc)

set.seed(102)
data.i <- aregImpute( ~ labsedcase+c2nd1+c2nd2+factor(intervention)+age+sex+vaccine+indexav+co_index, data=datac, n.impute=10)
data.nomiss <- list(datac, datac, datac, datac, datac, datac, datac, datac, datac, datac)

for(i in 1:10){
     data.nomiss[[i]]$labsedcase[is.na(data.nomiss[[i]]$labsedcase)] <- data.i$imputed$labsedcase[,i]
     data.nomiss[[i]]$c2nd1[is.na(data.nomiss[[i]]$c2nd1)] <- data.i$imputed$c2nd1[,i]
     data.nomiss[[i]]$c2nd2[is.na(data.nomiss[[i]]$c2nd2)] <- data.i$imputed$c2nd2[,i]
     data.nomiss[[i]]$sex[is.na(data.nomiss[[i]]$sex)] <- data.i$imputed$sex[,i]
     data.nomiss[[i]]$age[is.na(data.nomiss[[i]]$age)] <- data.i$imputed$age[,i]
     data.nomiss[[i]]$agegp <- 1*(data.nomiss[[i]]$age<=15)+1*(data.nomiss[[i]]$age<=5)
     data.nomiss[[i]] <- data.nomiss[[i]][data.nomiss[[i]]$d_contact==0&data.nomiss[[i]]$member>0,]
}

library(gee)

model1 <- model2 <- model3 <-  list(NA)

for(i in 1:10){
      model1[[i]] <- gee(labsedcase ~ factor(intervention)+factor(agegp)+sex+vaccine+factor(indexagegp)+indexsex+indexav+co_index ,
                         data=data.nomiss[[i]], id=factor(data.nomiss[[i]]$hhID), corstr = "exchangeable", family="binomial")
      model2[[i]] <- gee(c2nd1 ~ factor(intervention)+factor(agegp)+sex+vaccine+factor(indexagegp)+indexsex+indexav+co_index  ,
                         data=data.nomiss[[i]], id=factor(data.nomiss[[i]]$hhID), corstr = "exchangeable", family="binomial")
      model3[[i]] <- gee(c2nd2 ~ factor(intervention)+factor(agegp)+sex+vaccine+factor(indexagegp)+indexsex+indexav+co_index  ,
                         data=data.nomiss[[i]], id=factor(data.nomiss[[i]]$hhID), corstr = "exchangeable", family="binomial")
}

#
# Combine imputed results and summarise parameter estimates
#

combine.mi <- function(model, n.impute=10){
	betas <- matrix(c(model[[1]]$coef, model[[2]]$coef,
		model[[3]]$coef, model[[4]]$coef, model[[5]]$coef,
		model[[6]]$coef, model[[7]]$coef,
		model[[8]]$coef, model[[9]]$coef, model[[10]]$coef), byrow=FALSE, ncol=10)
	vars <- matrix(c(sqrt(diag(model[[1]][[20]])), sqrt(diag(model[[2]][[20]])),
		sqrt(diag(model[[3]][[20]])), sqrt(diag(model[[4]][[20]])), sqrt(diag(model[[5]][[20]])),
		sqrt(diag(model[[6]][[20]])), sqrt(diag(model[[7]][[20]])), sqrt(diag(model[[8]][[20]])),
		sqrt(diag(model[[9]][[20]])), sqrt(diag(model[[10]][[20]]))), byrow=FALSE, ncol=10)^2
	coef.names <- names(model[[1]]$coef)
	mean.coefs <- rowMeans(betas)
	Ubar <- rowMeans(vars)
	B <- rowSums((betas - mean.coefs)*(betas-mean.coefs) /
		(n.impute - 1))
	T <- (1 + 1/n.impute) * B + Ubar
	degf <- (n.impute - 1)*(1 + Ubar / ((1 + 1/n.impute)*B))*
		(1 + Ubar / ((1 + 1/n.impute)*B))
	output <- data.frame(OR = exp(mean.coefs),
		                   lowerCI = exp(mean.coefs - qt(0.975, df=degf)*sqrt(T)),
	                     upperCI = exp(mean.coefs + qt(0.975, df=degf)*sqrt(T)),
	                     row.names=coef.names)
  round(output,2)
}

appt8[c(1,4,7,9,11,14,16,18),c(2,5,8)] <- 1
cm1 <- combine.mi(model1, n.impute=10)
cm2 <- combine.mi(model2, n.impute=10)
cm3 <- combine.mi(model3, n.impute=10)
appt8[c(2,3,5,6,8,10,12,13,15,17,19),2:4] <- cbind(cm1$OR[-1],cm1$lowerCI[-1],cm1$upperCI[-1])
appt8[c(2,3,5,6,8,10,12,13,15,17,19),5:7] <- cbind(cm2$OR[-1],cm2$lowerCI[-1],cm2$upperCI[-1])
appt8[c(2,3,5,6,8,10,12,13,15,17,19),8:10] <- cbind(cm3$OR[-1],cm3$lowerCI[-1],cm3$upperCI[-1])

### n 
new <- data.nomiss[[1]]
dim(new)

appt8[,1] <- c(dim(new[new$intervention==1,])[1],dim(new[new$intervention==3,])[1],dim(new[new$intervention==4,])[1],
  dim(new[new$agegp==0,])[1], dim(new[new$agegp==1,])[1], dim(new[new$agegp==2,])[1],
  dim(new[new$sex==0,])[1], dim(new[new$sex==1,])[1],
  dim(new[new$vaccine==0,])[1],dim(new[new$vaccine==1,])[1],
  length(unique(new$hhID[new$indexagegp==0])),length(unique(new$hhID[new$indexagegp==1])),length(unique(new$hhID[new$indexagegp==2])),
  length(unique(new$hhID[new$indexsex==0])),length(unique(new$hhID[new$indexsex==1])),
  length(unique(new$hhID[new$indexav==0])),length(unique(new$hhID[new$indexav==1])),
  length(unique(new$hhID[new$co_index==0])),length(unique(new$hhID[new$co_index==1])))

rownames(appt8) <- c("control","hand","handmask","adult","child6-15","child<=5","female","male","not vac","vaccine",
              "adult index","child6-15 index","child<=5 index","female index","male index","no antiviral","antiviral","no co-index","co_index")
colnames(appt8) <- c("n","lab-OR","CI-low","CI-up","cd1-OR","CI-low","CI-up","cd2-OR","CI-low","CI-up")

# seems one household is missing here (n=330)
# b/c hhID=187, all contacts with d_contact==1

appt8
  
# End of script

