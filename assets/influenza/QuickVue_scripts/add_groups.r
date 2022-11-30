#
# R syntax to reproduce information for functions from:
#
# Cheng CKY, Cowling BJ, Chan KH, Fang VJ, Seto WH, et al.
# Factors affecting QuickVue Influenza A+B rapid test performance in
# the community setting.
# Diagnostic Microbiology and Infectious Disease, 2009; 65(1): 35-41.
#
# Last updated by Calvin Cheng and Vicky Fang
# Sep 08, 2009

# copied to tab 1

cdc$agegroup[cdc$age <=5 &!is.na(cdc$age)]<- 1
cdc$agegroup[cdc$age >5 & cdc$age <=10 &!is.na(cdc$age)]<-2
cdc$agegroup[cdc$age >10 & cdc$age <=15 &!is.na(cdc$age)]<-3
cdc$agegroup[cdc$age >15 & cdc$age <=30 &!is.na(cdc$age)]<-4
cdc$agegroup[cdc$age >30 & cdc$age <=50 &!is.na(cdc$age)]<-5
cdc$agegroup[cdc$age >50 &!is.na(cdc$age)]<-6

#07/06/07 add adult /child group
cdc$isadult[cdc$age <16 &!is.na(cdc$age)]<- 0
cdc$isadult[cdc$age >=16 &!is.na(cdc$age)]<- 1


cdc$onsettime[cdc$onsettime==9] <- NA
#Add entry onsetday for day after symptom onset
cdc$onsetday[cdc$onsettime ==1]<- 0 
cdc$onsetday[cdc$onsettime ==2 | cdc$onsettime ==3]<- 1
cdc$onsetday[cdc$onsettime ==4]<- 2

#Add entry fcs for WHO clinic definition fever + cough | sore throat 
cdc$fcs <- 1*(cdc$measure_fever==1 & (cdc$cough==1 | cdc$sthroat==1))

#071123 combine joint and muscle and add new entry jointmuscle
cdc$jointmuscle[cdc$sjoint ==1 |cdc$pmuscle==1] <-1
cdc$jointmuscle[is.na(cdc$jointmuscle)] <-0

#Add entry QVpos for QV positive to 1 and negative and na to 0 
cdc$QVpos[cdc$QVres ==1 |cdc$QVres ==2 & !is.na(cdc$QVres)] <- 1
cdc$QVpos[cdc$QVres ==3 |cdc$QVres ==4 | cdc$QVres ==9 & !is.na(cdc$QVres)] <- 0

#071123 Add entry QVposA and QVposB for only QV is A positive or B positive
cdc$QVposA[cdc$QVres ==1 & !is.na(cdc$QVres)] <- 1
cdc$QVposA[is.na(cdc$QVposA)] <- 0

cdc$QVposB[cdc$QVres ==2 & !is.na(cdc$QVres)] <- 1
cdc$QVposB[is.na(cdc$QVposB)] <- 0

#Add entry cultposA and cultposB for only culture is A positive or B positive
cdc$cultposA[cdc$culture =="A" & !is.na(cdc$culture)] <- 1
cdc$cultposA[is.na(cdc$cultposA)] <- 0

cdc$cultposB[cdc$culture =="B" & !is.na(cdc$culture)] <- 1
cdc$cultposB[is.na(cdc$cultposB)] <- 0

#Add entry cultpos for culture positive to 1 and negative to 0
cdc$cultpos[(cdc$culture =="A" |cdc$culture =="B") & !is.na(cdc$culture)] <- 1
cdc$cultpos[cdc$culture =="0" & !is.na(cdc$culture)] <- 0

#071204 add month column
cdc$month <- as.numeric(substr(as.character(cdc$date),1,1))

#071204 add clinic quickvue test group done column
abc <- table(cdc$clinic)
for (i in 1:dim(abc)){
  if (abc[i][[1]] < 30) abc[i][[1]] <- 1
  else if (abc[i][[1]] < 60) abc[i][[1]] <- 2
  else if (abc[i][[1]] < 90) abc[i][[1]] <- 3
  else abc[i][[1]] <- 4
}

checkgp <- function(input, tab){
for (i in 1:dim(tab)){
  if (input == names(tab)[i] )  return(tab[i][[1]])
 }
}  

for (i in 1:dim(cdc)[1]){
cdc$qvdonegp[i] <- checkgp(cdc$clinic[i], abc)
}

# END OF SCRIPTS
