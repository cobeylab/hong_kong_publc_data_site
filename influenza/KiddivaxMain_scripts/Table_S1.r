#
# R syntax to reproduce Table S1 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

dir <- "http://sph.hku.hk/data/KiddivaxMainV1/"

random <- read.csv(paste(dir, "randomcode.csv", sep=""))
arr <- read.csv(paste(dir, "ARR.csv", sep=""))
arr <- merge(arr,random,all.x=TRUE)

tab <- matrix(rep(NA,10*13), ncol=13, dimnames=list(
              c("Fever","Chills","Tired","Headache","Cough","Muscle pain","Swelling","Redness","Bruising","Pain"),
              c("TIV: Mild","%","Moderate","%","Severe","%","Placebo: Mild","%","Moderate","%","Severe","%","p-value")))

arr <- arr[arr$returned==1,-2]
arr <- data.frame(lapply(arr,function(x,...){x[is.na(x)] <- 0 ; x}))
arr$pain <- arr$bruise <- arr$red <- arr$swell <-
            arr$mpain <- arr$cough <- arr$headache <- arr$tired <- arr$chill <-arr$fever <-  NA
for (i in 1:nrow(arr)){
  for (j in 47:56){
   arr[i,j] <- max(arr[i,(j-47)*4+6:9],na.rm=TRUE)
  }
}

arr.new <- arr[c(1,46:56)]
total <- c(dim(arr.new[arr.new$intervention=="TIV",])[1],dim(arr.new[arr.new$intervention=="placebo",])[1])

##
fever <- table(arr.new$fever,arr.new$intervention=="placebo")
tab[1,1:6*2-1] <- c(fever[-1,]); tab[1,1:6*2] <- c(round(prop.table(fever,2),2)[-1,]); tab[1,13] <- round(fisher.test(fever)$p.value,2)

chill <- rbind(table(arr.new$chill,arr.new$intervention=="placebo"),c(0,0))
tab[2,1:6*2-1] <- c(chill[-1,]); tab[2,1:6*2] <- c(round(prop.table(chill,2),2)[-1,]); tab[2,13] <- round(fisher.test(chill)$p.value,2)

tired <- table(arr.new$tired,arr.new$intervention=="placebo")
tab[3,1:6*2-1] <- c(tired[-1,]); tab[3,1:6*2] <- c(round(prop.table(tired,2),2)[-1,]); tab[3,13] <- round(fisher.test(tired)$p.value,2)

headache <- table(arr.new$headache,arr.new$intervention=="placebo")
tab[4,1:6*2-1] <- c(headache[-1,]); tab[4,1:6*2] <- c(round(prop.table(headache,2),2)[-1,]); tab[4,13] <- round(fisher.test(headache)$p.value,2)

cough <- table(arr.new$cough,arr.new$intervention=="placebo")
tab[5,1:6*2-1] <- c(cough[-1,]); tab[5,1:6*2] <- c(round(prop.table(cough,2),2)[-1,]); tab[5,13] <- round(fisher.test(cough)$p.value,2)

mpain <- table(arr.new$mpain,arr.new$intervention=="placebo")
tab[6,1:6*2-1] <- c(mpain[-1,]); tab[6,1:6*2] <- c(round(prop.table(mpain,2),2)[-1,]); tab[6,13] <- round(fisher.test(mpain)$p.value,2)

swell <- table(arr.new$swell,arr.new$intervention=="placebo")
tab[7,1:6*2-1] <- c(swell[-1,]); tab[7,1:6*2] <- c(round(prop.table(swell,2),2)[-1,]); tab[7,13] <- round(fisher.test(swell)$p.value,2)

red <- rbind(table(arr.new$red,arr.new$intervention=="placebo"),c(0,0))[c(1,2,4,3),]
tab[8,1:6*2-1] <- c(red[-1,]); tab[8,1:6*2] <- c(round(prop.table(red,2),2)[-1,]); tab[8,13] <- round(fisher.test(red)$p.value,2)

bruise <- rbind(table(arr.new$bruise,arr.new$intervention=="placebo"),c(0,0))
tab[9,1:6*2-1] <- c(bruise[-1,]); tab[9,1:6*2] <- c(round(prop.table(bruise,2),2)[-1,]); tab[9,13] <- round(fisher.test(bruise)$p.value,2)

pain <- table(arr.new$pain,arr.new$intervention=="placebo")
tab[10,1:6*2-1] <- c(pain[-1,]); tab[10,1:6*2] <- c(round(prop.table(pain,2),2)[-1,]); tab[10,13] <- round(fisher.test(pain)$p.value,2)

tab

#
# End of script.
#

