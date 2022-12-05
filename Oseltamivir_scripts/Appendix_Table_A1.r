#
# R syntax to reproduce information for Appendix Table A1 from:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

source("http://www.hku.hk/bcowling/influenza/Oseltamivir_scripts/dataframe.r")

appt1 <- matrix(rep(NA,150),ncol=10)

# Imputing age by relation, occupation and schooling
for (i in 1:nrow(trans.c)){
    if (trans.c$hhID[i]==7001 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==7090 & trans.c$member[i]==3) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==7093 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==7148 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8148 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8148 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8148 & trans.c$member[i]==3) {trans.c$age[i] <-6}
    if (trans.c$hhID[i]==8273 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8273 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8283 & trans.c$member[i]==3) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8301 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==2) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==3) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==4) {trans.c$age[i] <-5}
    if (trans.c$hhID[i]==8317 & trans.c$member[i]==5) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8398 & trans.c$member[i]==1) {trans.c$age[i] <-18}
    if (trans.c$hhID[i]==8398 & trans.c$member[i]==2) {trans.c$age[i] <-18}
}
index.hh2 <- index.hh[c("hhID","av","agegp","male","vaccine","flu.type","chronic_disease","b_m0_score","fever",
                        "antibiotics", "antipyretic", "antihistamine", "steroid")]
trans.data <- merge(trans.c,index.hh2[c("hhID","av")],all.x=T)
trans.data$agegp <- cut(trans.data$age,c(-0.1,5,12,17,100))

# Total n
index.tami <- index.hh2[index.hh2$av=="tamiflu",]
index.ntami <- index.hh2[index.hh2$av=="none",]
contact.tami <- trans.data[trans.data$av=="tamiflu",]
contact.ntami <- trans.data[trans.data$av=="none",]

# male
appt1[1,c(1,3,6,8)] <- c(sum(index.tami$male), sum(index.ntami$male), sum(contact.tami$male), sum(contact.ntami$male))

# age group
appt1[2:5,c(1,3,6,8)] <- cbind(table(index.tami$agegp), table(index.ntami$agegp),
                               table(contact.tami$agegp), table(contact.ntami$agegp))

# vaccine status
appt1[6,c(1,3,6,8)] <-  c(sum(index.tami$vaccine), sum(index.ntami$vaccine),
                         sum(contact.tami$vaccine), sum(contact.ntami$vaccine,na.rm=T))

# flu type
appt1[7:8,c(1,3)] <- cbind(table(index.tami$flu.type),table(index.ntami$flu.type))

# chronic condition
appt1[9,c(1,3,6,8)] <- c(sum(index.tami$chronic_disease), sum(index.ntami$chronic_disease),
                         sum(contact.tami$chronic_disease,na.rm=T), sum(contact.ntami$chronic_disease,na.rm=T))

# baseline symptom score
appt1[10,c(1,3)] <- round(c(mean(index.tami$b_m0_score),mean(index.ntami$b_m0_score)),1)

# fever
appt1[11,c(1,3)] <- c(sum(index.tami$fever),sum(index.ntami$fever))

# anti-
appt1[12:15,c(1,3)] <- cbind(colSums(index.tami[10:13]),colSums(index.ntami[10:13]))

appt1[,2] <- round(appt1[,1]/nrow(index.tami)*100,1)
appt1[,4] <- round(appt1[,3]/nrow(index.ntami)*100,1)
appt1[,7] <- round(appt1[,6]/nrow(contact.tami)*100,1)
appt1[,9] <- round(appt1[,8]/nrow(contact.ntami)*100,1)

# function to calculate chi-square test p-value
chi.p.i <- function(x,y){
   a <- matrix(c(x,y,nrow(index.tami)-x,nrow(index.ntami)-y),ncol=2)
   round(chisq.test(a)$p.value,2)
}
chi.p.c <- function(x,y){
   a <- matrix(c(x,y,nrow(contact.tami)-x,nrow(contact.ntami)-y),ncol=2)
   round(chisq.test(a)$p.value,2)
}

appt1[1,c(5,10)] <- c(chi.p.i(appt1[1,1],appt1[1,3]), chi.p.c(appt1[1,6],appt1[1,8]))
appt1[2,5] <- round(chisq.test(appt1[2:5,c(1,3)])$p.value,2)
appt1[2,10] <- round(chisq.test(appt1[2:5,c(6,8)])$p.value,2)
appt1[6,c(5,10)] <- c(chi.p.i(appt1[6,1],appt1[6,3]), chi.p.c(appt1[6,6],appt1[6,8]))
appt1[7,5] <- chi.p.i(appt1[7,1],appt1[7,3])
appt1[9,c(5,10)] <- c(chi.p.i(appt1[9,1],appt1[9,3]), chi.p.c(appt1[9,6],appt1[9,8]))
appt1[10,5] <- round(t.test(index.tami$b_m0_score,index.ntami$b_m0_score)$p.value,2)
for(i in 11:15){
  appt1[i,5] <- chi.p.i(appt1[i,1],appt1[i,3])
}

# 95% CI for baseline symptom score
se1 <- sd(index.tami$b_m0_score)/sqrt(nrow(index.tami))
se2 <- sd(index.ntami$b_m0_score)/sqrt(nrow(index.ntami))
appt1[10,2] <- paste(round(mean(index.tami$b_m0_score)-1.96*se1,1),"-",round(mean(index.tami$b_m0_score)+1.96*se1,1),sep="")
appt1[10,4] <- paste(round(mean(index.ntami$b_m0_score)-1.96*se2,1),"-",round(mean(index.ntami$b_m0_score)+1.96*se2,1),sep="")

colnames(appt1) <- c("Oseltamivir.Index","%","No antiviral","%","p-value", "Oseltamivir.Contact","%","No antiviral","%","p-value")
rownames(appt1) <- c("Male","Age<=5","Age6-12","Age13-17","Age18+","vac","FluA","FluB","chronic","symptom score",
                    "Fever","Antibiotics","Antipyretics","Antihistamines","Steroid")
appt1

#
# End of script
#

