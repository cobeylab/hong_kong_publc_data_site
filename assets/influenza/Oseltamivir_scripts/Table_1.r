#
# R syntax to reproduce information for Table 1 from:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

source("http://www.hku.hk/bcowling/influenza/Oseltamivir_scripts/dataframe.r")  

tab1 <- matrix(rep(NA,120),ncol=5)

table(index.tami$av)
tamigp <- index.tami[index.tami$av=="tamiflu",]
nontamigp <- index.tami[index.tami$av=="none",]

tab1[1,c(1,3)] <- c(sum(tamigp$male),sum(nontamigp$male))
tab1[2:5,1] <- table(tamigp$agegp)
tab1[2:5,3] <- table(nontamigp$agegp)
tab1[6,c(1,3)] <- c(sum(tamigp$vaccine),sum(nontamigp$vaccine))
tab1[7:8,1] <- table(tamigp$flu.type)
tab1[7:8,3] <- table(nontamigp$flu.type)
tab1[9,c(1,3)] <- c(sum(tamigp$chronic_disease),sum(nontamigp$chronic_disease))
tab1[10,c(1,3)] <- round(c(mean(tamigp$b_m0_score), mean(nontamigp$b_m0_score)),1)
tab1[11:17,1] <- colSums(tamigp[c(19,5:10)])
tab1[11:17,3] <- colSums(nontamigp[c(19,5:10)])
tab1[18:20,1] <- c(sum(tamigp$onsettime<=2), sum(tamigp$onsettime>=3&tamigp$onsettime<=4), sum(tamigp$onsettime==5))
tab1[18:20,3] <- c(sum(nontamigp$onsettime<=2), sum(nontamigp$onsettime>=3&nontamigp$onsettime<=4), sum(nontamigp$onsettime==5))
tab1[21:24,1] <- colSums(tamigp[c(15,17,16,18)])
tab1[21:24,3] <- colSums(nontamigp[c(15,17,16,18)])

# calculte the proportion
tab1[c(1:9,11:24),2] <- round(tab1[c(1:9,11:24),1]/dim(tamigp)[1]*100,1)
tab1[c(1:9,11:24),4] <- round(tab1[c(1:9,11:24),3]/dim(nontamigp)[1]*100,1)

# function to calculate chi-square test p-value
chi.p <- function(x,y){
   a <- matrix(c(x,y,dim(tamigp)[1]-x,dim(nontamigp)[1]-y),ncol=2)
   round(chisq.test(a)$p.value,2)
}                     

tab1[1,5] <- chi.p(tab1[1,1],tab1[1,3])                
tab1[2,5] <- round(chisq.test(tab1[2:5,c(1,3)])$p.value,2)
tab1[6,5] <- chi.p(tab1[6,1],tab1[6,3])
tab1[7,5] <- round(chisq.test(tab1[7:8,c(1,3)])$p.value,2)
tab1[9,5] <- chi.p(tab1[9,1],tab1[9,3])
tab1[10,5] <- round(t.test(tamigp$b_m0_score,nontamigp$b_m0_score)$p.value,2)
for (i in 11:17){
   tab1[i,5] <- chi.p(tab1[i,1],tab1[i,3])
}
tab1[18,5] <- round(chisq.test(tab1[18:20,c(1,3)])$p.value,2)
for (i in 21:24){
   tab1[i,5] <- chi.p(tab1[i,1],tab1[i,3])
}

# calculate the 95% CI for baseline symptom score
tab1[10,2] <-  paste(round(mean(tamigp$b_m0_score)-1.96*sd(tamigp$b_m0_score)/sqrt(dim(tamigp)[1]),1),"-",
                     round(mean(tamigp$b_m0_score)+1.96*sd(tamigp$b_m0_score)/sqrt(dim(tamigp)[1]),1),sep="")
tab1[10,4] <-  paste(round(mean(nontamigp$b_m0_score)-1.96*sd(nontamigp$b_m0_score)/sqrt(dim(nontamigp)[1]),1),"-",
                     round(mean(nontamigp$b_m0_score)+1.96*sd(nontamigp$b_m0_score)/sqrt(dim(nontamigp)[1]),1),sep="")

colnames(tab1) <- c("Oseltamivir","%","No antiviral","%","p-value")
rownames(tab1) <- c("Male","Age<=5","Age6-12","Age13-17","Age18+","vac","FluA","FluB","chronic","symptom score",
                    "Fever","Headache","Sore throat","Cough","Muscle pain","Coryza","Phlegm","Delay<=24hr","Delay24-48hr",
                    "Delay>48hr","Antibiotics","Antipyretics","Antihistamines","Steroids")
tab1

#
# End of script
#
