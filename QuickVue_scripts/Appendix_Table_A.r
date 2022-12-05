#
# R syntax to reproduce information in Appendix Table A (QuickVue performance vs gold standard) from:
#
# Cheng CKY, Cowling BJ, Chan KH, Fang VJ, Seto WH, et al.
# Factors affecting QuickVue Influenza A+B rapid test performance in
# the community setting.
# Diagnostic Microbiology and Infectious Disease, 2009; 65(1): 35-41.
#
# Last updated by Calvin Cheng and Vicky Fang
# Sep 08, 2009


dir <- "http://sph.hku.hk/data/HongKongNPIpilotV3/"

cdc <- read.csv(paste(dir,"QuickVue_data.csv",sep=""))
source ("http://www.hku.hk/bcowling/influenza/QuickVue_scripts/add_groups.r")
source ("http://www.hku.hk/bcowling/influenza/QuickVue_scripts/qv_functions.r")

# extract only culture +ve A and culture +ve B

cultA <- cdc[cdc$cultposA ==1,]
cultB <- cdc[cdc$cultposB ==1,]

#function for getting the sensitivity for A +ve
snres <- function(cdc, entry, parameter) {
  tab <- cdc[c(entry==parameter & !is.na(entry)),]
  sn <- binom.test(x=length(tab$QVpos[tab$QVpos==1]), dim(tab)[1])

output <- data.frame(total = length(tab$QVpos), Value = c(sn[[1]] / sn[[2]]),
                     lowerCI = c(sn[[4]][1]), upperCI = c(sn[[4]][2]))
round(output, 2)
}

# for agegroup
a1 <- snres(cultA, cultA$agegroup, 1)
a2 <- snres(cultA, cultA$agegroup, 2)
a3 <- snres(cultA, cultA$agegroup, 3)
a4 <- snres(cultA, cultA$agegroup, 4)
a5 <- snres(cultA, cultA$agegroup, 5)
a6 <- snres(cultA, cultA$agegroup, 6)
out <- rbind (a1,a2,a3,a4,a5,a6)
row.names(out)<- c("age 0-5","5-10","10-15","16-30","31-50","50+")
apt1p11 <- out

a1 <- snres(cultB, cultB$agegroup, 1)
a2 <- snres(cultB, cultB$agegroup, 2)
a3 <- snres(cultB, cultB$agegroup, 3)
a4 <- snres(cultB, cultB$agegroup, 4)
a5 <- snres(cultB, cultB$agegroup, 5)
a6 <- snres(cultB, cultB$agegroup, 6)
out <- rbind (a1,a2,a3,a4,a5,a6)
apt1p12 <- out

#sex
a1 <- snres(cultA, cultA$male, 1)
a2 <- snres(cultA, cultA$male, 0)
out <- rbind (a1,a2)
row.names(out)<- c("male","female")
apt1p21 <- out

a1 <- snres(cultB, cultB$male, 1)
a2 <- snres(cultB, cultB$male, 0)
out <- rbind (a1,a2)
apt1p22 <- out

#symptom onset
a1 <- snres(cultA, cultA$onsettime, 1)
a2 <- snres(cultA, cultA$onsettime, 2)
a3 <- snres(cultA, cultA$onsettime, 3)
a4 <- snres(cultA, cultA$onsettime, 4)
a5 <- snres(cultA, cultA$onsettime, 5)
out <- rbind (a1,a2,a3,a4,a5)
row.names(out)<- c("onset <12hrs","12-24hrs","24-36hrs","36-48hrs",">48hrs")
apt1p31 <- out

a1 <- snres(cultB, cultB$onsettime, 1)
a2 <- snres(cultB, cultB$onsettime, 2)
a3 <- snres(cultB, cultB$onsettime, 3)
a4 <- snres(cultB, cultB$onsettime, 4)
a5 <- snres(cultB, cultB$onsettime, 5)
out <- rbind (a1,a2,a3,a4,a5)
apt1p32 <- out

#clinic experience
a1 <- snres(cultA, cultA$qvdonegp, 1)
a2 <- snres(cultA, cultA$qvdonegp, 2)
a3 <- snres(cultA, cultA$qvdonegp, 3)
a4 <- snres(cultA, cultA$qvdonegp, 4)
out <- rbind (a1,a2,a3,a4)
row.names(out)<- c("QV <30","30-60","60-90",">90")
apt1p41 <- out

a1 <- snres(cultB, cultB$qvdonegp, 1)
a2 <- snres(cultB, cultB$qvdonegp, 2)
a3 <- snres(cultB, cultB$qvdonegp, 3)
a4 <- snres(cultB, cultB$qvdonegp, 4)
out <- rbind (a1,a2,a3,a4)
apt1p42 <- out

apt1p1 <- rbind(apt1p11,apt1p21,apt1p31,apt1p41)
apt1p2 <- rbind(apt1p12,apt1p22,apt1p32,apt1p42)
cbind(apt1p1,apt1p2)

# END OF SCRIPTS


