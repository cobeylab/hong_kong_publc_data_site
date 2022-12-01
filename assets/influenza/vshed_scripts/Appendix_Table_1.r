
#
# R syntax to reproduce information for Appendix Table 1 from:
#
# Lau LLH, Cowling BJ, Fang VJ, Chan KH, Lau EHY, et al.
# Viral shedding and clinical illness in naturally acquired influenza virus infections
# JID, 2010 (in press).
#
# Last updated by Vicky Fang and Lincoln Lau
# Jan 30, 2010

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"
source("http://www.hku.hk/bcowling/influenza/Vshed_scripts/JID_dataframe.r")
demog <- read.csv(paste(dir, "adherence_m.csv", sep=""))
appt1 <- matrix(rep(NA,24),ncol=4)

appt1data <- merge(contact[c(1:2,9)], demog[c(1:5,15:17)], by=c("hhID","member"), all.x=T)
rev(table(appt1data$analyse))

# make the table

appt1[1,c(3,1)] <- table(appt1data$analyse[appt1data$male==1])
appt1[2,c(3,1)] <- table(appt1data$analyse[appt1data$age>=16])
appt1[3,c(3,1)] <- table(appt1data$analyse[appt1data$ever_smoke==1])
appt1[4,c(3,1)] <- table(appt1data$analyse[appt1data$vaccine08==1])
appt1[5,c(3,1)] <- table(appt1data$analyse[appt1data$vaccine07==1])
appt1[6,c(3,1)] <- table(appt1data$analyse[appt1data$vaccine06==1])

appt1[,2] <- round(appt1[,1]/nrow(appt1data[appt1data$analyse==1,]),2)
appt1[,4] <- round(appt1[,3]/nrow(appt1data[appt1data$analyse==0,]),2)

rownames(appt1) <- c("Male","Age>=16","Ever smoke","Vaccine07-08","Vaccine06-07","Vaccine05-06")
colnames(appt1) <- c("Included","%","Excluded","%")
appt1

#
# End of script
#
