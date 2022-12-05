
#
# R syntax to reproduce information for Table 1 from:
#
# Lau LLH, Cowling BJ, Fang VJ, Chan KH, Lau EHY, et al.
# Viral shedding and clinical illness in naturally acquired influenza virus infections
# JID, 2010 (in press).
#
# Last updated by Vicky Fang and Lincoln Lau
# Jan 30, 2010


dir <- "http://sph.hku.hk/data/HongKongNPIstudyV4/"
source("http://www.hku.hk/bcowling/influenza/Vshed_scripts/JID_dataframe.r")
clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))

tab1 <- matrix(rep(NA,28),ncol=4)

plotdata <- merge(plotdata,clinic[c(2,13)],all.x=T)
plotdata$flu.type <- "A"
plotdata$flu.type[plotdata$QVres==2] <- "B"
tab1data <- plotdata
table(tab1data$flu.type)

names(tab1data)[6] <- "day"
tab1data <- merge(tab1data, symp[c(1:3,5:11)], by=c("hhID", "member", "day"))

# make the table

tab1[1,c(1,3)] <- table(tab1data$flu.type[tab1data$rnose==1])
tab1[2,c(1,3)] <- table(tab1data$flu.type[tab1data$cough==1])
tab1[3,c(1,3)] <- table(tab1data$flu.type[tab1data$sthroat==1])
tab1[4,c(1,3)] <- table(tab1data$flu.type[tab1data$headache==1])
tab1[5,c(1,3)] <- table(tab1data$flu.type[tab1data$phlegm==1])
tab1[6,c(1,3)] <- table(tab1data$flu.type[tab1data$pmuscle==1])
tab1[7,c(1,3)] <- table(tab1data$flu.type[tab1data$fever==1])

tab1[,2] <- round(tab1[,1]/nrow(tab1data[tab1data$flu.type=="A",]),2)
tab1[,4] <- round(tab1[,3]/nrow(tab1data[tab1data$flu.type=="B",]),2)

rownames(tab1) <- c("Runny nose","Cough","Sore throat","Headache","Phlegm","Myalgia","Fever")
colnames(tab1) <- c("Flu A","%","Flu B","%")
tab1

#
# End of script
#
