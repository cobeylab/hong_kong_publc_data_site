#
# Script to reproduce information in Appendix Table 1 from:
#
# Cowling BJ, Chan KH, Fang VJ, Cheng CKY, Fung ROP, Wai W, et al.
# Facemasks and hand hygiene to prevent influenza transmission 
# in households, a randomized trial.
# Annals of Internal Medicine, 2009 (in press).
#
# Last updated by Vicky Fang and Ben Cowling
# Sep 08, 2009

dir <- "http://sph.hku.hk/data/HongKongNPIstudyV3/"
source("http://www.hku.hk/bcowling/influenza/NPImain_scripts/Analyzed_hh.r")

dat <- read.csv(paste(dir, "hchar_h.csv", sep=""))
clinic <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
dat$analyzed <- hculture$analyzed[hculture$member==0]
dat <- merge(dat,clinic[c("hhID","onsettime")],by="hhID",all.x=TRUE)
dat$onsettime <- dat$onsettime*12
dat <- dat[dat$analyzed==1,]

#function for getting the delay days
delayday <- function(dataset){
  temp1 <- as.matrix(dataset[dataset<=12])
  d0 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>12&dataset<=24])
  d1 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>24&dataset<=36])
  d2 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>36&dataset<=48])
  d3 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>48&dataset<=60])
  d4 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>60&dataset<=72])
  d5 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>72&dataset<=84])
  d6 <- dim(na.omit(temp1))[1]
  temp1 <- as.matrix(dataset[dataset>84&dataset<=96])
  d7 <- dim(na.omit(temp1))[1]
  result <- c(d0,d1,d2,d3,d4,d5,d6,d7)
  return (result)
 }
#end

appt1 <- matrix(rep(NA,48),ncol=6,byrow=FALSE)

# 1. the delay of symptom onset to random assignment
appt1[,1] <- delayday(dat$onsettime)
appt1[,2] <- round(appt1[,1]/dim(dat)[1]*100)

# 2. the delay of random assignment to intervention
between <-  dat$onset_v1_delay - dat$onsettime
appt1[,3] <- delayday(between)
appt1[,4] <- round(appt1[,3]/dim(dat)[1]*100)

# 3. the delay from symptom onset to intervention (in hours)
appt1[,5] <- delayday(dat$onset_v1_delay)
appt1[,6] <- round(appt1[,5]/dim(dat)[1]*100)

colnames(appt1) <- c("onset-random","%","random-intervention","%","onset-intervention","%")
rownames(appt1) <- c("0-12h","12-24h","24-36h","36-48h","48-60h","60-72h","72-84h","84-96h")
appt1

# End of script


