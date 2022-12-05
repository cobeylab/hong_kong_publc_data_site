#
# Script to reproduce information in Figure 1 (flow chart) from:
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

q1data <- read.csv(paste(dir, "clinicdat_h.csv", sep=""))
arm <- read.csv(paste(dir, "randomarm_407.csv", sep=""))
av <- read.csv(paste(dir, "antiviral_m.csv", sep=""))

## For randomized index subjects (n=407)

q1 <- q1data
arm$hhID <- as.numeric(substr(arm$hhID,5,7))

tab1 <- merge(arm[2:3],q1,by="scrID",all.x=TRUE)
tab1$onsettime[is.na(tab1$onsettime)] <- 9
tab1 <- merge(tab1,av[av$member==0,c(1,3)],by="hhID",all.x=TRUE)
tab1$av <- as.character(tab1$av)
tab1$av[is.na(tab1$av)] <- 0

tab1_c <- tab1[tab1$intervention==1,]
tab1_h <- tab1[tab1$intervention==3,]
tab1_m <- tab1[tab1$intervention==4,]

# Random allocated - total
dim(tab1)
# QuickVue result: A vs B
table(tab1$QVres)

# Random allocated - three interventions
c(dim(tab1_c)[1],dim(tab1_h)[1],dim(tab1_m)[1])

## Home visit (322 in total)
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""))
hchar$analyzed <- hculture$analyzed[hculture$member==0]

hchar$delay <- hchar$onset_v1_delay-hchar$onsettime
hchar <- hchar[c(1,2,3,10,11)] # hhID, intervention, hhsize, analyzed, delay
tab1.2 <- merge(tab1[-3],hchar,by="hhID")

tab1_c <- tab1.2[tab1.2$intervention==1,]
tab1_h <- tab1.2[tab1.2$intervention==3,]
tab1_m <- tab1.2[tab1.2$intervention==4,]
# received allocated intervention - number of households
c(dim(tab1_c)[1],dim(tab1_h)[1],dim(tab1_m)[1])
# received allocated intervention - number of contacts
c(sum(tab1_c$familysize)-dim(tab1_c)[1], sum(tab1_h$familysize)-dim(tab1_h)[1], sum(tab1_m$familysize)-dim(tab1_m)[1])

# received allocated intervention - household size (IQR)
round(quantile(tab1_c$familysize,c(0.5,0.25,0.75)))
round(quantile(tab1_h$familysize,c(0.5,0.25,0.75)))
round(quantile(tab1_m$familysize,c(0.5,0.25,0.75)))

## Analyzed hhs (259 in total)
analyze <- tab1.2[tab1.2$analyzed==1,]

tab1_c <- analyze[analyze$intervention==1,]
tab1_h <- analyze[analyze$intervention==3,]
tab1_m <- analyze[analyze$intervention==4,]
# analyzed hhs - number of households
c(dim(tab1_c)[1],dim(tab1_h)[1],dim(tab1_m)[1])
# analyzed hhs - number of contacts
c(sum(tab1_c$familysize)-dim(tab1_c)[1], sum(tab1_h$familysize)-dim(tab1_h)[1], sum(tab1_m$familysize)-dim(tab1_m)[1])

# analyzed hhs - household size (IQR)
round(quantile(tab1_c$familysize,c(0.5,0.25,0.75)))
round(quantile(tab1_h$familysize,c(0.5,0.25,0.75)))
round(quantile(tab1_m$familysize,c(0.5,0.25,0.75)))

# Exclude from analysis

nanalyze <- tab1.2[tab1.2$analyzed==0,]

tab1_c <- nanalyze[nanalyze$intervention==1,]
tab1_h <- nanalyze[nanalyze$intervention==3,]
tab1_m <- nanalyze[nanalyze$intervention==4,]
# exclude from analysis - number of households
c(dim(tab1_c)[1],dim(tab1_h)[1],dim(tab1_m)[1])
# exclude from analysis - household size (IQR)
round(quantile(tab1_c$familysize,c(0.5,0.25,0.75)))
round(quantile(tab1_h$familysize,c(0.5,0.25,0.75)))
round(quantile(tab1_m$familysize,c(0.5,0.25,0.75)))

# End of script
