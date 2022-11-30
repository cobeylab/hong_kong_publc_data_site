#
# Script to reproduce information in Figure 1 (study flowchart) from:
#
# Cowling BJ, Fung ROP, Cheng KY, Fang VJ, Chan KH, Seto WH, et al.
# Preliminary findings of a randomized trial of non-pharmaceutical
# interventions to prevent influenza transmission in households.
# PLoS ONE, 2008; 3(5):e2101.
# http://dx.doi.org/10.1371/journal.pone.0002101
#
# Last updated by Vicky Fang and Ben Cowling
# January 6, 2009

# Totally 944 index subjects recruited from clinics,
# 198 randomized,
# -70 who refused to participate,
# 128 allocated to control arm,
# 35 allocated to mask arm,
# 35 allocated to hand hygiene arm

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

qv <- read.csv(paste(dir, "clinicdat_h.csv", sep=""), header=TRUE)
hchar <- read.csv(paste(dir, "hchar_h.csv", sep=""), header=TRUE)
hculture <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
baseflu <- read.csv(paste(dir, "adherence_m.csv", sep=""), header=TRUE)

# Step 1: Enrolment

# QuickVue +ve among 198 index subjects
dim(qv[!is.na(qv$QVres) & (qv$QVres==1 | qv$QVres==2),])[1]   # QV: +ve
dim(qv[!is.na(qv$QVres) & qv$QVres==1,])[1]                   # QV: A
dim(qv[!is.na(qv$QVres) & qv$QVres==2,])[1]                   # QV: B


# Step 2: Allocation

# Received allocated intervention
control <- hchar[hchar$intervention==1,]              # Control arm
dim(control)[1]
sum(control$familysize) - dim(control)[1]
median(control$familysize)

mask <- hchar[hchar$intervention==2,]                 # Mask arm
dim(mask)[1]
sum(mask$familysize) - dim(mask)[1]
median(mask$familysize)

hand <- hchar[hchar$intervention==3,]                 # Hand hygiene arm
dim(hand)[1]
sum(hand$familysize) - dim(hand)[1]
median(hand$familysize)


# Step 3: Analysis

hc0 <- hculture[hculture$visit==0 & hculture$member==0,]            # Extract home culture results for index subhects (visit 0)
for (j in 1:nrow(hc0)){
     if( (!is.na(hc0$culture[j])&as.character(hc0$culture[j])>0) )
         hc0$cultpos0[j] <- 1
         else  hc0$cultpos0[j] <- 0
}

hc1 <- hculture[hculture$visit==1 & hculture$member==0,]            # Extract home culture results for index subhects (visit 1)
for (j in 1:nrow(hc1)){
     if( (!is.na(hc1$culture[j])&as.character(hc1$culture[j])>0) )
         hc1$cultpos1[j] <- 1
         else  hc1$cultpos1[j] <- 0
}

hc01 <- merge(hc0[,c(1,2,6)],hc1[,c(1,6)],by="hhID",all.x=T)        # Define 'baseline'(=1) when V0 & V1 culture -ve
hc01$baseline <- 1*(!is.na(hc01$cultpos0) & hc01$cultpos0==0 & !is.na(hc01$cultpos1) & hc01$cultpos1==0)  
hc_comb01 <- hc01[,c(1,5)]
analysis <- merge(hchar,hc_comb01,by="hhID",all.x=T)


# Analysed clusters:
control_als <- analysis[analysis$baseline==0 & analysis$intervention==1,]
dim(control_als)[1]
sum(control_als$familysize) - dim(control_als)[1]           # Participants

mask_als <- analysis[analysis$baseline==0 & analysis$intervention==2,]
dim(mask_als)[1]
sum(mask_als$familysize) - dim(mask_als)[1]                 # Participants

hand_als <- analysis[analysis$baseline==0 & analysis$intervention==3,]
dim(hand_als)[1]
sum(hand_als$familysize) - dim(hand_als)[1]                 # Participants


# Confirmation of influenza infection in index cases will be added in a subsequent script, when RT-PCR results are
# made available in the online dataset.

# End of script
