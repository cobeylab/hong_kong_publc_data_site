#
# Script to reproduce information in the manuscript text of:
#
# Cowling BJ, Fung ROP, Cheng KY, Fang VJ, Chan KH, Seto WH, et al.
# Preliminary findings of a randomized trial of non-pharmaceutical
# interventions to prevent influenza transmission in households.
# PLoS ONE, 2008; 3(5):e2101.
# http://dx.doi.org/10.1371/journal.pone.0002101
#
# Last updated by Vicky Fang and Ben Cowling
# January 5, 2009


#
# RESULTS | Main outcomes | Overall laboratory-confirmed SAR

dir <- "http://sph.hku.hk/data/HongKongNPIpilotV2/"

hc <- read.csv(paste(dir, "home_culture.csv", sep=""), header=TRUE)
mark <- data.frame(hhID = unique(baseflu$hhID))
hc <- merge(hc,mark,by="hhID",all.y=TRUE)
hc <- hc[order(hc$hhID,hc$member,hc$visit),]
for (i in 1:nrow(hc)){
     if(!is.na(hc$PCR[i]) & (is.na(hc$culture[i])|(!is.na(hc$culture[i])&hc$culture[i]==0)) ) hc$culture[i] <- hc$PCR[i]
}
hc <- hc[-5]

hculture <- data.frame(hhID = baseflu$hhID, member = baseflu$member)

hc.temp <- reshape(hc, timevar="visit", idvar=c("hhID","member"), direction="wide", v.names="culture")
hculture <- merge(hculture,hc.temp, by=c("hhID","member"), all.x=TRUE)
names(hculture) <- c("hhID","member","V0","V1","V2","V3","V4")

## exd_index: none of V0/V1 culture is A/B; d_index: both V0 and V1 culture is 0; Contact exclusion: V1 culture is A/B
for (i in 1:nrow(hculture)){
     if(hculture$member[i]==0 & ( (is.na(hculture$V0[i]) & is.na(hculture$V1[i])) | (is.na(hculture$V0[i]) & hculture$V1[i]==0)
				  | (is.na(hculture$V1[i]) & hculture$V0[i]==0) | (hculture$V0[i]==0 & hculture$V1[i]==0) )) 
              {hculture$exd_index[i]<-1}    else {hculture$exd_index[i]<-0}
     if(hculture$member[i]==0 &  !is.na(hculture$V0[i]) & !is.na(hculture$V1[i]) & hculture$V0[i]==0 & hculture$V1[i]==0)
              {hculture$d_index[i]<-1}      else {hculture$d_index[i]<-0}
     if(hculture$member[i]!=0 & ( !is.na(hculture$V1[i]) & (hculture$V1[i]=="A" | hculture$V1[i]=="B") ))
              {hculture$exd_contact[i]=1}   else{hculture$exd_contact[i]=0}
}

## Define the household which should be excluded as long as index in this hh should be excluded
exd_index <- hculture[hculture$member==0,c(1,8)]   # for calculating secondary cases
d_index <- hculture[hculture$member==0,c(1,9)]     # for calculating SAR

dim(hculture)
hculture <- merge(hculture[,-8], exd_index)
hculture <- merge(hculture[,-8], d_index)
dim(hculture)

for (i in 1:nrow(hculture)){
    if ( hculture$exd_index[i]==1 | hculture$exd_contact[i] ==1)
       {hculture$exclude[i] <-1}
       else  {hculture$exclude[i] <-0}
}

## Define secondary cases
for (i in 1:nrow(hculture)){
    if ( hculture$member[i] != 0 & hculture$exclude[i] == 0 & !is.na(hculture$V1[i]) & ( (hculture$V2[i] !=0 & !is.na(hculture$V2[i])) |
               (hculture$V3[i] !=0 & !is.na(hculture$V3[i])) | (hculture$V4[i] !=0 & !is.na(hculture$V4[i])) ) )
              {hculture$labsedcase[i] <- 1}
	 else {hculture$labsedcase[i] <- 0}
}

total <- dim(hculture[hculture$member!=0&hculture$d_index==0,])[1]
sum(hculture$labsedcase)/total
binom.test(sum(hculture$labsedcase),total)$conf[1:2]


#
# RESULTS | Main outcomes | Within-household correlation

hculture.cc <- merge(hculture,hchar[,1:3])
hculture.cc$intervention <- factor(hculture.cc$intervention, levels=1:3, labels=c("control", "mask", "hand"))
names(hculture.cc)[dim(hculture.cc)[2]]<- "hhsize"
hculture.cc <- hculture.cc[hculture.cc$member!=0&hculture.cc$d_index==0,]
hculture.cc$hhsize <- hculture.cc$hhsize -1
hculture.cc$hhsize[hculture.cc$member>1] <- 0

## Function 'icc' to estimate the intracluster correlation coefficient
icc <- function(data_id,data_arm,data_clustersize,data_success){
 M <- length(data_id)
 Y <- sum(data_success)
 P <- Y/M
 K <- length(unique(data_id))

 arm <- c("control","mask","hand")   # interventions
 Yij <- Mij <- Pij <- hhID <- matrix(nrow=3,ncol=length(unique(data_id[data_arm==arm[1]])))
 for (i in 1:length(arm)){
      hhID[i,1:length(unique(data_id[data_arm==arm[i]]))] <- as.character(unique(data_id[data_arm==arm[i]]))
 }

 for (i in 1:length(arm)){
      for (j in 1:length(unique(data_id[data_arm==arm[i]]))){
        Yij[i,j] <- sum(data_success[data_id==hhID[i,j]])
        Mij[i,j] <- sum(data_clustersize[data_id==hhID[i,j]])
        Pij[i,j] <- Yij[i,j]/Mij[i,j]
      }
 }

 Y_i <- M_i <- P_i <- m_Ai <- rho_i <- rep(NA,length(arm))
 output <- data.frame(Y_i,M_i,m_Ai,P_i,rho_i)
 for (i in 1:length(arm)){
     output$Y_i[i] <- sum(na.exclude(Yij[i,]))
     output$M_i[i] <- sum(na.exclude(Mij[i,]))
     output$P_i[i] <- output$Y_i[i]/output$M_i[i]
     output$m_Ai[i] <- sum(na.exclude(Mij[i,]^2))/output$M_i[i]
 }

 m_0 <- (M-sum(output$m_Ai))/(K-2)
 MSW <- sum(rowSums(Mij*Pij*(1-Pij),na.rm = TRUE)) / (M-K)
 MSC <- sum(rowSums(Mij*(Pij-output$P_i)*(Pij-output$P_i),na.rm = TRUE)) / (K-2)
 rho <- (MSC-MSW)/(MSC+(m_0-1)*MSW)
 output$rho_i <- rho
 output
}

icc(hculture.cc$hhID,hculture.cc$intervention,hculture.cc$hhsize,hculture.cc$labsedcase)$rho_i[1] 


#
# RESULTS | Ancillary analyses | Self-reported face mask use

aftervisit <- read.csv(paste(dir, "adherence_m.csv", sep=""), header=TRUE)
aftervisit <- merge(aftervisit,hchar[,1:2],by="hhID",all.x=TRUE)

# For mask group
index_m <- aftervisit[aftervisit$intervention==2&aftervisit$member==0,]
dim(index_m[!is.na(index_m$when_mask)&index_m$when_mask<=2,])[1]/dim(index_m)[1]
contact_m <- aftervisit[aftervisit$intervention==2&aftervisit$member>0,]
dim(contact_m[!is.na(contact_m$when_mask)&contact_m$when_mask<=2,])[1]/dim(contact_m)[1]

# For control group
index_c <- aftervisit[aftervisit$intervention==1&aftervisit$member==0,]
dim(index_c[!is.na(index_c$when_mask)&index_c$when_mask<=2,])[1]/dim(index_c)[1]
contact_c <- aftervisit[aftervisit$intervention==1&aftervisit$member>0,]
dim(contact_c[!is.na(contact_c$when_mask)&contact_c$when_mask<=2,])[1]/dim(contact_c)[1]

# For hand hygiene group
index_h <- aftervisit[aftervisit$intervention==3&aftervisit$member==0,]
dim(index_h[!is.na(index_h$when_mask)&index_h$when_mask<=2,])[1]/dim(index_h)[1]
contact_h <- aftervisit[aftervisit$intervention==3&aftervisit$member>0,]
dim(contact_h[!is.na(contact_h$when_mask)&contact_h$when_mask<=2,])[1]/dim(contact_h)[1]


#
# RESULTS | Ancillary analyses | Measured use of face masks

mask <- aftervisit[aftervisit$intervention==2,]
mask_index_m <- mask[mask$member==0,]
quantile(mask_index_m$mask_usage,probs=c(0.25,0.5,0.75),na.rm=TRUE)
mask_contact_m <- mask[mask$member>0,]
quantile(mask_contact_m$mask_usage,probs=c(0.25,0.5,0.75),na.rm=TRUE)


#
# RESULTS | Ancillary analyses | Self-reported hand hygiene

# For hand hygiene group
index_h <- aftervisit[aftervisit$intervention==3&aftervisit$member==0,]
dim(index_h[!is.na(index_h$washhand)&index_h$washhand<=2,])[1]/dim(index_h)[1]
contact_h <- aftervisit[aftervisit$intervention==3&aftervisit$member>0,]
dim(contact_h[!is.na(contact_h$washhand)&contact_h$washhand<=2,])[1]/dim(contact_h)[1]

# For control group
index_c <- aftervisit[aftervisit$intervention==1&aftervisit$member==0,]
dim(index_c[!is.na(index_c$washhand)&index_c$washhand<=2,])[1]/dim(index_c)[1]
contact_c <- aftervisit[aftervisit$intervention==1&aftervisit$member>0,]
dim(contact_c[!is.na(contact_c$washhand)&contact_c$washhand<=2,])[1]/dim(contact_c)[1]

# For mask group
index_m <- aftervisit[aftervisit$intervention==2&aftervisit$member==0,]
dim(index_m[!is.na(index_m$washhand)&index_m$washhand<=2,])[1]/dim(index_m)[1]
contact_m <- aftervisit[aftervisit$intervention==2&aftervisit$member>0,]
dim(contact_m[!is.na(contact_m$washhand)&contact_m$washhand<=2,])[1]/dim(contact_m)[1]


#
# RESULTS | Ancillary analyses | Measured use of alcohol hand rub and liquid hand soap

soap <- read.csv(paste(dir, "adherence_h.csv", sep=""), header=TRUE)

quantile(soap$alcohol_usage,probs=c(0.25,0.5,0.75),na.rm=TRUE)
quantile(soap$soap_usage,probs=c(0.25,0.5,0.75),na.rm=TRUE)

quantile(aftervisit$smallgel_usage[aftervisit$member==0],probs=c(0.25,0.5,0.75),na.rm=TRUE)
quantile(aftervisit$smallgel_usage[aftervisit$member>0],probs=c(0.25,0.5,0.75),na.rm=TRUE)


# End of script
