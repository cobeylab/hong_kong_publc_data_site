#
# Script to reproduce information in Appendix Table 10 from:
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
dat$analyzed <- hculture$analyzed[hculture$member==0]
dat <- dat[dat$analyzed==1,]

# get those only participants who began intervention within 36 hrs of symptom onset
hhID36 <- dat$hhID[dat$onset_v1_delay <=36]

# Proportion of subjects who report wash hand/wear mask often/always during the f/u period - according to Q4

aftervisit <- read.csv(paste(dir, "adherence_m.csv", sep=""))
aftervisit <- merge(aftervisit,housechar[1:2],by="hhID",all.x=TRUE)
aftervisit$analyzed <- hculture$analyzed
aftervisit <- aftervisit[aftervisit$analyzed==1,]

# for only 36 hours
aftervisit <- aftervisit[aftervisit$hhID %in% hhID36,]
out <-function(value){
return( round(value *100,0))
}

appt10 <- matrix(rep(NA,98),ncol=14,byrow=FALSE)

# For control group
index_c <- aftervisit[aftervisit$intervention==1&aftervisit$member==0,]
appt10[c(3,6),1] <- c(out(dim(index_c[!is.na(index_c$washhand)&index_c$washhand<=2,])[1]/dim(index_c)[1]),
                      out(dim(index_c[!is.na(index_c$mask)&index_c$mask<=2,])[1]/dim(index_c)[1]))
contact_c <- aftervisit[aftervisit$intervention==1&aftervisit$member>0,]
appt10[c(3,6),2] <- c(out(dim(contact_c[!is.na(contact_c$washhand)&contact_c$washhand<=2,])[1]/dim(contact_c)[1]),
                      out(dim(contact_c[!is.na(contact_c$mask)&contact_c$mask<=2,])[1]/dim(contact_c)[1]))

# For hand hygiene group
index_h <- aftervisit[aftervisit$intervention==3&aftervisit$member==0,]
appt10[c(3,6),3] <- c(out(dim(index_h[!is.na(index_h$washhand)&index_h$washhand<=2,])[1]/dim(index_h)[1]),
                      out(dim(index_h[!is.na(index_h$mask)&index_h$mask<=2,])[1]/dim(index_h)[1]))
contact_h <- aftervisit[aftervisit$intervention==3&aftervisit$member>0,]
appt10[c(3,6),6] <- c(out(dim(contact_h[!is.na(contact_h$washhand)&contact_h$washhand<=2,])[1]/dim(contact_h)[1]),
                      out(dim(contact_h[!is.na(contact_h$mask)&contact_h$mask<=2,])[1]/dim(contact_h)[1]))

# For mask+hh group
index_m <- aftervisit[aftervisit$intervention==4&aftervisit$member==0,]
appt10[c(3,6),9] <- c(out(dim(index_m[!is.na(index_m$washhand)&index_m$washhand<=2,])[1]/dim(index_m)[1]),
                      out(dim(index_m[!is.na(index_m$mask)&index_m$mask<=2,])[1]/dim(index_m)[1]))
contact_m <- aftervisit[aftervisit$intervention==4&aftervisit$member>0,]
appt10[c(3,6),12] <- c(out(dim(contact_m[!is.na(contact_m$washhand)&contact_m$washhand<=2,])[1]/dim(contact_m)[1]),
                       out(dim(contact_m[!is.na(contact_m$mask)&contact_m$mask<=2,])[1]/dim(contact_m)[1]))

# For liquid soap use
# For control group
index_c <- aftervisit[aftervisit$intervention==1&aftervisit$member==0,]
contact_c <- aftervisit[aftervisit$intervention==1&aftervisit$member>0,]
appt10[1,1:2] <- c(out(dim(index_c[!is.na(index_c$soap)&index_c$soap<=2,])[1]/dim(index_c)[1]),
                   out(dim(contact_c[!is.na(contact_c$soap)&contact_c$soap<=2,])[1]/dim(contact_c)[1]))

# For hand hygiene group
index_h <- aftervisit[aftervisit$intervention==3&aftervisit$member==0,]
contact_h <- aftervisit[aftervisit$intervention==3&aftervisit$member>0,]
appt10[1,c(3,6)] <- c(out(dim(index_h[!is.na(index_h$soap)&index_h$soap<=2,])[1]/dim(index_h)[1]),
                      out(dim(contact_h[!is.na(contact_h$soap)&contact_h$soap<=2,])[1]/dim(contact_h)[1]))

# For mask+hh group
index_m <- aftervisit[aftervisit$intervention==4&aftervisit$member==0,]
contact_m <- aftervisit[aftervisit$intervention==4&aftervisit$member>0,]
appt10[1,c(9,12)] <- c(out(dim(index_m[!is.na(index_m$soap)&index_m$soap<=2,])[1]/dim(index_m)[1]),
                       out(dim(contact_m[!is.na(contact_m$soap)&contact_m$soap<=2,])[1]/dim(contact_m)[1]))

# For handrub use
# For control group
index_c <- aftervisit[aftervisit$intervention==1&aftervisit$member==0,]
contact_c <- aftervisit[aftervisit$intervention==1&aftervisit$member>0,]
appt10[2,1:2] <- c(out(dim(index_c[!is.na(index_c$handrub)&index_c$handrub<=2,])[1]/dim(index_c)[1]),
                   out(dim(contact_c[!is.na(contact_c$handrub)&contact_c$handrub<=2,])[1]/dim(contact_c)[1]))

# For hand hygiene group
index_h <- aftervisit[aftervisit$intervention==3&aftervisit$member==0,]
contact_h <- aftervisit[aftervisit$intervention==3&aftervisit$member>0,]
appt10[2,c(3,6)] <- c(out(dim(index_h[!is.na(index_h$handrub)&index_h$handrub<=2,])[1]/dim(index_h)[1]),
                      out(dim(contact_h[!is.na(contact_h$handrub)&contact_h$handrub<=2,])[1]/dim(contact_h)[1]))

# For mask+hh group
index_m <- aftervisit[aftervisit$intervention==4&aftervisit$member==0,]
contact_m <- aftervisit[aftervisit$intervention==4&aftervisit$member>0,]
appt10[2,c(9,12)] <- c(out(dim(index_m[!is.na(index_m$handrub)&index_m$handrub<=2,])[1]/dim(index_m)[1]),
                       out(dim(contact_m[!is.na(contact_m$handrub)&contact_m$handrub<=2,])[1]/dim(contact_m)[1]))


# -------------------- Average soap/handrub/mask used (median+IQR) --------------------------------###

soap <- read.csv(paste(dir, "adherence_h.csv", sep="")) 
handrub <- mask <- aftervisit

soap$analyzed <- hculture$analyzed[hculture$member==0]
soap <- soap[soap$analyzed==1,]
soap <- soap[soap$hhID %in% hhID36,]

soap <- merge(soap,housechar[1:2],by="hhID",all.x=TRUE)
soap$usage <- soap$soap_given-soap$soap_remain
handrub$usage <- handrub$smallgel_given-handrub$smallgel_remain
mask$usage <- mask$given_mask-mask$remain_mask


# For hand hygiene group
soap_h <- soap[soap$intervention==3,]
appt10[4,3:5] <- round(quantile(soap_h$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)

handrub_index_h <- handrub[handrub$intervention==3&handrub$member==0,]
appt10[5,3:5] <- round(quantile(handrub_index_h$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)
handrub_contact_h <- handrub[handrub$intervention==3&handrub$member>0,]
appt10[5,6:8] <- round(quantile(handrub_contact_h$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)

# For mask+hh group
soap_m <- soap[soap$intervention==4,]
appt10[4,9:11] <- round(quantile(soap_m$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)

handrub_index_m <- handrub[handrub$intervention==4&handrub$member==0,]
appt10[5,9:11] <- round(quantile(handrub_index_m$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)
handrub_contact_m <- handrub[handrub$intervention==4&handrub$member>0,]
appt10[5,12:14] <- round(quantile(handrub_contact_m$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)

mask_index_m <- mask[mask$intervention==4&mask$member==0,]
appt10[7,9:11] <- round(quantile(mask_index_m$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)
mask_contact_m <- mask[mask$intervention==4&mask$member>0,]
appt10[7,12:14] <- round(quantile(mask_contact_m$usage,probs=c(0.5,0.25,0.75),na.rm=TRUE),1)

colnames(appt10) <- c("C_index","C_contact","H_index","H_i_IQRlow","H_i_IQRup","H_contact","H_c_IQRlow","H_c_IQRup",
                       "M_index","M_i_IQRlow","M_i_IQRup","M_contact","M_c_IQRlow","M_c_IQRup")
rownames(appt10) <- c("liquid soap","hand rub","good hh","median soap","median rub","mask","median mask")
appt10

# End of script

