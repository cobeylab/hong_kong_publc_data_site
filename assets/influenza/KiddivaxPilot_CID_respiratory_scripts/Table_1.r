#
# R syntax to reproduce information for Table 1 from:
#
# Cowling BJ, Fang VJ, Nishiura H, et al.
# Increased risk of noninfluenza respiratory virus infections
# sssociated with receipt of inactivated influenza vaccine
# CID, 2012.
#
# Last updated by Fang VJ and Cowling BJ.
# January 24, 2013
#

require(chron)
dir <- "http://sph.hku.hk/data/KiddivaxPilotV2/"
sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
hchar <- read.csv(paste(dir, "housechar_h.csv", sep=""))

sero <- merge(sero[sero$member==0&sero$hhID!=9200&sero$hhID!=9219,],random,all.x=T)
sero <- merge(sero,demog[demog$member==0,c("hhID","male","age")],all.x=T)
sero <- merge(sero,hchar[c("hhID","familysize")],all.x=T)

table(sero$intervention)

tab <- matrix(NA,ncol=4,nrow=6)
colnames(tab) <- c("TIV","%","Placebo","%"); rownames(tab) <- c("Age:6-8","9-11","12-15","Female","Follow-up days","Mean household size")

tab[1:3,1] <- table(cut(sero$age[sero$intervention=="TIV"],c(0,8,11,15)))
tab[1:3,3] <- table(cut(sero$age[sero$intervention=="placebo"],c(0,8,11,15)))
tab[4,1] <- sum(sero$male[sero$intervention=="TIV"]==0)
tab[4,3] <- sum(sero$male[sero$intervention=="placebo"]==0)
tab[1:4,2] <- round(tab[1:4,1]/69,2);  tab[1:4,4] <- round(tab[1:4,3]/46,2)

sero$day0 <- as.numeric(dates(as.character(sero$start),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$day1 <- as.numeric(dates(as.character(sero$end),format="d/m/y")-dates("1/1/2009",format="d/m/y"))
sero$len <- sero$day1-sero$day0

tab[5,1] <- round(median(sero$len[sero$intervention=="TIV"])); tab[5,3] <- round(median(sero$len[sero$intervention=="placebo"]))
tab[6,1] <- round(mean(sero$familysize[sero$intervention=="TIV"]),1); tab[6,3] <- round(mean(sero$familysize[sero$intervention=="placebo"]),1)

tab


# End of script

