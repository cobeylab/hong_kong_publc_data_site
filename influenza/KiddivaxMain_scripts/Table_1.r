#
# R syntax to reproduce Table 1 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

dir <- "http://sph.hku.hk/data/KiddivaxMainV1/"

random <- read.csv(paste(dir, "randomcode.csv", sep=""))
demog <- read.csv(paste(dir, "demographic.csv", sep=""))
hchar <- read.csv(paste(dir, "housechar.csv", sep=""))

tab <- matrix(rep(NA,8*5), ncol=5, dimnames=list(c("n","Female","age6-8","age9-11","age12-17","vac0809","num member","flatsize"),
                                                 c("Missing","TIV","%/(sd)","placebo", "%/(sd)")))

demog <- merge(demog,random,all.x=T); index <- demog

# Number of subjects in each intervention group
tab[1,c(2,4)] <- rev(table(index$intervention))
# sex
tab[2,1] <- sum(is.na(index$male))
tab[2,c(2,4)] <- rev(table(index$intervention[index$male==0]))
tab[2,c(3,5)] <- round(prop.table(table(index$male,index$intervention),2)[1,2:1],2)
# agegp
index$agegp <- cut(index$age,c(0,8,11,17))
tab[3,1] <- sum(is.na(index$age))
tab[3:5,c(2,4)] <- table(index$agegp,index$intervention)[,2:1]
tab[3:5,c(3,5)] <- round(prop.table(table(index$agegp,index$intervention),2)[,2:1],2)
# vax history
tab[6,1] <- sum(is.na(index$vac0809))
tab[6,c(2,4)] <- rev(table(index$intervention[index$vac0809==1]))
tab[6,c(3,5)] <- round(prop.table(table(index$vac0809,index$intervention),2)[2,2:1],2)

### household characteristics
hchar <- merge(hchar,random,all.x=T)
# family size
tab[7,1] <- 0
tab[7,2:3] <- round(c(mean(hchar$familysize[hchar$intervention=="TIV"]),sd(hchar$familysize[hchar$intervention=="TIV"])),1)
tab[7,4:5] <- round(c(mean(hchar$familysize[hchar$intervention=="placebo"]),sd(hchar$familysize[hchar$intervention=="placebo"])),1)
# flat size
tab[8,1] <- sum(is.na(hchar$house_size))
tab[8,2:3] <- round(c(mean(hchar$house_size[hchar$intervention=="TIV"],na.rm=T),sd(hchar$house_size[hchar$intervention=="TIV"],na.rm=T)),0)
tab[8,4:5] <- round(c(mean(hchar$house_size[hchar$intervention=="placebo"]),sd(hchar$house_size[hchar$intervention=="placebo"])),0)

tab

#
# End of script.
#


