#
# R syntax to reproduce information for Table 1 from:
#
# Klick B, Nishiura H, Ng S, et al.
# Transmissibility of seasonal and pandemic influenza in a cohort
# of households in Hong Kong in 2009
# Epidemiology. 2011 Nov;22(6):793-6.
#
# Last updated by Klick B, Fang VJ and Cowling BJ.
# December 22, 2011
#

# NOTE: The numbers of sH3N2 infection in age group 0-14 and 15-59 should be 14 and 13, while in the paper the two numbers were switched (typo).

dir <- "http://sph.hku.hk/data/KiddivaxPilotV1/"

demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
random <- read.csv(paste(dir, "randomcode_h.csv", sep=""))

demog <- demog[demog$hhID!=9108&demog$hhID!=9172,]

sero$ph1.post.titer[(sero$hhID==9186|sero$hhID==9194|sero$hhID==9213|sero$hhID==9218)&sero$member==0] <- 20

#convert ph1 titer to numeric & NA for Positive PositiveOne & 5 for Negative
sero$ph1.mid.titer[sero$ph1.mid.titer=="Negative"] <-5
sero$ph1.post.titer[sero$ph1.post.titer=="Negative" | sero$ph1.post.titer=="PositiveOne"] <-5

sero$ph1 <- 1*(as.numeric(as.character(sero$ph1.post.titer))/as.numeric(as.character(sero$ph1.mid.titer))>=4&as.numeric(as.character(sero$ph1.post.titer))>=40)
sero$ph1[is.na(sero$ph1)&sero$ph1.post.titer=="Positive"] <- 1
sero$ph1[is.na(sero$ph1)&sero$ph1.post.titer==5] <- 0
sero$sh1 <- 1*(sero$sh1.posts/sero$sh1.mids>=4&sero$sh1.posts>=40)
sero$sh3 <- 1*(sero$sh3.posts/sero$sh3.mids>=4&sero$sh3.posts>=40)

## classifly the 12 cross antibody titer rise.
sero$sh1[sero$hhID==9101&sero$member==0] <- 0
sero$sh1[sero$hhID==9110&sero$member==2] <- 0
sero$sh3[sero$hhID==9122&sero$member==4] <- 0
sero$sh1[sero$hhID==9128&sero$member==3] <- 0
sero$sh1[sero$hhID==9161&sero$member==1] <- 0
sero$sh1[sero$hhID==9167&sero$member==3] <- 0
sero$sh3[sero$hhID==9167&sero$member==3] <- 0
sero$sh1[sero$hhID==9193&sero$member==0] <- 0
sero$sh3[sero$hhID==9193&sero$member==0] <- 0
sero$sh1[sero$hhID==9195&sero$member==3] <- 0
sero$sh1[sero$hhID==9203&sero$member==3] <- 0
sero$sh3[sero$hhID==9206&sero$member==3] <- 0
sero$ph1[sero$hhID==9208&sero$member==1] <- 0
sero$sh1[sero$hhID==9211&sero$member==0] <- 0
sero$sh3[sero$hhID==9214&sero$member==3] <- 0
##

demog <- merge(demog,random,all.x=T)
all2 <- merge(demog,sero[c("hhID","member","ph1","sh1","sh3")],by=c("hhID","member"),all.x=T)

#
tab <-matrix(rep(NA,6*7), ncol=7,
              dimnames=list(c("Age:0-14","15-59","60+","Male","Female","vax0809"),c("Num","pH1","%","sH1","%","sH3","%")))

# age
all2$agegp <- cut(all2$age,c(-0.1,14,59,100))
tab[1:3,1] <- table(all2$agegp)
tab[1:3,2] <- table(all2$agegp[all2$ph1==1&!is.na(all2$ph1)])
tab[1:3,4] <- table(all2$agegp[all2$sh1==1&!is.na(all2$sh1)])
tab[1:3,6] <- table(all2$agegp[all2$sh3==1&!is.na(all2$sh3)])
tmp <- table(all2$agegp[!is.na(all2$ph1)]); tab[1:3,3] <- round(tab[1:3,2]/tmp,2)
tmp <- table(all2$agegp[!is.na(all2$sh1)]); tab[1:3,5] <- round(tab[1:3,4]/tmp,2)
tmp <- table(all2$agegp[!is.na(all2$sh3)]); tab[1:3,7] <- round(tab[1:3,6]/tmp,2)

# sex
tab[4:5,1] <- rev(table(all2$male))
tab[4:5,2] <- rev(table(all2$male[all2$ph1==1&!is.na(all2$ph1)]))
tab[4:5,4] <- rev(table(all2$male[all2$sh1==1&!is.na(all2$sh1)]))
tab[4:5,6] <- rev(table(all2$male[all2$sh3==1&!is.na(all2$sh3)]))
tmp <- rev(table(all2$male[!is.na(all2$ph1)])); tab[4:5,3] <- round(tab[4:5,2]/tmp,2)
tmp <- rev(table(all2$male[!is.na(all2$sh1)])); tab[4:5,5] <- round(tab[4:5,4]/tmp,2)
tmp <- rev(table(all2$male[!is.na(all2$sh3)])); tab[4:5,7] <- round(tab[4:5,6]/tmp,2)

# vac0809
all2$vaccine <- pmax(all2$vac09.pre,all2$vac09.mid,1*(all2$intervention=="TIV"&all2$member==0),na.rm=T)
tab[6,1] <- sum(all2$vaccine==1)
tab[6,2] <- sum(all2$vaccine[all2$ph1==1&!is.na(all2$ph1)]==1)
tab[6,4] <- sum(all2$vaccine[all2$sh1==1&!is.na(all2$sh1)]==1)
tab[6,6] <- sum(all2$vaccine[all2$sh3==1&!is.na(all2$sh3)]==1)
tmp <- sum(all2$vaccine[!is.na(all2$ph1)]==1); tab[6,3] <- round(tab[6,2]/tmp,2)
tmp <- sum(all2$vaccine[!is.na(all2$sh1)]==1); tab[6,5] <- round(tab[6,4]/tmp,2)
tmp <- sum(all2$vaccine[!is.na(all2$sh3)]==1); tab[6,7] <- round(tab[6,6]/tmp,2)

tab

# End of script.

