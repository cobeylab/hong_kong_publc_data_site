#
# R syntax to reproduce information for Table 2 from:
#
# Ng S, Ip DKM, Fang VJ, et al.
# The effect of age and recent influenza vaccination history on the immunogenicity and efficacy 
# of 2009-10 seasonal trivalent inactivated influenza vaccination in children
# PLoS ONE 2013 (in press).
#
# Last updated by Ng S, Fang VJ and Cowling BJ.
# March 8, 2013
#

source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_PLoSONE_scripts/dataframe.r")

dat$hxgp <-NA
dat$hxgp[dat$vac0708==1 & dat$vac0809==1] <-4
dat$hxgp[dat$vac0708==0 & dat$vac0809==1] <-3
dat$hxgp[dat$vac0708==1 & dat$vac0809==0] <-2
dat$hxgp[dat$vac0708==0 & dat$vac0809==0] <-1

tab <- matrix(rep(NA, 10*7), nrow=7, dimnames=list(c("6-8", "male", "chron", NA,"9-17","male", "chron"),c("n", "%","n", "%","n", "%","n", "%","n", "%")))

for (i in 1:4) {
                tab[1,(i-1)*2+1]<- dim(dat[dat$age%in%6:8 & dat$vac0910.pre==1 & dat$hxgp%in%i,])[1]
                tab[2,(i-1)*2+1]<- dim(dat[dat$male%in%1 & dat$age%in%6:8 & dat$vac0910.pre==1 & dat$hxgp%in%i,])[1]
                tab[3,(i-1)*2+1] <-dim(dat[dat$chron2%in%1 & dat$age%in%6:8 & dat$vac0910.pre==1 & dat$hxgp%in%i,])[1]
                tab[2:3,(i-1)*2+2]  <- round( tab[2:3,(i-1)*2+1] / tab[1,(i-1)*2+1] ,2)
                tab[1,9] <-dim( dat[dat$age%in%6:8 & dat$vac0910.pre==1 & is.na(dat$hxgp),])[1]
                tab[1,10]  <- round( tab[1,9] / 170 ,2)
                
                tab[5,(i-1)*2+1]<- dim(dat[dat$age%in%9:17 & dat$vac0910.pre==1 & dat$hxgp%in%i,])[1]
                tab[6,(i-1)*2+1]<- dim(dat[dat$male%in%1 & dat$age%in%9:17 & dat$vac0910.pre==1 & dat$hxgp%in%i,])[1]
                tab[7,(i-1)*2+1] <-dim(dat[dat$chron2%in%1 & dat$age%in%9:17 & dat$vac0910.pre==1 & dat$hxgp%in%i,])[1]
                tab[6:7,(i-1)*2+2]  <- round( tab[6:7,(i-1)*2+1] / tab[5,(i-1)*2+1] ,2)
                tab[5,9] <-dim(dat[dat$age%in%9:17 & dat$vac0910.pre==1 & is.na(dat$hxgp),])[1]
                tab[5,10] <-round(tab[5,9]/309,2)
                }

tab

#
# End of script.
#
