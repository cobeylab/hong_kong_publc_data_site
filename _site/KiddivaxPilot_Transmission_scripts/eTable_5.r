#
# R syntax to reproduce information for eTable 5 from:
#
# Klick B, Nishiura H, Ng S, et al.
# Transmissibility of seasonal and pandemic influenza in a cohort
# of households in Hong Kong in 2009
# Epidemiology. 2011 Nov;22(6):793-6.
#
# Last updated by Klick B, Fang VJ and Cowling BJ.
# January 3, 2012
#

# REMARK: Please note that each MCMC loop may take hours of running time, depending on the computer speed.

dir <- "http://sph.hku.hk/data/KiddivaxPilotV1/"
source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_Transmission_scripts/MCMC_function.r")

sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
code <- read.csv(paste(dir, "randomcode_h.csv", sep=""))

tab <- matrix(rep(NA,9*12),ncol=12,
              dimnames=list(c("All:pH1","sH1","sH3","<=1:160:pH1","sH1","sH3","<=1:40:pH1","sH1","sH3"),
                            c("CPI-child","CI_low","CI_up","CPI-adult","CI_low","CI_up","SAP-child","CI_low","CI_up","SAP-adult","CI_low","CI_up")))

sero$ph1.post.titer[(sero$hhID==9186|sero$hhID==9194|sero$hhID==9213|sero$hhID==9218)&sero$member==0] <- 20

#convert ph1 titer to numeric & NA for Positive PositiveOne & 5 for Negative
sero$ph1.mid.titer[sero$ph1.mid.titer=="Negative"] <-5
sero$ph1.post.titer[sero$ph1.post.titer=="Negative" | sero$ph1.post.titer=="PositiveOne"] <-5

sero$ph1s <- 1*(as.numeric(as.character(sero$ph1.post.titer))/as.numeric(as.character(sero$ph1.mid.titer))>=4&as.numeric(as.character(sero$ph1.post.titer))>=40)
sero$ph1s[is.na(sero$ph1s)&sero$ph1.post.titer=="Positive"] <- 1
sero$ph1s[is.na(sero$ph1s)&sero$ph1.post.titer==5] <- 0
sero$sh1s <- 1*(sero$sh1.posts/sero$sh1.mids>=4&sero$sh1.posts>=40)
sero$sh3s <- 1*(sero$sh3.posts/sero$sh3.mids>=4&sero$sh3.posts>=40)

c_data<-merge(x=sero,y=demog,by=c("hhID","member"),all.x=T,all.y=T)

c_data<-subset(c_data, hhID!=9108)
c_data<-subset(c_data, hhID!=9172)


###################################
#### We start MCMC for eTable5 ####
###################################

###############
### p H1N1 ####
###############

set.seed(35)

pH1_data<-c_data
pH1_data$d1_pmax<-pmax(pH1_data$ph1s,na.rm=TRUE)
pH1_data$d2_pmax<-pmax(pH1_data$ph1s,na.rm=FALSE)
pH1_data$ph1_miss<-ifelse(pH1_data$d1_pmax==1|!is.na(pH1_data$d2_pmax),0,1)
pH1_data$ph1_miss<-ifelse(is.na(pH1_data$ph1_miss),1,pH1_data$ph1_miss)
pH1_data$ph1_inf<-ifelse(pH1_data$d1_pmax==1|!is.na(pH1_data$d2_pmax),pH1_data$d1_pmax,NA)
pH1_data$child<-ifelse(pH1_data$age<=14,1,0)
pH1_data$adult<-ifelse(pH1_data$age>14,1,0)
pH1_data$ch_inf<-pH1_data$ph1_inf*pH1_data$child
pH1_data$ad_inf<-pH1_data$ph1_inf*pH1_data$adult

imp.inf<-rbinom(length(pH1_data$d1_pmax),1,.1)

pH1_data$imput_inf<-ifelse(pH1_data$ph1_miss==1,imp.inf,pH1_data$ph1_inf)
pH1_data$imp_ch_inf<-pH1_data$imput_inf*pH1_data$child
pH1_data$imp_ad_inf<-pH1_data$imput_inf*pH1_data$adult

t.ch<-as.vector(by(pH1_data$child, pH1_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(pH1_data$imp_ch_inf, pH1_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(pH1_data$adult, pH1_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(pH1_data$imp_ad_inf, pH1_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(pH1_data$child, pH1_data$hhID, length))

pH1_data$n_ch_inf<-rep(inf.ch, time=len_fam)
pH1_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
mcmc.num <- 15504
ph1_sum<-mcmc.fun.1(num=mcmc.num,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="ph1", mcmc.data=pH1_data, inf.data=inf)
write.csv(ph1_sum, "ph1_sum_modif.csv")

###############
### s H1N1 ####
###############

set.seed(35)

sH1_data<-c_data
sH1_data$d1_pmax<-pmax(sH1_data$sh1s,na.rm=TRUE)
sH1_data$d2_pmax<-pmax(sH1_data$sh1s,na.rm=FALSE)
sH1_data$sh1_miss<-ifelse(sH1_data$d1_pmax==1|!is.na(sH1_data$d2_pmax),0,1)
sH1_data$sh1_miss<-ifelse(is.na(sH1_data$sh1_miss),1,sH1_data$sh1_miss)
sH1_data$sh1_inf<-ifelse(sH1_data$d1_pmax==1|!is.na(sH1_data$d2_pmax),sH1_data$d1_pmax,NA)
sH1_data$child<-ifelse(sH1_data$age<=14,1,0)
sH1_data$adult<-ifelse(sH1_data$age>14,1,0)
sH1_data$ch_inf<-sH1_data$sh1_inf*sH1_data$child
sH1_data$ad_inf<-sH1_data$sh1_inf*sH1_data$adult

imp.inf<-rbinom(length(sH1_data$d1_pmax),1,.1)

sH1_data$imput_inf<-ifelse(sH1_data$sh1_miss==1,imp.inf,sH1_data$sh1_inf)
sH1_data$imp_ch_inf<-sH1_data$imput_inf*sH1_data$child
sH1_data$imp_ad_inf<-sH1_data$imput_inf*sH1_data$adult

t.ch<-as.vector(by(sH1_data$child, sH1_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(sH1_data$imp_ch_inf, sH1_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(sH1_data$adult, sH1_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(sH1_data$imp_ad_inf, sH1_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(sH1_data$child, sH1_data$hhID, length))

sH1_data$n_ch_inf<-rep(inf.ch, time=len_fam)
sH1_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
sh1_sum<-mcmc.fun.2(num=15512,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="sh1", mcmc.data=sH1_data, inf.data=inf)
write.csv(sh1_sum, "sh1_sum_modif_beta1.csv")

################
#### s H3N2 ####
################

set.seed(35)

sH3_data<-c_data
sH3_data$d1_pmax<-pmax(sH3_data$sh3s,na.rm=TRUE)
sH3_data$d2_pmax<-pmax(sH3_data$sh3s,na.rm=FALSE)
sH3_data$sh3_miss<-ifelse(sH3_data$d1_pmax==1|!is.na(sH3_data$d2_pmax),0,1)
sH3_data$sh3_miss<-ifelse(is.na(sH3_data$sh3_miss),1,sH3_data$sh3_miss)
sH3_data$sh3_inf<-ifelse(sH3_data$d1_pmax==1|!is.na(sH3_data$d2_pmax),sH3_data$d1_pmax,NA)
sH3_data$child<-ifelse(sH3_data$age<=14,1,0)
sH3_data$adult<-ifelse(sH3_data$age>14,1,0)
sH3_data$ch_inf<-sH3_data$sh3_inf*sH3_data$child
sH3_data$ad_inf<-sH3_data$sh3_inf*sH3_data$adult

imp.inf<-rbinom(length(sH3_data$d1_pmax),1,.1)

sH3_data$imput_inf<-ifelse(sH3_data$sh3_miss==1,imp.inf,sH3_data$sh3_inf)
sH3_data$imp_ch_inf<-sH3_data$imput_inf*sH3_data$child
sH3_data$imp_ad_inf<-sH3_data$imput_inf*sH3_data$adult

t.ch<-as.vector(by(sH3_data$child, sH3_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(sH3_data$imp_ch_inf, sH3_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(sH3_data$adult, sH3_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(sH3_data$imp_ad_inf, sH3_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(sH3_data$child, sH3_data$hhID, length))

sH3_data$n_ch_inf<-rep(inf.ch, time=len_fam)
sH3_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
sh3_sum<-mcmc.fun.1(num=15504,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="sh3", mcmc.data=sH3_data, inf.data=inf)
write.csv(sh3_sum, "sh3_sum_modif_beta1.csv")

######################################
### p H1N1  -- equal or under 160 ####
######################################

set.seed(49)

pH1_data<-subset(c_data,as.numeric(as.character(ph1.mid.titer))<161)

pH1_data$d1_pmax<-pmax(pH1_data$ph1s,na.rm=TRUE)
pH1_data$d2_pmax<-pmax(pH1_data$ph1s,na.rm=FALSE)
pH1_data$ph1_miss<-ifelse(pH1_data$d1_pmax==1|!is.na(pH1_data$d2_pmax),0,1)
pH1_data$ph1_miss<-ifelse(is.na(pH1_data$ph1_miss),1,pH1_data$ph1_miss)
pH1_data$ph1_inf<-ifelse(pH1_data$d1_pmax==1|!is.na(pH1_data$d2_pmax),pH1_data$d1_pmax,NA)
pH1_data$child<-ifelse(pH1_data$age<=14,1,0)
pH1_data$adult<-ifelse(pH1_data$age>14,1,0)
pH1_data$ch_inf<-pH1_data$ph1_inf*pH1_data$child
pH1_data$ad_inf<-pH1_data$ph1_inf*pH1_data$adult

imp.inf<-rbinom(length(pH1_data$d1_pmax),1,.1)

pH1_data$imput_inf<-ifelse(pH1_data$ph1_miss==1,imp.inf,pH1_data$ph1_inf)
pH1_data$imp_ch_inf<-pH1_data$imput_inf*pH1_data$child
pH1_data$imp_ad_inf<-pH1_data$imput_inf*pH1_data$adult

t.ch<-as.vector(by(pH1_data$child, pH1_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(pH1_data$imp_ch_inf, pH1_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(pH1_data$adult, pH1_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(pH1_data$imp_ad_inf, pH1_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(pH1_data$child, pH1_data$hhID, length))

pH1_data$n_ch_inf<-rep(inf.ch, time=len_fam)
pH1_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
ph1_sum<-mcmc.fun.1(num=15509,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="ph1", mcmc.data=pH1_data, inf.data=inf)
write.csv(ph1_sum, "ph1_sum_modif_under160.csv")

#####################################
### s H1N1 -- equal or under 160 ####
#####################################

set.seed(35)

sH1_data<-subset(c_data,sh1.mids<161)

sH1_data$d1_pmax<-pmax(sH1_data$sh1s,na.rm=TRUE)
sH1_data$d2_pmax<-pmax(sH1_data$sh1s,na.rm=FALSE)
sH1_data$sh1_miss<-ifelse(sH1_data$d1_pmax==1|!is.na(sH1_data$d2_pmax),0,1)
sH1_data$sh1_miss<-ifelse(is.na(sH1_data$sh1_miss),1,sH1_data$sh1_miss)
sH1_data$sh1_inf<-ifelse(sH1_data$d1_pmax==1|!is.na(sH1_data$d2_pmax),sH1_data$d1_pmax,NA)
sH1_data$child<-ifelse(sH1_data$age<=14,1,0)
sH1_data$adult<-ifelse(sH1_data$age>14,1,0)
sH1_data$ch_inf<-sH1_data$sh1_inf*sH1_data$child
sH1_data$ad_inf<-sH1_data$sh1_inf*sH1_data$adult

imp.inf<-rbinom(length(sH1_data$d1_pmax),1,.1)

sH1_data$imput_inf<-ifelse(sH1_data$sh1_miss==1,imp.inf,sH1_data$sh1_inf)
sH1_data$imp_ch_inf<-sH1_data$imput_inf*sH1_data$child
sH1_data$imp_ad_inf<-sH1_data$imput_inf*sH1_data$adult

t.ch<-as.vector(by(sH1_data$child, sH1_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(sH1_data$imp_ch_inf, sH1_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(sH1_data$adult, sH1_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(sH1_data$imp_ad_inf, sH1_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(sH1_data$child, sH1_data$hhID, length))

sH1_data$n_ch_inf<-rep(inf.ch, time=len_fam)
sH1_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
sh1_sum<-mcmc.fun.2(num=15508,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="sh1", mcmc.data=sH1_data, inf.data=inf)
write.csv(sh1_sum, "sh1_sum_modif_beta1_under160.csv")

######################################
### s H3N2  -- equal or under 160 ####
######################################

set.seed(55)

sH3_data<-subset(c_data,sh3.mids<161)

sH3_data$d1_pmax<-pmax(sH3_data$sh3s,na.rm=TRUE)
sH3_data$d2_pmax<-pmax(sH3_data$sh3s,na.rm=FALSE)
sH3_data$sh3_miss<-ifelse(sH3_data$d1_pmax==1|!is.na(sH3_data$d2_pmax),0,1)
sH3_data$sh3_miss<-ifelse(is.na(sH3_data$sh3_miss),1,sH3_data$sh3_miss)
sH3_data$sh3_inf<-ifelse(sH3_data$d1_pmax==1|!is.na(sH3_data$d2_pmax),sH3_data$d1_pmax,NA)
sH3_data$child<-ifelse(sH3_data$age<=14,1,0)
sH3_data$adult<-ifelse(sH3_data$age>14,1,0)
sH3_data$ch_inf<-sH3_data$sh3_inf*sH3_data$child
sH3_data$ad_inf<-sH3_data$sh3_inf*sH3_data$adult

imp.inf<-rbinom(length(sH3_data$d1_pmax),1,.1)

sH3_data$imput_inf<-ifelse(sH3_data$sh3_miss==1,imp.inf,sH3_data$sh3_inf)
sH3_data$imp_ch_inf<-sH3_data$imput_inf*sH3_data$child
sH3_data$imp_ad_inf<-sH3_data$imput_inf*sH3_data$adult

t.ch<-as.vector(by(sH3_data$child, sH3_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(sH3_data$imp_ch_inf, sH3_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(sH3_data$adult, sH3_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(sH3_data$imp_ad_inf, sH3_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(sH3_data$child, sH3_data$hhID, length))

sH3_data$n_ch_inf<-rep(inf.ch, time=len_fam)
sH3_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
sh3_sum<-mcmc.fun.1(num=15507,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="sh3", mcmc.data=sH3_data, inf.data=inf)
write.csv(sh3_sum, "sh3_sum_modif_beta1_under160.csv")

####################################
### p H1N1 -- equal or under 40 ####
####################################

set.seed(48)

pH1_data<-subset(c_data,as.numeric(as.character(ph1.mid.titer))<41)

pH1_data$d1_pmax<-pmax(pH1_data$ph1s,na.rm=TRUE)
pH1_data$d2_pmax<-pmax(pH1_data$ph1s,na.rm=FALSE)
pH1_data$ph1_miss<-ifelse(pH1_data$d1_pmax==1|!is.na(pH1_data$d2_pmax),0,1)
pH1_data$ph1_miss<-ifelse(is.na(pH1_data$ph1_miss),1,pH1_data$ph1_miss)
pH1_data$ph1_inf<-ifelse(pH1_data$d1_pmax==1|!is.na(pH1_data$d2_pmax),pH1_data$d1_pmax,NA)
pH1_data$child<-ifelse(pH1_data$age<=14,1,0)
pH1_data$adult<-ifelse(pH1_data$age>14,1,0)
pH1_data$ch_inf<-pH1_data$ph1_inf*pH1_data$child
pH1_data$ad_inf<-pH1_data$ph1_inf*pH1_data$adult

imp.inf<-rbinom(length(pH1_data$d1_pmax),1,.1)

pH1_data$imput_inf<-ifelse(pH1_data$ph1_miss==1,imp.inf,pH1_data$ph1_inf)
pH1_data$imp_ch_inf<-pH1_data$imput_inf*pH1_data$child
pH1_data$imp_ad_inf<-pH1_data$imput_inf*pH1_data$adult

t.ch<-as.vector(by(pH1_data$child, pH1_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(pH1_data$imp_ch_inf, pH1_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(pH1_data$adult, pH1_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(pH1_data$imp_ad_inf, pH1_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(pH1_data$child, pH1_data$hhID, length))

pH1_data$n_ch_inf<-rep(inf.ch, time=len_fam)
pH1_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
ph1_sum<-mcmc.fun.1(num=15511,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="ph1", mcmc.data=pH1_data, inf.data=inf)
write.csv(ph1_sum, "ph1_sum_modif_under40.csv")

####################################
### s H1N1 -- equal or under 40 ####
####################################

set.seed(35)

sH1_data<-subset(c_data,sh1.mids<41)

sH1_data$d1_pmax<-pmax(sH1_data$sh1s,na.rm=TRUE)
sH1_data$d2_pmax<-pmax(sH1_data$sh1s,na.rm=FALSE)
sH1_data$sh1_miss<-ifelse(sH1_data$d1_pmax==1|!is.na(sH1_data$d2_pmax),0,1)
sH1_data$sh1_miss<-ifelse(is.na(sH1_data$sh1_miss),1,sH1_data$sh1_miss)
sH1_data$sh1_inf<-ifelse(sH1_data$d1_pmax==1|!is.na(sH1_data$d2_pmax),sH1_data$d1_pmax,NA)
sH1_data$child<-ifelse(sH1_data$age<=14,1,0)
sH1_data$adult<-ifelse(sH1_data$age>14,1,0)
sH1_data$ch_inf<-sH1_data$sh1_inf*sH1_data$child
sH1_data$ad_inf<-sH1_data$sh1_inf*sH1_data$adult

imp.inf<-rbinom(length(sH1_data$d1_pmax),1,.1)

sH1_data$imput_inf<-ifelse(sH1_data$sh1_miss==1,imp.inf,sH1_data$sh1_inf)
sH1_data$imp_ch_inf<-sH1_data$imput_inf*sH1_data$child
sH1_data$imp_ad_inf<-sH1_data$imput_inf*sH1_data$adult

t.ch<-as.vector(by(sH1_data$child, sH1_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(sH1_data$imp_ch_inf, sH1_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(sH1_data$adult, sH1_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(sH1_data$imp_ad_inf, sH1_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(sH1_data$child, sH1_data$hhID, length))

sH1_data$n_ch_inf<-rep(inf.ch, time=len_fam)
sH1_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
sh1_sum<-mcmc.fun.2(num=15508,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="sh1", mcmc.data=sH1_data, inf.data=inf)
write.csv(sh1_sum, "sh1_sum_modif_beta1_under40.csv")

#####################################
### s H3N2  -- equal or under 40 ####
#####################################

set.seed(35)

sH3_data<-subset(c_data,sh3.mids<41)

sH3_data$d1_pmax<-pmax(sH3_data$sh3s,na.rm=TRUE)
sH3_data$d2_pmax<-pmax(sH3_data$sh3s,na.rm=FALSE)
sH3_data$sh3_miss<-ifelse(sH3_data$d1_pmax==1|!is.na(sH3_data$d2_pmax),0,1)
sH3_data$sh3_miss<-ifelse(is.na(sH3_data$sh3_miss),1,sH3_data$sh3_miss)
sH3_data$sh3_inf<-ifelse(sH3_data$d1_pmax==1|!is.na(sH3_data$d2_pmax),sH3_data$d1_pmax,NA)
sH3_data$child<-ifelse(sH3_data$age<=16,1,0)
sH3_data$adult<-ifelse(sH3_data$age>16,1,0)
sH3_data$ch_inf<-sH3_data$sh3_inf*sH3_data$child
sH3_data$ad_inf<-sH3_data$sh3_inf*sH3_data$adult

imp.inf<-rbinom(length(sH3_data$d1_pmax),1,.1)

sH3_data$imput_inf<-ifelse(sH3_data$sh3_miss==1,imp.inf,sH3_data$sh3_inf)
sH3_data$imp_ch_inf<-sH3_data$imput_inf*sH3_data$child
sH3_data$imp_ad_inf<-sH3_data$imput_inf*sH3_data$adult

t.ch<-as.vector(by(sH3_data$child, sH3_data$hhID, sum, na.rm=TRUE))
inf.ch<-as.vector(by(sH3_data$imp_ch_inf, sH3_data$hhID, sum, na.rm=TRUE))
t.ad<-as.vector(by(sH3_data$adult, sH3_data$hhID, sum, na.rm=TRUE))
inf.ad<-as.vector(by(sH3_data$imp_ad_inf, sH3_data$hhID, sum, na.rm=TRUE))

len_fam<-as.vector(by(sH3_data$child, sH3_data$hhID, length))

sH3_data$n_ch_inf<-rep(inf.ch, time=len_fam)
sH3_data$n_ad_inf<-rep(inf.ad, time=len_fam)

inf<-data.frame(x1=inf.ch, x2=t.ch, x3=inf.ad, x4=t.ad)

# MCMC loop
sh3_sum<-mcmc.fun.1(num=15504,st.q_h.ch=3/4, st.q_h.ad=3/4, st.q_c.ch=7/8, st.q_c.ad=7/8, var.name="sh3", mcmc.data=sH3_data, inf.data=inf)
write.csv(sh3_sum, "sh3_sum_modif_beta1_under40.csv")

#################
#### eTable5 ####
#################

# ph1 everyone
d1<-read.csv("ph1_sum_modif.csv")
tab[1,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[1,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[1,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[1,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[1,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[1,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[1,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[1,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# sh1 everyone
d1<-read.csv("sh1_sum_modif_beta1.csv")
tab[2,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[2,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[2,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[2,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[2,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[2,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[2,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[2,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# sh3 everyone
d1<-read.csv("sh3_sum_modif_beta1.csv")
tab[3,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[3,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[3,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[3,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[3,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[3,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[3,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[3,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# ph1 under or equal 1:160
d1<-read.csv("ph1_sum_modif_under160.csv")
tab[4,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[4,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[4,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[4,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[4,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[4,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[4,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[4,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# sh1 under or equal 1:160
d1<-read.csv("sh1_sum_modif_beta1_under160.csv")
tab[5,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[5,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[5,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[5,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[5,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[5,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[5,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[5,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# sh3 under or equal 1:160
d1<-read.csv("sh3_sum_modif_beta1_under160.csv")
tab[6,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[6,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[6,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[6,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[6,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[6,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[6,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[6,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# ph1 under or equal 1:40
d1<-read.csv("ph1_sum_modif_under40.csv")
tab[7,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[7,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[7,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[7,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[7,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[7,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[7,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[7,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# sh1 under or equal 1:40
d1<-read.csv("sh1_sum_modif_beta1_under40.csv")
tab[8,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[8,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[8,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[8,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[8,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[8,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[8,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[8,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

# sh3 under or equal 1:40
d1<-read.csv("sh3_sum_modif_beta1_under40.csv")
tab[9,1] <- 1-mean(d1$q_c.ch.v[10000:15000]); tab[9,2:3] <- 1-sort(d1$q_c.ch.v[10000:15000])[c(4875,125)]
tab[9,4] <- 1-mean(d1$q_c.ad.v[10000:15000]); tab[9,5:6] <- 1-sort(d1$q_c.ad.v[10000:15000])[c(4875,125)]
tab[9,7] <- 1-mean(d1$q_h.ch.v[10000:15000]); tab[9,8:9] <- 1-sort(d1$q_h.ch.v[10000:15000])[c(4875,125)]
tab[9,10] <- 1-mean(d1$q_h.ad.v[10000:15000]); tab[9,11:12] <- 1-sort(d1$q_h.ad.v[10000:15000])[c(4875,125)]

tab <- round(tab,2)
tab

# End of script.
