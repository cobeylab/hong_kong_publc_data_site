#
# R syntax to reproduce information for Table 1 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 11, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe.r")

hchar <- read.csv(paste(dir, "housechar_h.csv", sep=""))

tab1 <-matrix(rep(NA,16*4), ncol=4,
              dimnames=list(c("vaccinees", "male", "age6-8","age9-11","age12-15",
                            "vac08","contacts", "male", "age<15", "age16-45","age>45",
                            "vac08", "vac09", "HH", "no member", "flatsize"),
                            c("TIV","% /(sd)","placebo", "%")))

index.pre <- merge(index.pre,hchar,all.x=T)

index.pre$gp[index.pre$intervention=="TIV"] <-"1"
index.pre$gp[index.pre$intervention=="placebo"] <-"2"
member.pre$gp[member.pre$intervention=="TIV"] <-"1"
member.pre$gp[member.pre$intervention=="placebo"] <-"2"

# index
tab1[1,c(1,3)] <-table(index.pre$gp)
tab1[2,c(1,3)] <-table(index.pre$gp[index.pre$male==1])
tab1[c(3:5),c(1,3)] <-table(cut(index.pre$age,c(5,8,11,16)), index.pre$gp)
tab1[6,c(1,3)] <-table(index.pre$gp[index.pre$vac08==1])

# contacts
tab1[7,c(1,3)] <-table(member.pre$gp)
tab1[8,c(1,3)] <-table(member.pre$gp[member.pre$male==1])
tab1[c(9:11),c(1,3)] <-table(cut(member.pre$age,c(0,15,45,100)), member.pre$gp)
tab1[12,c(1,3)] <-table(member.pre$gp[member.pre$vac08==1])
tab1[13,c(1,3)] <-table(member.pre$gp[member.pre$vac09==1])

# Households
tab1[14,c(1,3)] <-table(index.pre$gp)
tab1[15,1] <-round(mean(index.pre$familysize[index.pre$gp==1],na.rm=T),1)
tab1[15,2] <-round(sd(index.pre$familysize[index.pre$gp==1],na.rm=T),1)
tab1[15,3] <-round(mean(index.pre$familysize[index.pre$gp==2],na.rm=T),1)
tab1[15,4] <-round(sd(index.pre$familysize[index.pre$gp==2],na.rm=T),1)
tab1[16,1] <-round(mean(index.pre$house_size[index.pre$gp==1],na.rm=T)/10.7639104167,1) #convert to sq meter
tab1[16,2] <-round(sd(index.pre$house_size[index.pre$gp==1],na.rm=T)/10.7639104167,1)  #convert to sq meter
tab1[16,3] <-round(mean(index.pre$house_size[index.pre$gp==2],na.rm=T)/10.7639104167,1)  #convert to sq meter
tab1[16,4] <-round(sd(index.pre$house_size[index.pre$gp==2],na.rm=T)/10.7639104167,1)   #convert to sq meter

# add %
tab1[c(2:6),2] <-round(tab1[c(2:6),1]/tab1[1,1]*100,0)
tab1[c(8:13),2] <-round(tab1[c(8:13),1]/tab1[7,1]*100,0)

tab1[c(2:6),4] <-round(tab1[c(2:6),3]/tab1[1,3]*100,0)
tab1[c(8:13),4] <-round(tab1[c(8:13),3]/tab1[7,3]*100,0)

tab1

#
# End of script.
#







