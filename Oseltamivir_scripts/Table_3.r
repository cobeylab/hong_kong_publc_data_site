#
# R syntax to reproduce information for Table 3 from:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

source("http://www.hku.hk/bcowling/influenza/Oseltamivir_scripts/dataframe.r")  
require(chron)
require("survival")

tab3 <- matrix(rep(NA,140),ncol=10)
tab3[c(3,7,8,10,13),c(2,7)] <- 1

# oseltamivir start day
tamiflu <- merge(av[av$member==0&av$av=="tamiflu",],hchar[c("hhID","clinic_date","clinic_day")],by="hhID",all.x=TRUE)
tamiflu$onset.tami <- dates(as.character(tamiflu$dayfrom),format="d/m/y")-dates(as.character(tamiflu$clinic_date),format="d/m/y")+tamiflu$clinic_day
tamiflu$onset.tami[tamiflu$onset.tami<2] <- 1
tamiflu$onset.tami[tamiflu$onset.tami>=2] <- 2


# add other variables
vshed$exact <- 3
vshed$timeL[vshed$timeL==0] <- 0.01
duration <- merge(vshed,tamiflu[c("hhID","onset.tami")],all.x=TRUE)
duration$onset.tami[is.na(duration$onset.tami)] <- 0
duration <- merge(duration,clinic7[c("hhID","agegp","male","vaccine","b_m0_score","flu.type")],all.x=TRUE)

# Regression
d07 <- duration[duration$hhID<8000,]
d08 <- duration[duration$hhID>8000,]

survtest07 <- survreg(Surv(time=timeL,time2=timeR, event=exact, type="interval") ~ factor(onset.tami)+factor(as.character(agegp))
                     +male+vaccine+b_m0_score+factor(flu.type), data = d07, dist="weibull")
output <- data.frame(coef=coef(survtest07), se=sqrt(diag(vcov(survtest07)))[-10])[-1,]
output$AF <- exp(output$coef)
output$lower.CI <- exp(output$coef - 1.96*output$se)
output$upper.CI <- exp(output$coef + 1.96*output$se)
output$p.value <- 2*pnorm(-1*abs(output$coef/output$se))
tab3[c(1,2,4:6,9,11,12,14),2:5] <- as.matrix(round(output, 2)[,-(1:2)])

tab3[,1] <- c(table(d07$onset.tami)[c(2,3,1)],table(d07$agegp),table(d07$male),table(d07$vaccine),NA,table(d07$flu.type))

##

survtest08 <- survreg(Surv(time=timeL,time2=timeR, event=exact, type="interval") ~ factor(onset.tami)+factor(as.character(agegp))
                     +male+vaccine+b_m0_score+factor(flu.type), data = d08, dist="weibull")
output <- data.frame(coef=coef(survtest08), se=sqrt(diag(vcov(survtest08)))[-10])[-1,]
output$AF <- exp(output$coef)
output$lower.CI <- exp(output$coef - 1.96*output$se)
output$upper.CI <- exp(output$coef + 1.96*output$se)
output$p.value <- 2*pnorm(-1*abs(output$coef/output$se))
tab3[c(1,2,4:6,9,11,12,14),7:10] <- as.matrix(round(output, 2)[,-(1:2)])

tab3[,6] <- c(table(d08$onset.tami)[c(2,3,1)],table(d08$agegp),table(d08$male),table(d08$vaccine),NA,table(d08$flu.type))

colnames(tab3) <- c("2007 (n)","AF","95%CI_low","95%CI_up","p-value","2008 (n)","AF","95%CI_low","95%CI_up","p-value")
rownames(tab3) <- c("Oseltamivir<2d","Oseltamivir>=2d","No antiviral","Age<5","Age6-12","Age13-17","Age18+",
                    "Female","Male","Not vaccine","Vaccine","Baseline score","FluA","FluB")

tab3

#
# End of script
#

