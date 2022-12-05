#
# R syntax to reproduce information for Table 2 from:
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

tab2 <- matrix(rep(NA,405),ncol=15)
tab2[c(4,8,9,11,14,16,18,20,22,24,26),c(2,7,12)] <- 1

# oseltamivir start day

tamiflu <- merge(av[av$member==0&av$av=="tamiflu",],hchar[c("hhID","clinic_date","clinic_day")],by="hhID",all.x=TRUE)
tamiflu$onset.tami <- dates(as.character(tamiflu$dayfrom),format="d/m/y")-dates(as.character(tamiflu$clinic_date),format="d/m/y")+tamiflu$clinic_day
tamiflu$onset.tami[tamiflu$onset.tami>=2] <- 2
tamiflu$onset.tami <- tamiflu$onset.tami+1

# all symptoms (n=384)
allsymp <- merge(output.all,tamiflu[c(1,8)],all.x=T)
allsymp$onset.tami[is.na(allsymp$onset.tami)] <- 0
allsymp$year <- 1*(allsymp$hhID>8000)

# AFT-model

survtest <- survreg(Surv(time=all_time,all_event) ~ factor(onset.tami)+factor(as.character(agegp))+male+vaccine+b_m0_score+factor(flu.type)
                          +chronic_disease+antibiotics+antipyretic+antihistamine+steroid+year, data = allsymp, dist="weibull")
output <- data.frame(coef=coef(survtest), se=sqrt(diag(vcov(survtest)))[-18])[-1,]
output$AF <- exp(output$coef)
output$lower.CI <- exp(output$coef - 1.96*output$se)
output$upper.CI <- exp(output$coef + 1.96*output$se)
output$p.value <- 2*pnorm(-1*abs(output$coef/output$se))
tab2[c(1:3,5:7,10,12,13,15,17,19,21,23,25,27),2:5] <- as.matrix(round(output, 2)[,-(1:2)])

tab2[,1] <- c(table(allsymp$onset.tami)[c(2:4,1)], table(allsymp$agegp), 
       table(allsymp$male), table(allsymp$vaccine), NA, table(allsymp$flu.type), table(allsymp$chronic_disease),  
       table(allsymp$antibiotics), table(allsymp$antipyretic), table(allsymp$antihistamine), table(allsymp$steroid), 
       table(allsymp$year))


# fever  (n=321)
feversymp <- merge(output.fever,tamiflu[c(1,8)],all.x=T)
feversymp$onset.tami[is.na(feversymp$onset.tami)] <- 0
feversymp$year <- 1*(feversymp$hhID>8000)

# AFT-model

survtest.f <- survreg(Surv(time=f_time,f_event) ~ factor(onset.tami)+factor(as.character(agegp))+male+vaccine+b_m0_score+factor(flu.type)
                          +chronic_disease+antibiotics+antipyretic+antihistamine+steroid+year, data = feversymp, dist="weibull")
output <- data.frame(coef=coef(survtest.f), se=sqrt(diag(vcov(survtest.f)))[-18])[-1,]
output$AF <- exp(output$coef)
output$lower.CI <- exp(output$coef - 1.96*output$se)
output$upper.CI <- exp(output$coef + 1.96*output$se)
output$p.value <- 2*pnorm(-1*abs(output$coef/output$se))
tab2[c(1:3,5:7,10,12,13,15,17,19,21,23,25,27),7:10] <- as.matrix(round(output, 2)[,-(1:2)])

tab2[,6] <- c(table(feversymp$onset.tami)[c(2:4,1)], table(feversymp$agegp), 
       table(feversymp$male), table(feversymp$vaccine), NA, table(feversymp$flu.type), table(feversymp$chronic_disease),  
       table(feversymp$antibiotics), table(feversymp$antipyretic), table(feversymp$antihistamine), table(feversymp$steroid), 
       table(feversymp$year))


# respiratory symptoms (n=380)
rsymp <- merge(output.rs,tamiflu[c(1,8)],all.x=T)
rsymp$onset.tami[is.na(rsymp$onset.tami)] <- 0
rsymp$year <- 1*(rsymp$hhID>8000)

# AFT-model

survtest.rs <- survreg(Surv(time=r_time,r_event) ~ factor(onset.tami)+factor(as.character(agegp))+male+vaccine+b_m0_score+factor(flu.type)
                          +chronic_disease+antibiotics+antipyretic+antihistamine+steroid+year, data = rsymp, dist="weibull")
output <- data.frame(coef=coef(survtest.rs), se=sqrt(diag(vcov(survtest.rs)))[-18])[-1,]
output$AF <- exp(output$coef)
output$lower.CI <- exp(output$coef - 1.96*output$se)
output$upper.CI <- exp(output$coef + 1.96*output$se)
output$p.value <- 2*pnorm(-1*abs(output$coef/output$se))
tab2[c(1:3,5:7,10,12,13,15,17,19,21,23,25,27),12:15] <- as.matrix(round(output, 2)[,-(1:2)])

tab2[,11] <- c(table(rsymp$onset.tami)[c(2:4,1)], table(rsymp$agegp), 
       table(rsymp$male), table(rsymp$vaccine), NA, table(rsymp$flu.type), table(rsymp$chronic_disease),  
       table(rsymp$antibiotics), table(rsymp$antipyretic), table(rsymp$antihistamine), table(rsymp$steroid), 
       table(rsymp$year))

colnames(tab2) <- c("All(n)","AF","95%CI_low","95%CI_up","p-value","Fever(n)","AF","95%CI_low","95%CI_up","p-value",
                    "Resp_symp(n)","AF","95%CI_low","95%CI_up","p-value")
rownames(tab2) <- c("Oseltamivir<1d","Oseltamivir1-2d","Oseltamivir>2d","No antiviral","Age<5","Age6-12","Age13-17","Age18+",
                    "Female","Male","No vaccine","Vaccine","Basline score","FluA","FluB","No chronic disease","Chronic disease",
                    "No antibiotic","Antibiotic","No antipyretic","Antipyretic","No antihistamine","Antihistamine",
                    "No steroid","Steroid","Year2007","Year2008")
tab2

#
# End of script
#











