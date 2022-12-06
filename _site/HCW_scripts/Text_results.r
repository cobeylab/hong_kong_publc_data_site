#
# R syntax to reproduce information in the manuscript text from:
#
# Zhou Y, Ng DMW, Seto WH, et al.
# Seroprevalence of antibody to pandemic influenza A (H1N1) 2009 among
# healthcare works after the first wave in Hong Kong.
# Journal of Hospital Infection, 2011.
#
# Last updated by Zhou Y and Fang VJ.
# September 20, 2011

dir <- "http://sph.hku.hk/data/"
data <- read.csv(paste(dir,"hcws2010.csv",sep=""))

unvax <- data[data$vacpH1N1==0,]

# "Among 599 HCWs who reported that they had not received pH1N1 vaccine, 74 (12%) had pH1N1 antibody titer >=1:40 by VN."
round(c(sum(unvax$VNantibody==1),sum(unvax$VNantibody==1)/nrow(unvax)),2)

# "In a random sample of 59/599 tested by HAI, 9 (15%) had antibody titer >=1:40."
temp <- unvax[!is.na(unvax$HIantibody)&!(unvax$confirmedH1N1==2&!is.na(unvax$confirmedH1N1)),]
round(c(sum(temp$HIantibody==1),sum(temp$HIantibody==1)/nrow(temp)),2)

# "Among the 599 HCWs, 19 (3.2%) reported laboratory-confirmed pH1N1 infection during the first wave,
round(nrow(unvax[!is.na(unvax$confirmedH1N1)&unvax$confirmedH1N1==2,])/c(1,nrow(unvax)),3)
# and 58% (95% CI: 34-80) of those 19 had antibody titer >=1:40 by VN,
temp <- nrow(unvax[!is.na(unvax$confirmedH1N1)&unvax$confirmedH1N1==2&unvax$VNantibody==1,]);
round(c(temp/19,binom.test(temp,19)$conf[1:2]),2)
# and 75% (95% CI: 49-91) had antibody titer >=1:40 by HAI."
temp <- nrow(unvax[!is.na(unvax$confirmedH1N1)&unvax$confirmedH1N1==2&unvax$HIantibody==1,]);
round(c(temp/19,binom.test(temp,19)$conf[1:2]),2)

# "Among the 574 HCWs who did not report laboratory-confirmed pH1N1 infection, 11% (95% CI: 8.5-14) had antibody titer >=1:40 by VN."
temp <- unvax[!is.na(unvax$confirmedH1N1)&unvax$confirmedH1N1<=1,];
round(c(nrow(temp[temp$VNantibody==1,])/nrow(temp),binom.test(nrow(temp[temp$VNantibody==1,]),nrow(temp))$conf[1:2]),2)

# "Of the 599 HCWs, 338 (57%) reported ILI since July 2009 and 19% (95% CI: 15-23) of those HCWs had antibody titer >=1:40 by VN
temp <- unvax[!is.na(unvax$ILI)&unvax$ILI==1,];
round(c(nrow(temp),nrow(temp[temp$VNantibody==1,])/nrow(temp),binom.test(nrow(temp[temp$VNantibody==1,]),nrow(temp))$conf[1:2]),2)
# versus 4.3% (95% CI: 2.2-7.6) of the 255 HCWs who did nit report ILI during the pandemic."
temp <- unvax[!is.na(unvax$ILI)&unvax$ILI==0,];
round(c(nrow(temp),nrow(temp[temp$VNantibody==1,])/nrow(temp),binom.test(nrow(temp[temp$VNantibody==1,]),nrow(temp))$conf[1:2]),3)

# End of script.
