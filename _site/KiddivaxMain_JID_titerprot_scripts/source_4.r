#
# R syntax to generate information used in:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

#
# To create a dataframe on predicted titer for plotting estimated protection in each subjects
#
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_2.r")

# Read in data
kdata <- read.csv(paste(dir, "serology.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))

kdata2 <- kdata[kdata$member==0 & kdata$prevax!="" & kdata$postvax!="" & !is.na(kdata$postvax.B.Brisbane),]
kdata2 <-merge(kdata2,random, by="hhID",all.x=T)
kdata2 <-merge(kdata2, demog[,c("hhID", "member", "age")], all.x=T)

kdata2$postday <-as.numeric (as.Date( as.character(kdata2$postvax), format="%d/%m/%Y") -   as.Date( as.character(kdata2$prevax), format="%d/%m/%Y"))

kdata3 <- kdata2

w1<-model.b1.age$coef[c(2,4)] #coef[1]: caltime ; coef[2]: caltime*(age<9)
w0<-model.b0$coef[2]

mrow <- 361

for ( i in 1) {
                temp <-data.frame(hhID=rep(kdata3$hhID[i],each=mrow),TIV=1*(kdata3$intervention[i]=="TIV"), age=rep(kdata3$age[i], each=mrow), postvax.b=NA, t=(1:mrow)-1, postday=rep(kdata3$postday[i],each=mrow))
                temp2 <-log2(runif(1,kdata3$postvax.B.Brisbane[i]-5, (kdata3$postvax.B.Brisbane[i] * 2)-5 ))
                days <- temp$postday[1]:(mrow-1)
                bdays <- 0:temp$postday[1]
                temp$postvax.b[days+1] <- temp2 +  (days-days[1]) * (w0 * (temp$TIV[1] %in% 0) +   w1[1] * (temp$TIV[1]%in% 1) +  w1[2] * (temp$TIV[1]%in% 1) * (temp$age[1]<9)  )
                temp$postvax.b[max(bdays+1):1] <- temp2 -  bdays * (w0 * (temp$TIV[1] %in% 0) +    w1[1] * (temp$TIV[1]%in% 1) +  w1[2] * (temp$TIV[1]%in% 1) * (temp$age[1]<9))
                pred <-temp
                }

for ( i in 2:nrow(kdata3)) {
                temp <-data.frame(hhID=rep(kdata3$hhID[i],each=mrow),TIV=1*(kdata3$intervention[i]=="TIV"), age=rep(kdata3$age[i], each=mrow), postvax.b=NA, t=(1:mrow)-1, postday=rep(kdata3$postday[i],each=mrow))
                temp2 <-log2(runif(1,kdata3$postvax.B.Brisbane[i]-5, (kdata3$postvax.B.Brisbane[i] * 2)-5 ))
                days <- temp$postday[1]:(mrow-1)
                bdays <- 0:temp$postday[1]
                temp$postvax.b[days+1] <- temp2 +  (days-days[1]) * (w0 * (temp$TIV[1] %in% 0) +   w1[1] * (temp$TIV[1]%in% 1) +  w1[2] * (temp$TIV[1]%in% 1) * (temp$age[1]<9)  )
                temp$postvax.b[max(bdays+1):1] <- temp2 -  bdays * (w0 * (temp$TIV[1] %in% 0) +   w1[1] * (temp$TIV[1]%in% 1) +  w1[2] * (temp$TIV[1]%in% 1) * (temp$age[1]<9))
                pred<-rbind(pred, temp)
                            }

Pred.B.titer <-pred


#
# End of script.
#


