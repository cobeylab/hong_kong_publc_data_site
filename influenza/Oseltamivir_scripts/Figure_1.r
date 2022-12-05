#
# R syntax to reproduce information for Figure 1 from:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

source("http://www.hku.hk/bcowling/influenza/Oseltamivir_scripts/dataframe.r")
library(survival)

windows(width=8,height=3)
layout(matrix(1:3,nrow=1),width=c(1,1),heights=c(1,1))
plot(survfit(Surv(all_time, all_event)~av, data=output.all),main="(a) All Symptoms",xlim=c(0,10),lty=c(2:1),lwd=2,
     xlab="Time from symptom onset \n to alleviation, days", ylab="Proportion of patients with any symptoms", axes=F )
axis(2, at=seq(0,1,0.2),las=1)
axis(1, at=seq(0,11,2))
legend(5,0.4,c("No antiviral","Oseltamivir"),lty=c(2:1),lwd=2,cex=0.85)
text(2.5,0.005,"log-rank p=0.01")

plot(survfit(Surv(f_time, f_event)~av, data=output.fever),main="(b) Fever",xlim=c(0,10),lty=c(2:1),lwd=2,
     xlab="Time from symptom onset \n to alleviation, days", ylab="Proportion of patients with fever", axes=F )
axis(2, at=seq(0,1,0.2),las=1)
axis(1, at=seq(0,11,2))
legend(5,0.4,c("No antiviral","Oseltamivir"),lty=c(2:1),lwd=2,cex=0.85)
text(2.5,0.005,"log-rank p=0.13")

plot(survfit(Surv(r_time, r_event)~av, data=output.rs),main="(c) Respiratory Symptoms",xlim=c(0,10),lty=c(2:1),lwd=2,
     xlab="Time from symptom onset \n to alleviation, days", ylab="Proportion of patients with respiratory symptoms", axes=F )
axis(2, at=seq(0,1,0.2),las=1)
axis(1, at=seq(0,11,2))
legend(5,0.4,c("No antiviral","Oseltamivir"),lty=c(2:1),lwd=2,cex=0.85)
text(2.5,0.005,"log-rank p=0.03")

#
# End of script
#

