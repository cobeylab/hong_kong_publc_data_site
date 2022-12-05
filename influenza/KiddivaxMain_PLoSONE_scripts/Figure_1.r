#
# R syntax to reproduce information for Figure 1 from:
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

my.col <- c(rgb(0.4,0.4,0.7, alpha=0.8), rgb(0.5,0.5,0.5, alpha=0.8))

# Function for plot

fun <-function (agegp1,agegp2, ycol, v07, v08, plab1,plab2) {          ## 2 variables for each age groups
  plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
  points(jitter(rep(1, length(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07])), factor=10),log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[1])
  points(jitter(rep(2.5,length(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07])), factor=7), log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[2])
  points(jitter(rep(4.5,length(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07])), factor=3),log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[1])
  points(jitter(rep(6,length(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07])), factor=2),log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[2] )
  
  segments(1.5,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           1.5,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(1.5,median(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  segments(3.2,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           3.2,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(3.2,median(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  segments(5.2,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           5.2,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(5.2,median(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  segments(6.5,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           6.5,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(6.5,median(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  {text(5.2, log(10240)+0.8, plab1, cex=0.9)}
  {text(6.5, log(10240)+0.8, plab2, cex=0.9)}
  
  axis(1,at=c(0,1,2.5,4.5,6,7), labels=c("", "Y", "O", "Y", "O", ""))
  axis(2, at=c(log(5*2^(0:12))), labels=c(5*2^(0:12)), las=1)
  segments(3.5,0, 3.5,log(20480), lty=3)
}

fun2 <-function (agegp1,agegp2, ycol, v07, v08, plab1,plab2) {          ## 2 variables for each age groups
  plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
  points(jitter(rep(1, length(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07])), factor=10),log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[1])
  points(jitter(rep(2.5,length(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07])), factor=7), log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[2])
  points(jitter(rep(4.5,length(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07])), factor=3),log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[1])
  points(jitter(rep(6,length(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07])), factor=2),log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]), col=my.col[2] )
  
  segments(1.5,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           1.5,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(1.5,median(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  segments(3.2,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           3.2,quantile(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(3.2,median(log(dat[,ycol][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  segments(5.2,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           5.2,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(5.2,median(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp1 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  segments(6.5,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.25),
           6.5,quantile(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T,0.75))
  points(6.5,median(log(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp2 & dat$vac0809==v08 & dat$vac0708==v07]),na.rm=T),pch="-")
  
  {text(5.2-0.2, log(10240)+0.8, plab1, cex=0.9)}
  {text(6.5, log(10240)+0.8, plab2, cex=0.9)}
  
  axis(1,at=c(0,1,2.5,4.5,6,7), labels=c("", "Y", "O", "Y", "O", ""))
  axis(2, at=c(log(5*2^(0:12))), labels=rep("",13), las=1)
  segments(3.5,0, 3.5,log(20480), lty=3)
}

# Function for p-value
getp <- function (agegp, ycol, v07, v08) {
   p.value <- wilcox.test(dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp & dat$vac0809==0 & dat$vac0708==0],
                          dat[,ycol+1][dat$intervention=="TIV" & dat$age%in%agegp & dat$vac0809==v08 & dat$vac0708==v07])$p.value              
   if(p.value>=0.01) p <- paste("p=",round(p.value,2),sep="")            
   else p <- "p<0.01"  
   p
}

# Create plot

windows( width=3.5*5, height=3*6)
layout(matrix(1:30, nrow=6,byrow=F), widths=c(rep(1.8,5)), heights=c(rep(1.7,6)))

par(mar=c(2,1.5,0,0), oma=c(0,0,0,0), mgp=c(1,0.5,0))

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(5.5,(log(5)+log(20480))/2,"Seasonal\nA(H1N1)", font=2)

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(5.5,(log(5)+log(20480))/2,"Seasonal\nA(H3N2)", font=2)

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(5.5,(log(5)+log(20480))/2,"B\n/Florida", font=2)

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(5.5,(log(5)+log(20480))/2,"B\n/Brisbane", font=2)

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(5.5,(log(5)+log(20480))/2,"Pandemic\nA(H1N1)", font=2)

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(3.5,log(10),"Not received TIV\nin 2007-08 or 2008-09", font=2)
mtext("PRE     POST", side=1, line=1.2, cex=0.8)

fun(6:8,9:12, 29, 0,0,"Ref", "Ref")
fun(6:8,9:12, 33, 0,0,"Ref", "Ref")
fun(6:8,9:12, 41, 0,0,"Ref", "Ref")
fun(6:8,9:12, 37, 0,0,"Ref", "Ref")
fun(6:8,9:12, 25, 0,0,"Ref", "Ref")

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(3.5,log(10),"Received TIV\nin 2007-08 only", font=2)
mtext("PRE     POST", side=1, line=1.2, cex=0.8)

fun2(6:8,9:12, 29, 1,0, getp(6:8,29,1,0),getp(9:17,29,1,0))
fun2(6:8,9:12, 33, 1,0, getp(6:8,33,1,0),getp(9:17,33,1,0))
fun2(6:8,9:12, 41, 1,0, getp(6:8,41,1,0),getp(9:17,41,1,0))
fun2(6:8,9:12, 37, 1,0, getp(6:8,37,1,0),getp(9:17,37,1,0))
fun2(6:8,9:12, 25, 1,0, getp(6:8,25,1,0),getp(9:17,25,1,0))

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(3.5,log(10),"Received TIV\nin 2008-09 only", font=2)
mtext("PRE     POST", side=1, line=1.2, cex=0.8)

fun2(6:8,9:12, 29, 0,1, getp(6:8,29,0,1),getp(9:17,29,0,1))
fun2(6:8,9:12, 33, 0,1, getp(6:8,33,0,1),getp(9:17,33,0,1))
fun2(6:8,9:12, 41, 0,1, getp(6:8,41,0,1),getp(9:17,41,0,1))
fun2(6:8,9:12, 37, 0,1, getp(6:8,37,0,1),getp(9:17,37,0,1))
fun2(6:8,9:12, 25, 0,1, getp(6:8,25,0,1),getp(9:17,25,0,1))

plot(0,0,type="n", ylim=c(log(5),log(20480)+0.5), xlim=c(0, 7), main="", xlab="", ylab="", axes=F)
text(3.5,log(10),"Received TIV\nin 2007-08 & 2008-09", font=2)
mtext("PRE     POST", side=1, line=1.2, cex=0.8)

fun2(6:8,9:12, 29, 1,1, getp(6:8,29,1,1),getp(9:17,29,1,1))
fun2(6:8,9:12, 33, 1,1, getp(6:8,33,1,1),getp(9:17,33,1,1))
fun2(6:8,9:12, 41, 1,1, getp(6:8,41,1,1),getp(9:17,41,1,1))
fun2(6:8,9:12, 37, 1,1, getp(6:8,37,1,1),getp(9:17,37,1,1))
fun2(6:8,9:12, 25, 1,1, getp(6:8,25,1,1),getp(9:17,25,1,1))

#
# End of script.
#

