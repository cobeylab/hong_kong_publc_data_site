#
# R syntax to reproduce information for Figure S3 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 13, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe.r")

my.col <-c(gray(0.2),gray(0.4),gray(0.7),gray(0.9))
title<-c("Fever","Shivering/Chills","Feeling Tired","Headache","Cough","Muscle Pain","Swelling", "Redness", "Bruising", "Pain/soreness")

################################################################################
# Plot Daily Adverse Reaction                                                  #
################################################################################

mat.ar1 <-function (p,q,r,s) matrix(c(
dim(arr2[arr2[,p]==2 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],
dim(arr2[arr2[,p]==1 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],

dim(arr2[arr2[,q]==2 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],
dim(arr2[arr2[,q]==1 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],

dim(arr2[arr2[,r]==2 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],
dim(arr2[arr2[,r]==1 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],

dim(arr2[arr2[,s]==2 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1],
dim(arr2[arr2[,s]==1 & arr2$intervention=="TIV",])[1]/dim(arr2[arr2$intervention=="TIV",])[1]

), ncol=4, byrow=F,)

mat.ar0 <-function (p,q,r,s) matrix(c(
dim(arr2[arr2[,p]==2 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],
dim(arr2[arr2[,p]==1 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],


dim(arr2[arr2[,q]==2 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],
dim(arr2[arr2[,q]==1 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],

dim(arr2[arr2[,r]==2 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],
dim(arr2[arr2[,r]==1 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],


dim(arr2[arr2[,s]==2 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1],
dim(arr2[arr2[,s]==1 & arr2$intervention=="placebo",])[1]/dim(arr2[arr2$intervention=="placebo",])[1]

), ncol=4, byrow=F,)

windows(height=13, width=20)
layout(matrix(1:12, ncol=4, byrow=TRUE))
par(mar=c(2,3,2,0), oma=c(0.5,1.5,0,0))

for ( i in c(8,9,7,6,0,5,2,1,3,4)){
barplot(mat.ar1(2+(i+1)*4,3+(i+1)*4,4+(i+1)*4,5+(i+1)*4),ylim=c(0,0.5),col=my.col[c(1,2)],beside=F,space=c(0,1.5,1.5,1.5,1.5),axes=F)
barplot(mat.ar0(2+(i+1)*4,3+(i+1)*4,4+(i+1)*4,5+(i+1)*4),ylim=c(0,0.5),col=my.col[c(3,4)],beside=F,space=c(1,1.5,1.5,1.5,1.5),axes=F,add=T)
axis(2,at=0:5*0.1,labels=paste(0:5*10,"%",sep=""),las=1)
axis(1,at=c(1,3.5,6,8.5),labels=c("Day 0", "Day 1", "Day 2", "Day 3"))
mtext(title[i+1],side=3, line=-1,cex=0.9,font=2)
mtext("Proportion",side=2,line=3, cex=0.8)
}
#plot legend manually
plot(-100,-100,xlab=F,ylab=F,axes=F,xlim=c(0,6),ylim=c(0,0.25))
#adjust x coordinates
x <- 0
y <- -0.05
points(x+0.2,y+0.24,pch=22,bg=my.col[2],cex=3)
text(x+0.3,y+0.24,"mild",pos=4)
points(x+0.2,y+0.215,pch=22,bg=my.col[1],cex=3)
text(x+0.3,y+0.215,"moderate",pos=4)
text(x+3.05,y+0.235,"TIV")

text(x+2.5,y+0.23,as.symbol("}"),cex=2)

points(x+0.2,y+0.15,pch=22,bg=my.col[4],cex=3)
text(x+0.3,y+0.15,"mild",pos=4)
points(x+0.2,y+0.125,pch=22,bg=my.col[3],cex=3)
text(x+0.3,y+0.125,"moderate",pos=4)
text(x+3.45,y+0.145,"placebo")

text(x+2.5,y+0.14,as.symbol("}"),cex=2)

##-------------------------------------------------------------
m <-function (x) matrix( c(
dim(arr2[arr2[,x]==0 & arr2$intervention=="TIV",])[1],
dim(arr2[(arr2[,x]==1 | arr2[,x]==2 )& arr2$intervention=="TIV",])[1],

dim(arr2[arr2[,x]==0 & arr2$intervention=="placebo",])[1],
dim(arr2[(arr2[,x]==1 | arr2[,x]==2 )& arr2$intervention=="placebo",])[1] )
,ncol=2,byrow=F)

sig <-matrix(rep(NA,40),ncol=40)
for (i in 1:40){
sig[1,i] <-round(fisher.test(m(5+i))$p.value,2)  }

#
# End of script.
#


