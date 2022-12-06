#
# R syntax to reproduce information for Figure S1 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# November 8, 2010

dir <- "http://sph.hku.hk/data/KiddivaxPilotV1/"
qm <- read.csv(paste(dir, "QMHisolate.csv", sep=""))

wk <-c(0,rep(NA,25))
  for (i in 1:12){
      wk[i+1] <-wk[i]+ nrow(qm[qm$year==2008 & qm$month==i,])}
  for (i in 1:12){
      wk[13+i]<-wk[12+i]+nrow(qm[qm$year==2009 & qm$month==i,])}

mth <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec",
         "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec" )

m <- matrix(rep(NA,dim(qm)[1]*3), nrow=3)
m[1,] <-qm[,5]/qm[,10]
m[2,] <-qm[,6]/qm[,10]
m[3,] <-qm[,9]/qm[,10]
m <- m[,1:104]

windows(width=7, height=3.9)
par(mar=c(0,0.5,0.5,0), oma=c(0,1,0,0))
my.col <-c(gray(0.3),gray(0.6),gray(0.8),rgb(0,0,1,alpha=0.1),rgb(1,0,0.5,alpha=0.1))
plot(-100.-100,ylim=c(-0.05,0.4),ylab="",xlab="",xlim=c(0,104),axes=F)

polygon(c(wk[11]-1,wk[11]-1,wk[12]-1,wk[12]-1),c(0.45,0,0,0.45),col=my.col[4],border=NA)
polygon(c(wk[12]-0.5,wk[12]-0.5,wk[13],wk[13]),c(0.45,0,0,0.45),col=my.col[4],border=NA)
polygon(c(wk[16],wk[16],wk[17],wk[17]),c(0.45,0,0,0.45),col=my.col[4],border=NA)
polygon(c(wk[21],wk[21],wk[23],wk[23]),c(0.45,0,0,0.45),col=my.col[4],border=NA)

barplot(m,ylim=c(-0.1,0.45),col=my.col[c(1,2,3)],beside=F,space=0,axes=F,add=T)
axis(1,line=-2.6,at=wk,labels=c(rep("",26)))
axis(1,line=-2.6,at=wk[c(1,13,25)],labels=c(rep("",3)),tcl=-1)
axis(2, pos=-0.1,at=c(0,0.45), labels=c(rep("",2)),col.ticks="white")
axis(2, pos=-0.1,at=c(1:4*0.1), labels=c("10%","20%","30%","40%"),las=1,cex.axis=0.8, hadj=0.5, tcl=-0.2)
mtext("Proportion of positive isolates",side=2,line=0.5,adj=0.7,cex=0.8)

for (i in 1:24){
text((wk[i]+wk[i+1])/2,-0.010,mth[i],cex=0.7,adj=0.5)
}

text(104/4,-0.03,2008,cex=0.8)
text(104/4*3,-0.03,2009,cex=0.8)

#adjust x coordinates
x<-5
y<-0.1
#plot legend manually
points(x+0.2,y+0.26,pch=22,bg=my.col[3],cex=2)
text(x+0.3,y+0.26,"Pandemic A(H1N1)",pos=4,cex=0.8)
points(x+0.2,y+0.245,pch=22,bg=my.col[2],cex=2)
text(x+0.3,y+0.245,"Seasonal A",pos=4,cex=0.8)
points(x+0.2,y+0.23,pch=22,bg=my.col[1],cex=2)
text(x+0.3,y+0.23,"Seasonal B",pos=4,cex=0.8)

text((wk[12]+wk[11])/2-1,0.4165,"blood\ndraw\n+\nvaccine/\nplacebo\n",cex=0.5,pos=1)
text((wk[12]+wk[13])/2-0.5,0.425,"blood\ndraw",cex=0.5,pos=1)
text((wk[16]+wk[17])/2,0.425,"blood\ndraw",cex=0.5,pos=1)
text((wk[21]+wk[23])/2-0.5,0.425,"blood\ndraw",cex=0.5,pos=1)

#
# End of script.
#

