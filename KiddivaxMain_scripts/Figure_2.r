#
# R syntax to reproduce Figure 2 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

require(graphics)

dir <- "http://sph.hku.hk/data/KiddivaxMainV1/"
qm <- read.csv(paste(dir, "qmdata_w_pH1.csv", sep=""))

qm$date <- as.Date(qm$date,format="%m/%d/%Y")
qm$yearmonth <- as.numeric(paste(substr(qm$date,1,4),substr(qm$date,6,7),sep=""))
qm <- qm[qm$yearmonth>=200801&qm$yearmonth<=201012,]

wk <- c(0,cumsum(table(qm$yearmonth)))

mth <- rep(c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),3)

m <- matrix(rep(NA,dim(qm)[1]*3), nrow=3)
m[1,] <-qm[,4]/qm[,5]
m[2,] <-qm[,3]/qm[,5]
m[3,] <-qm[,2]/qm[,5]

# plot

windows(width=7,height=4)
par(mar=c(0,1.5,0.5,0), oma=c(0,0.7,0,0))
my.col <- c(gray(0.3),gray(0.6),gray(0.8),rgb(0,0,1,alpha=0.1),rgb(1,0,0.5,alpha=0.1))
plot(-100.-100,ylim=c(-0.05,0.4),ylab="",xlab="",xlim=c(69,156),axes=F)

polygon(c(wk[21],wk[21],wk[26]+1,wk[26]+1),c(0.45,0,0,0.45),col=my.col[4],border=NA)
polygon(c(wk[22],wk[22],wk[27]-1,wk[27]-1),c(0.36,0,0,0.36),col=my.col[4],border=NA)
polygon(c(wk[29]-2,wk[29]-2,wk[30]-2,wk[30]-2),c(0.45,0,0,0.45),col=my.col[4],border=NA)
polygon(c(wk[33]-2,wk[33]-2,wk[37]-2,wk[37]-2),c(0.45,0,0,0.45),col=my.col[4],border=NA)

m[,1:69] <- 0
barplot(m,ylim=c(-0.1,0.45),xlim=c(69,156),col=my.col[c(1,2,3)],beside=F,space=0,axes=F,add=T)
axis(1,line=-2.7,at=wk[17:37],labels=c(rep("",21)))
axis(1,line=-2.7,at=wk[c(25,37)],labels=c(rep("",2)),tcl=-1)
axis(2, pos=69,at=c(0,0.45), labels=c(rep("",2)),col.ticks="white",tcl=-4)
axis(2, pos=69,at=c(1:4*0.1), labels=c("10%","20%","30%","40%"),las=1,cex.axis=0.7, hadj=0.5, tcl=-0.2)
mtext("Proportion of positive isolates",side=2,line=0.8,adj=0.7,cex=0.9)


for (i in 17:36){
text((wk[i]+wk[i+1])/2,-0.010,mth[i],cex=0.7,adj=0.5)
}

text(85,-0.03,2009,cex=0.8)
text(104/4*5,-0.03,2010,cex=0.8)

#-------------------------------------------------------------
#adjust x coordinates
x<-73
y<-0.1
#plot legend manually
points(x+0.2,y+0.26,pch=22,bg=my.col[3],cex=1.7)
text(x+0.3,y+0.26,"Pandemic A(H1N1)",pos=4,cex=0.6)
points(x+0.2,y+0.245,pch=22,bg=my.col[2],cex=1.7)
text(x+0.3,y+0.245,"Seasonal A(H3N2)",pos=4,cex=0.6)
points(x+0.2,y+0.23,pch=22,bg=my.col[1],cex=1.7)
text(x+0.3,y+0.23,"Seasonal B",pos=4,cex=0.6)

text((wk[22]+wk[26])/2-1,0.4165,"blood draw\n+\nvaccine/placebo\n",cex=0.6,pos=1)
text((wk[22]+wk[27])/2-1,0.36,"blood draw",cex=0.6,pos=1)
text((wk[29]+wk[30])/2-2,0.4165,"blood\ndraw",cex=0.6,pos=1)
text((wk[33]+wk[37])/2-2,0.4165,"blood\ndraw",cex=0.6,pos=1)

#
# End of script.
#


