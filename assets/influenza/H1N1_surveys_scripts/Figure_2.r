#
# Script to reproduce information in Figure 2 from:
#
# Cowling BJ, Ng DMW, Ip DKM, Liao Q, Lam WWT, Wu JT, et al.
# Community Psychological and Behavioral Responses through the First Wave
# of the 2009 Influenza A(H1N1) Pandemic in Hong Kong
# JID 2010; 202(6):867-876.
#
# Last updated by Diane Ng, Vicky Fang and Ben Cowling.
# October 26, 2010


raw.cases <- read.csv("http://sph.hku.hk/data/cum_casesv1.csv")

require(IDPmisc)
require(fields)

survey <- data.frame(dates=c(1, 14, 42, 55, 69, 79, 89, 100, 115, 134, 155, 176, 197),  
# The 7th date should be 84 (jitter a bit for plot)
   
# panel 2
  anxiety=c(1.83, 1.82, 1.79, 1.78, 1.78, 1.79, 1.79, 1.81, 1.80, 1.81, 1.85, 1.84, 1.85),

# panel 3
  worry=c(56.3, 50.1, 49.8, 47.3, 42.0, 40.5, 45.5, 43.1, 45.7, 44.3, 42.7, 40.0, 40.2),
  absolute=c(21.3, 10.6, 10.7, 14.9, 10.5, 12.5, 13.1, 9.7, 12.0, 11.6, 11.7, 10.1, 10.4),
  relative=c(4.6, 5.4, 6.2, 7, 4.6, 5.2, 5.1, 3.6, 6.2, 5.5, 5.9, 3.8, 3.2),
  SARS=c(57.2, 24.4, 16.4, 11.6, 11.1, 15.3, 14.1, 14.2, 13.8, 13.2, 15.3, 13.2, 13.7),  

# panel 4
  droplets=c(93.1, 94.5, 92.8, 93.5, 92.6, NA, 92.7, 93.6, 93.1, 93.6, 93.5, 93, 92),
  dhc=c(62.1, 63.6, 64.3, 63.0, 59.7, NA, 61.1, 62, 62.5, 61.3, 60.3, 64.5, 59.4),
  ihc=c(69.6, 71.3, 69.4, 69.4, 68.9, NA, 67.9, 65.3, 66.8, 66.7, 64.3, 66.5, 67.5),
  of=c(72.3, 72, 70.5, 71.9, 69.9, NA, 71.6, 69.7, 79, 77, 76.3, 75.2, 78),
  cold=c(NA, 23.8, 24.6, 19.2, 24.1, NA, 23.8, 22.8, 24.1, 27.6, 29.1, 32.9, 36.8),

# panel 5
  handwash=c(65.1, 66.6, 62.7, 62.3, 61.7, 70.4, 63.3, 64.5, 60.8, 58.8, 58.8, 56.6, 57.0),
  handwashhome=c(NA, NA, 84.9, 85.5, 86.6, NA, 86.2, 84.8, 84, 83.8, 82.9, 85.6, 84.1),
  touch=c(11, 16.3, 13.5, 12.8, 10.8, 16.1, 11.0, 12.3, 11.1, 10.6, 10.9, 12.2, 10.6),
  soup=c(74.3, NA, 76.1, 76.5, 75.6, NA, 72.1, 72, 72.9, 73.5, 70.3, 74.3, 72.4),
  handwashtouch=c(NA, NA, 61.7, 43.8, 44, NA, 44.6, 44.5, 44.7, 41.7, 41.1, 41.6, 39.9),
  clean=c(NA, NA, 11.8, 17.4, 15, NA, 13, 13.7, 8, 9, 8.9, 7, 6.9),

# panel 6
  covermouth=c(75, 72.4, 73.6, 69.8, 72.1, NA, 70.8, 71.9, 66.3, 70.4, 68.9, 72.2, 71.3),
  mask=c(7.8, NA, 9.0, 10.8, 9.3, 8.7, 7.9, 8.7, 8.4, 9.4, 9.8, 8.4, 6.7),  
  eatout=c(NA, NA, 6.5, 10.6, 5.7, 6.5, 5.4, 5, 3.5, 4.6, 4.1, 2.7, 3.2),
  publictransport=c(NA, NA, 2.2, 5.3, 4.2, 2.6, 2.2, 3, 2.1, 1.5, 1.6, 1.3, 1.6),
  crowded=c(NA, NA, 20.4, 29.9, 20.1, 13.1, 17.3, 19.1, 15.2, 14.8, 11.4, 9.7, 10.3),
  reschedule=c(NA, NA, 17.3, 17.9, 12.9, NA, 12.2, 12.3, 7.4, 4.9, 6, 3.4, 2.5))

cases <- data.frame(day=1:206)
cases$total <- rep(0, 206)

for (i in 1:177){
cases$total[raw.cases$day[i]] <- raw.cases$Total[i]
}

# You may need to save the graph in a local place.
png("d:/fig2.png", width=15, height=24, units="in", res=200, bg="transparent")
layout(matrix(1:6, ncol=1), heights=c(1.6,0.8,1.35,1.35,1.35,1.35))

#1 epi curve
par(mar=c(1,7,1,0.5))
barplot(cases$total,
xlim=c(-4, 250), ylim=c(0,800), 
col=gray(0.5), border=gray(0.5), width=0.6, space=2/3,
axes=FALSE, axisnames=FALSE,
ylab="", cex.lab=1)
axis(1, pos=0, at=c(-6,1,32,62,93,124,154,185,215), labels=rep("", 9), cex.axis=1)
axis(2, pos=-6, las=1, cex.axis=1.2)
mtext("Number of cases", side=2, line=0.5, cex=1.2)
title(main="(a)", cex.main=1.2)

az <- 0.3
cex.t <- 1

#school closure 12 June = 43
t1 <- 43
y1 <- 180
Arrows(t1,y1,t1,40, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t1-3,t1),c(y1,y1),lwd=1)
text (t1-18.5, y = y1+10, "Closure of kindergartens and ", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)
text (t1-20.5, y = y1-10, "primary schools begins (12 Jun)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)

#designated fever clinics open 13 Jun = 45
t2 <- 45
y2 <- 240
Arrows(t2,y2,t2,40, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t2-3,t2),c(y2,y2),lwd=1)
text (t2-11.5, y = y2+10, "8 designated", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)
text (t2-17.5, y = y2-10, "flu clinics open (13 Jun)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)

#adjusted hospitalization criteria 29 Jun = 60
t3 <- 60
y3 <- 300
Arrows(t3,y3,t3,83, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t3-3,t3),c(y3,y3),lwd=1)
text (t3-24.5, y = y3, "Adjusted hospitalization criteria (29 Jun)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)

#first death 16 Jul = 77
t5 <- 77
y5 <- 360
Arrows(t5,y5,t5,118, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t5-3,t5),c(y5,y5),lwd=1)
text (t5-14.5, y = y5+10, "First death (16 Jul)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)
text (t5-19.5, y = y5-10, "42-year-old Philippine sailor", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)
     
#revised guidelines for antiviral treatment 20 Jul = 81
t6 <- 81
y6 <- 420
Arrows(t6,y6,t6,150, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t6-3,t6),c(y6,y6),lwd=1)
text (t6-28.5, y = y6, "Revised guidelines for antiviral treatment (20 Jul)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)
    
#first local death 31 Jul = 92
t7 <- 92
y7 <- 480
Arrows(t7,y7,t7,229, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t7-3,t7),c(y7,y7),lwd=1)
text (t7-16.5, y = y7+10, "First local death (31 Jul)",adj = NULL, pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)
text (t7-28.5, y = y7-10, "58-year-old local male with chronic heart disease",adj = NULL, pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)

#flu clinics resumed service 17 Aug = 109
t8 <- 109
y8 <- 540
Arrows(t8,y8,t8,328, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t8-3,t8),c(y8,y8),lwd=1)
text (t8-41, y = y8, "Designated flu clinics resume service for chronic disease patients (17 Aug)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)

#school reopen 1 Sep = 124
t8 <- 124
y8 <- 600
Arrows(t8,y8,t8,332, size = az, width = 2.8, open = F,sh.adj=1,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo="Black", verbose=T)
lines(c(t8-3,t8),c(y8,y8),lwd=1)
text (t8-17, y = y8, "Schools reopen (1 Sep)", adj = NULL,
     pos = NULL, offset = 0.4, vfont = NULL, cex = cex.t, col = NULL, font = NULL)

# containment phase: to 29 Jun (0-90)
tc1 <- 0
tc2 <- 60
ytc <- 725
lines(c(tc1,tc2),c(ytc,ytc),col=gray(0.7))
lines(c(tc1,tc1),c(ytc-5,ytc+5),col=gray(0.7))
lines(c(tc2,tc2),c(ytc-5,ytc+5),col=gray(0.7))
text (30, y =ytc+30, "Containment phase (1 May to 29 Jun)", adj = NULL,
pos = NULL, offset = 0.5, vfont = NULL,cex = cex.t, col =gray(0.4), font = NULL)

# mitigation phase: 10 Jun to (41- ) 
tm1 <- 41
tm2 <- 206
ytm <- 700
Arrows(tm1,ytm,tm2,ytm,size=az,width=2.8,open=T,sh.adj=0,sh.lwd=1,sh.lty=1,h.lwd=1,h.col.bo=gray(0.7),sh.col=gray(0.7),verbose=T)
lines(c(tm1,tm1),c(ytm-5,ytm+5),col=gray(0.7))
text (130, y =ytm+25, "Mitigation phase (starting from 10 Jun)", adj = NULL,
pos = NULL, offset = 0.5, vfont = NULL,cex = cex.t, col = gray(0.4), font = NULL)

#2 anxiety level
par(mar=c(1,7,1,0.5))
plot(NA, type="n", axes=FALSE, xlim=c(-4, 250), ylim=c(1, 4),
ylab="", cex.lab=1)

lines(survey$dates, survey$anxiety, type="b", pch=4, cex=1.1)
axis(1, pos=1, at=c(-6,1,32,62,93,124,154,185,215), labels=rep("", 9), cex.axis=1)
axis(2, pos=-6, at=1:4, las=1, cex.axis=1.2)
mtext("STAI score", side=2, line=0.5, cex=1.2)
title(main="(b)", cex.main=1.2)

#3 worry, absolute & relative susceptibility, severity compared to SARS
par(mar=c(1,7,1,0.5))
plot(NA, type="n", axes=FALSE, xlim=c(-4, 250), ylim=c(0, 60), xlab="", ylab="", cex.lab=1)
lines(survey$dates, survey$worry, lwd=1.5, type="b", pch=4, cex=1.1)
lines(survey$dates, survey$absolute, lty=2, lwd=1.5, type="b", pch=17, cex=1.4)
lines(survey$dates, survey$relative, lty=3, lwd=1.5, type="b", pch=18, cex=1.8)
lines(survey$dates, survey$SARS, lty=4, lwd=1.5, type="b", pch=15, cex=1.1)

axis(1, pos=0, at=c(-6,1,32,62,93,124,154,185,215), labels=rep("", 9), cex.axis=1)
axis(2, pos=-6, at=0:6*10,
labels=paste(0:6*10, "%", sep=""), las=1, cex.axis=1.2)
mtext("Proportion", side=2, line=0.5, cex=1.2)
text(200,survey$worry[13],"Worry if developed ILI",cex=1.5,adj = c(0,0))
text(200,survey$SARS[13],"Perceived severity compared to SARS",cex=1.5,adj = c(0,0))
text(200,survey$absolute[13],"Absolute susceptibility",cex=1.5,adj = c(0,0))
text(200,survey$relative[13],"Relative susceptibility",cex=1.5,adj = c(0,0))
title(main="(c)", cex.main=1.2)

#4 knowledge
par(mar=c(1,7,1,0.5))
plot(NA, type="n", axes=FALSE, xlim=c(-4, 250), ylim=c(0, 100), xlab="", ylab="", cex.lab=1)
lines(survey$dates[-c(6)], survey$droplets[-c(6)], lwd=1.5, type="b", pch=4, cex=1.1)
lines(survey$dates[-c(6)], survey$dhc[-c(6)],lty=2, lwd=1.5, type="b", pch=17, cex=1.4)
lines(survey$dates[-c(6)], survey$ihc[-c(6)],lty=3, lwd=1.5, type="b", pch=18, cex=1.8)
lines(survey$dates[-c(6)], survey$of[-c(6)],lty=4, lwd=1.5, type="b", pch=15, cex=1.1)
lines(survey$dates[-c(6)], survey$cold[-c(6)],lty=5, lwd=1.5, type="b", pch=2, cex=1.1)
axis(1, pos=0, at=c(-6,1,32,62,93,124,154,185,215), labels=rep("", 9), cex.axis=1)
axis(2, pos=-6, at=0:5/5*100,
labels=paste(0:5/5*100, "%", sep=""), las=1, cex.axis=1.2)
mtext("Proportion", side=2, line=0.5, cex=1.2)
text(200,survey$droplets[13],"Via droplets",cex=1.5,adj = c(0,0))
text(200,survey$of[13],"Via oral-feacal",cex=1.5,adj = c(0,0))
text(200,survey$ihc[13],"Via indirect contact",cex=1.5,adj = c(0,0))
text(200,survey$dhc[13],"Via direct contact",cex=1.5,adj = c(0,0))
text(200,survey$cold[13],"Via cold weather",cex=1.5,adj = c(0,0)) 
title(main="(d)", cex.main=1.2)

#5 personal hygiene
par(mar=c(1,7,1,0.5))
plot(NA, type="n", axes=FALSE, xlim=c(-4, 250), ylim=c(0, 100), xlab="", ylab="", cex.lab=1)
lines(survey$dates[-c(2)], survey$handwash[-c(2)], lwd=1.5, type="b", pch=4, cex=1.1)
lines(survey$dates[-c(6)], survey$handwashhome[-c(6)], lty=2, lwd=1.5, type="b", pch=17, cex=1.4)
lines(survey$dates[-c(2,6)], survey$soup[-c(2,6)], lty=4, lwd=1.5, type="b", pch=18, cex=1.8)
lines(survey$dates[-c(1,2,6)], survey$handwashtouch[-c(1,2,6)], lty=5, lwd=1.5, type="b", pch=15, cex=1.1)
lines(survey$dates[-c(6)], survey$clean[-c(6)], lty=6, lwd=1.5, type="b", pch=2, cex=1.1)
text(200,survey$handwashhome[13],"Washed hands after returning home",cex=1.5,adj = c(0,0))
text(200,survey$soup[13],"Used liquid soup when washing hands",cex=1.5,adj = c(0,0))
text(200,survey$handwash[13],"Washed hands after sneezing, 
coughing or touching nose",cex=1.5,adj = c(0,0.5))
text(200,survey$handwashtouch[13],"Washed hands after 
touching common objects",cex=1.5,adj = c(0,0.5))
text(200,survey$clean[13],"Cleaned or disinfected home",cex=1.5,adj = c(0,0)) 
axis(1, pos=0, at=c(-6,1,32,62,93,124,154,185,215), labels=rep("", 9), cex.axis=1)
axis(2, pos=-6, at=0:5/5*100,
labels=paste(0:5/5*100, "%", sep=""), las=1, cex.axis=1.2)
mtext("Proportion", side=2, line=0.5, cex=1.2) 
title(main="(e)", cex.main=1.2) 
 
#6 social distancing
par(mar=c(1,7,1,0.5))
plot(NA, type="n", axes=FALSE, xlim=c(-4, 250), ylim=c(0, 80), xlab="", ylab="", cex.lab=1)
lines(survey$dates[-c(6)], survey$covermouth[-c(6)], lwd=1.5, type="b", pch=4, cex=1.1)
lines(survey$dates[-c(2)], survey$mask[-c(2)], lty=2, lwd=1.5, type="b", pch=17, cex=1.4)
lines(survey$dates[-c(1,2)], survey$eatout[-c(1,2)], lty=3, lwd=1.5, type="b", pch=18, cex=1.8)
lines(survey$dates[-c(1,2)], survey$publictransport[-c(1,2)], lty=4, lwd=1.5, type="b", pch=15, cex=1.1)
lines(survey$dates[-c(1,2)], survey$crowded[-c(1,2)], lty=5, lwd=1.5, type="b", pch=5, cex=1.1)
lines(survey$dates[-c(1,2,6)], survey$reschedule[-c(1,2,6)], lty=6, lwd=1.5, type="b", pch=24, cex=1.1)

legend(180,55,legend=c("Covered mouth when sneezing or coughing","Avoided crowded places","Rescheduled travel plans","Wore mask","Avoided eating out","Avoided public transport"),lty=c(1,5,6,2,3,4),lwd=1.2,bty="n",cex=1.5) 

axis(1, pos=0, at=c(-6,1,32,62,93,124,154,185,215), labels=c("","May 1","Jun 1",
"Jul 1","Aug 1","Sep 1","Oct 1","Nov 1","Dec 1"), cex.axis=1.2)
axis(2, pos=-6, at=0:5/5*80,
labels=paste(0:5/5*80, "%", sep=""), las=1, cex.axis=1.2)
mtext("Proportion", side=2, line=0.5, cex=1.2)

title(main="(f)", cex.main=1.2)

dev.off()
  
# End of script
