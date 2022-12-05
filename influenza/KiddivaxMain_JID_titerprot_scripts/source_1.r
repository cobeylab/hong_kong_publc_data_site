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


# rescale factor obtained from the proxy for influenza A and B

dir <- "http://sph.hku.hk/data/KiddivaxMainV2/"

kdata <- read.csv(paste(dir, "serology.csv", sep=""))
swab <- read.csv(paste(dir, "swab.csv", sep=""))
random <- read.csv(paste(dir, "randomcode.csv", sep=""))
fluproxy <- read.csv("http://www.hku.hk/bcowling/influenza/influenza_proxy_1998to2013.csv")

##
time0 <-"31/07/2009"
fluproxy <- fluproxy[614:679,]
fluproxy$ti <- as.numeric(as.Date(fluproxy$Week.ending,format="%d/%m/%Y")-as.Date(time0,format="%d/%m/%Y"))-1

proxy <- data.frame(ti=fluproxy$ti,proxy=fluproxy$B.proxy)
s.proxy <-data.frame(ti=1:max(proxy$ti), s.proxy=predict(smooth.spline(proxy$ti, proxy$proxy, spar=0.5),1:max(proxy$ti))$y)

proxy.p <- data.frame(ti=fluproxy$ti,proxy.p=fluproxy$A.H1N1pdm.proxy+0.005)
s.proxy.p <-data.frame(ti=1:max(proxy.p$ti), s.proxy.p=predict(smooth.spline(proxy.p$ti, proxy.p$proxy.p, spar=0.6),1:max(proxy.p$ti))$y)

dat <- kdata[kdata$member==0,]
dat <- merge(dat,random, by="hhID", all.x=T)
dat <- dat[dat$intervention=="placebo",]

r <-((sum(s.proxy.p$s.proxy.p[90:500])/sum(s.proxy$s.proxy[90:500]) )    /  (prop.table(table(dat$post.season.pH1/dat$postvax.pH1 >=4)) /prop.table(table(dat$post.season.B.Brisbane/dat$postvax.B.Brisbane >=4))))[2]


#
# End of script.
#

