#
# R syntax to reproduce information for Figure 2 from:
#
# Ng S, Cowling BJ, Fang VJ, Chan KH, Ip DKM, et al.
# Effects of oseltamivir treatment on duration of clinical illness
# and viral shedding and household transmission of influenza virus.
# CID, 2010; 50:.
#
# Last updated by Vicky Fang and Sophie Ng
# February 5, 2010

require(chron)
source("http://www.hku.hk/bcowling/influenza/Oseltamivir_scripts/dataframe.r")  
source("http://www.hku.hk/bcowling/influenza/Serial_scripts/functions.r")

turnbull.plot <- function(turnbull.object, line.type=1, line.col=1, line.width=1, max.time=100, CI.col=gray(0.8))
{
    #
    # function to add horizontal lines
    # vertical lines will be added only for equivalence classes with length 0
    #
    obj <- turnbull.object
    m <- length(obj$cum.pi)
        lines(x=c(0, obj$eq.timeL[1]), y=rep(1, 2), lty=line.type, col=line.col, lwd=line.width)
        for(i in 1:(m-1)){
            lines(c(obj$eq.timeR[i], obj$eq.timeL[i+1]), rep(1-obj$cum.pi[i], 2), lty=line.type,
                col=line.col, lwd=line.width)
        }
        lines(x=c(obj$eq.timeL[m], max.time), y=rep(0, 2), lty=line.type, col=line.col, lwd=line.width)
        if(obj$eq.timeL[1]==obj$eq.timeR[1]) lines(rep(obj$eq.timeL[1], 2),
            c(1, 1-obj$cum.pi[1]), lty=line.type, col=line.col, lwd=line.width)
        for(i in 2:m){
            if(obj$eq.timeL[i]==obj$eq.timeR[i]) lines(rep(obj$eq.timeL[i], 2),
                c(1-obj$cum.pi[i-1], 1-obj$cum.pi[i]), lty=line.type,
                    col=line.col, lwd=line.width)
        }
    invisible(0)
}

# oseltamivir start day
tamiflu <- merge(av[av$member==0&av$av=="tamiflu",],hchar[c("hhID","clinic_date","clinic_day")],by="hhID",all.x=TRUE)
tamiflu$onset.tami <- dates(as.character(tamiflu$dayfrom),format="d/m/y")-dates(as.character(tamiflu$clinic_date),format="d/m/y")+tamiflu$clinic_day
tamiflu$onset.tami <- tamiflu$onset.tami+1

data <- merge(vshed,clinic7[c("hhID","onsettime")],all.x=T)
data <- merge(data,tamiflu[c("hhID","onset.tami")],all.x=T)
data$onset.tami[is.na(data$onset.tami)] <- 0
data$trunct <- floor(data$onsettime/2)
data$timeL[data$timeL==0] <- 0.01

data.g1 <- data[data$onset.tami!=0,]  # tamiflu
data.g2 <- data[data$onset.tami==0,]
data.g1$exact <- data.g2$exact <- 3
data.g1$exact[data.g1$timeL==data.g1$timeR] <- 1
data.g2$exact[data.g2$timeL==data.g2$timeR] <- 1

g2.turnbull <- make.turnbull(data.g2)

# Analysis by subgroups
sg11 <- data.g1[data.g1$onset.tami==1,]
sg12 <- data.g1[data.g1$onset.tami==2,]

sg11.turnbull <- make.turnbull(sg11)
sg12.turnbull <- make.turnbull(sg12)

windows(width=6,height=4.5)

plot(0, xlim=c(0,12), ylim=c(0,1), type="n",
  xlab="Days since symptom onset", ylab="Proportion of cases with viral shedding", bty="n", axes=FALSE)

turnbull.plot(g2.turnbull, line.type=1, line.col="black", line.width=2, max.time=50)

turnbull.plot(sg12.turnbull, line.type=2, line.col="blue", line.width=2, max.time=12)
turnbull.plot(sg11.turnbull, line.type=4, line.col="red", line.width=2, max.time=12)

axis(1, pos=0, lwd=1, cex.axis=1.0, at=3*0:5)
axis(2, pos=0, lwd=1, cex.axis=1.0, las=1)
legend(6,1.0,legend=c("no antiviral","oseltamivir 24-48h after onset","oseltamivir within 24h of onset"),
       lty=c(1,2,4),col=c("black","blue","red"),cex=0.75)

#
# End of script
#

