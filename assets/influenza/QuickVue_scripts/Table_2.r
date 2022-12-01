#
# R syntax to reproduce information in Table 2 (QuickVue performance vs gold standard) from:
#
# Cheng CKY, Cowling BJ, Chan KH, Fang VJ, Seto WH, et al.
# Factors affecting QuickVue Influenza A+B rapid test performance in
# the community setting.
# Diagnostic Microbiology and Infectious Disease, 2009; 65(1): 35-41.
#
# Last updated by Calvin Cheng and Vicky Fang
# Sep 08, 2009


dir <- "http://sph.hku.hk/data/HongKongNPIpilotV3/"

cdc <- read.csv(paste(dir,"QuickVue_data.csv",sep=""))
source ("http://www.hku.hk/bcowling/influenza/QuickVue_scripts/add_groups.r")

# Sn for flu A and B having qPCR result

cdc$logqPCR <- log(cdc$qPCR,10)
cdc$logqPCR[cdc$logqPCR == -Inf] <- 0

cdc$qPCRgp <-cut(cdc$logqPCR, breaks=c(-1 ,0, 4, 5 , 6, 10))
cdc$qPCRgp <- as.numeric(cdc$qPCRgp)
pcrava <- cdc[!is.na(cdc$qPCR),]

#function for outputing the result
sn <- function(gp, gpQVpos){
# handle only 1 case in table
   new1 <- vector(length =2)
   new1[1] <- length(gpQVpos) - sum(gpQVpos) # total -ve 
   new1[2] <- sum(gpQVpos)    #total +ve
   sn <- binom.test(x= new1[2], dim(gp)[1]) 

   output <- data.frame(Total = dim(gp)[1], Value = c(sn[[1]] / sn[[2]]),
                        lowerCI = c(sn[[4]][1]), upperCI = c(sn[[4]][2]))
   round(output, 2)
}

# for influenza A
pcrcultposA <- pcrava[pcrava$cultposA ==1,]
pcrcultposA1 <- pcrcultposA[pcrcultposA$qPCRgp==1,]
pcrcultposA2 <- pcrcultposA[pcrcultposA$qPCRgp==2,]
pcrcultposA3 <- pcrcultposA[pcrcultposA$qPCRgp==3,]
pcrcultposA4 <- pcrcultposA[pcrcultposA$qPCRgp==4,]
pcrcultposA5 <- pcrcultposA[pcrcultposA$qPCRgp==5,]

 a0<- sn(pcrcultposA, pcrcultposA$QVpos) # overall influenza A
 a1<- sn(pcrcultposA1, pcrcultposA1$QVpos)
 a2<- sn(pcrcultposA2, pcrcultposA2$QVpos)
 a3<- sn(pcrcultposA3, pcrcultposA3$QVpos)
 a4<- sn(pcrcultposA4, pcrcultposA4$QVpos)
 a5<- sn(pcrcultposA5, pcrcultposA5$QVpos)
 out <- rbind(a0,a1,a2,a3,a4,a5)
 row.names(out)<- c("Overall","Undetectable","<10^4 copies/ml",
 "10^4 to 10^5 copies/ml","10^5 to 10^6 copies/ml","> 10^6 copies/ml")
 table2p1 <- out

#for influenza B
pcrcultposB <- pcrava[pcrava$cultposB ==1,] # overall influenza B
pcrcultposB1 <- pcrcultposB[pcrcultposB$qPCRgp==1,]
pcrcultposB2 <- pcrcultposB[pcrcultposB$qPCRgp==2,]
pcrcultposB3 <- pcrcultposB[pcrcultposB$qPCRgp==3,]
pcrcultposB4 <- pcrcultposB[pcrcultposB$qPCRgp==4,]
pcrcultposB5 <- pcrcultposB[pcrcultposB$qPCRgp==5,]

a0<- sn(pcrcultposB, pcrcultposB$QVpos)
a1<- sn(pcrcultposB1, pcrcultposB1$QVpos)
a2<- sn(pcrcultposB2, pcrcultposB2$QVpos)
a3<- sn(pcrcultposB3, pcrcultposB3$QVpos)
a4<- sn(pcrcultposB4, pcrcultposB4$QVpos)
a5<- sn(pcrcultposB5, pcrcultposB5$QVpos)
out <- rbind(a0,a1,a2,a3,a4,a5)
table2p2 <- out

cbind(table2p1,table2p2)


# END OF SCRIPTS

