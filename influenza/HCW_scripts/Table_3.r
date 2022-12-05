#
# R syntax to reproduce information for Table 3 from:
#
# Zhou Y, Ng DMW, Seto WH, et al.
# Seroprevalence of antibody to pandemic influenza A (H1N1) 2009 among
# healthcare works after the first wave in Hong Kong.
# Journal of Hospital Infection, 2011.
#
# Last updated by Zhou Y and Fang VJ.
# September 20, 2011

dir <- "http://sph.hku.hk/data/"
data <- read.csv(paste(dir,"hcws2010.csv",sep=""))
tab3 <- matrix(NA,ncol=6,nrow=5)

hcws <- data[data$vacpH1N1==0,c("Record_ID","age","VNantibody")]
hcws$agegp <- as.character(factor(cut(hcws$age,c(0,24,34,44,54,100)),labels=c(1,2,3,4,5)))

tab3[,1] <- table(hcws$agegp,hcws$VNantibody)[,2]
tab3[,2] <- table(hcws$agegp)
tab3[,3] <- round(tab3[,1]/tab3[,2],2)
for(i in 1:5){
  tab3[i,4:5] <- round(binom.test(tab3[i,1],tab3[i,2])$conf[1:2],2)
}

tab3[1,6] <- round(chisq.test(matrix(c(20,8,114-20,49-8),ncol=2))$p.value,2)
tab3[2,6] <- round(chisq.test(matrix(c(15,25,130-15,125-25),ncol=2))$p.value,2)
tab3[3,6] <- round(chisq.test(matrix(c(13,21,122-13,162-21),ncol=2))$p.value,2)
tab3[4,6] <- round(chisq.test(matrix(c(4,14,81-4,190-14),ncol=2))$p.value,2)
tab3[5,6] <- round(fisher.test(matrix(c(1,6,19-1,72-6),ncol=2))$p.value,2)

colnames(tab3) <- c("n","N","%","CI_low","CI_up","P-value")
rownames(tab3) <- c("Age:19-24","25-34","35-44","45-54","55-64")

tab3

# End of script.

