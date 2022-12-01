#
# R syntax to reproduce information for Table 1 from:
#
# Zhou Y, Ng DMW, Seto WH, et al.
# Seroprevalence of antibody to pandemic influenza A (H1N1) 2009 among
# healthcare works after the first wave in Hong Kong.
# Journal of Hospital Infection, 2011.
#
# Last updated by Zhou Y and Fang VJ.
# September 15, 2011

dir <- "http://sph.hku.hk/data/"

data <- read.csv(paste(dir,"hcws2010.csv",sep=""))
tab1 <- matrix(NA,ncol=6,nrow=39)

unvax <- data[data$vacpH1N1==0,]    # 599 healthcare works who had not received pH1N1 vaccine

# Age group
tab1[1:6,1] <- table(cut(unvax$age,c(18.5,24,34,44,54,64)),exclude=NULL)
tab1[1:5,2] <- table(cut(unvax$age[unvax$VNantibody==1],c(18.5,24,34,44,54,64)))
tab1[1,6] <- round(fisher.test(cbind(tab1[1:5,2],tab1[1:5,1]-tab1[1:5,2]))$p.value,2) # p-value

# Sex
tab1[8:7,1] <- table(unvax$male); tab1[8:7,2] <- table(unvax$male[unvax$VNantibody==1])
tab1[7,6] <- round(chisq.test(cbind(tab1[7:8,2],tab1[7:8,1]-tab1[7:8,2]))$p.value,2)

# Occupation
tab1[9:13,1] <- table(unvax$occupation)
tab1[9:13,2] <- table(unvax$occupation[unvax$VNantibody==1])
tab1[9,6] <- round(fisher.test(cbind(tab1[9:13,2],tab1[9:13,1]-tab1[9:13,2]))$p.value,2)

# Department
tab1[14:20,1] <- table(unvax$department,exclude=NULL)
tab1[14:19,2] <- table(unvax$department[unvax$VNantibody==1])
tab1[14,6] <- round(fisher.test(cbind(tab1[14:19,2],tab1[14:19,1]-tab1[14:19,2]))$p.value,2)

# Contact with patients
tab1[21:24,1] <- table(cut(unvax$Numpatients,c(-0.1,0,5,200)),exclude=NULL)
tab1[21:23,2] <- table(cut(unvax$Numpatients[unvax$VNantibody==1],c(-0.1,0,5,200)))
tab1[21,6] <- round(chisq.test(cbind(tab1[21:23,2],tab1[21:23,1]-tab1[21:23,2]))$p.value,2)

# Acute care hospital or not
tab1[25:26,1] <- table(substr(unvax$Record_ID,1,1))
tab1[25:26,2] <- table(substr(unvax$Record_ID[unvax$VNantibody==1],1,1))
tab1[25,6] <- round(chisq.test(cbind(tab1[25:26,2],tab1[25:26,1]-tab1[25:26,2]))$p.value,2)

# No. of children
tab1[27:30,1] <- table(cut(unvax$Numchildren,c(-0.1,0,1,3)),exclude=NULL)
tab1[27:29,2] <- table(cut(unvax$Numchildren[unvax$VNantibody==1],c(-0.1,0,1,3)))
tab1[27,6] <- round(fisher.test(cbind(tab1[27:29,2],tab1[27:29,1]-tab1[27:29,2]))$p.value,2)

# Received 2009-10 vaccine
tab1[31:33,1] <- table(unvax$vac0910,exclude=NULL)
tab1[31:32,2] <- table(unvax$vac0910[unvax$VNantibody==1])
tab1[31,6] <- round(chisq.test(cbind(tab1[31:32,2],tab1[31:32,1]-tab1[31:32,2]))$p.value,2)

# Received 2008-09 vaccine
tab1[34:36,1] <- table(unvax$vac0809,exclude=NULL)
tab1[34:35,2] <- table(unvax$vac0809[unvax$VNantibody==1])
tab1[34,6] <- round(chisq.test(cbind(tab1[34:35,2],tab1[34:35,1]-tab1[34:35,2]))$p.value,2)

# Received 2007-08 vaccine
tab1[37:39,1] <- table(unvax$vac0708,exclude=NULL)
tab1[37:38,2] <- table(unvax$vac0708[unvax$VNantibody==1])
tab1[37,6] <- round(chisq.test(cbind(tab1[37:38,2],tab1[37:38,1]-tab1[37:38,2]))$p.value,2)

# proportion of VNpos (95% CI)
tab1[,3] <- signif(tab1[,2]/tab1[,1]*100,digits = 2)    # %
for(i in 1:nrow(tab1)){
   if(!is.na(tab1[i,2])) tab1[i,4:5] <- signif(binom.test(tab1[i,2],tab1[i,1])$conf[1:2]*100,digits = 2)  
}

tab1.o <- tab1[,-2]
colnames(tab1.o) <- c("Num","Antibody >=1:40 (%)","CI_low","CI_up","P-value")
rownames(tab1.o) <- c("Age:19-24","25-34","35-44","45-54","55-64","Unknown","Male","Female","Occupation:Doctor","Nurse","Clinical sp",
                      "Non-clinical sp","Other","Department:Medicine","Surgery","Emergency","Paediatrics","Other clinical","Non-clinical",
                      "Unknown","Contact with patients per day:0","1-5",">=6","Unknown","Acute care hosp","Non-acute care hosp",
                      "Num of children:0","1",">=2","Unknown","Received 09/10 vaccine:No","Yes","Unknown","Received 08/09 vaccine:No",
                      "Yes","Unknown","Received 07/08 vaccine:No","Yes","Unknown")
tab1.o

# End of script.

