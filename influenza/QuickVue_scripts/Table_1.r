#
# R syntax to reproduce information in Table 1 (QuickVue performance vs gold standard) from:
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
source ("http://www.hku.hk/bcowling/influenza/QuickVue_scripts/qv_functions.r")

############################
#2 by 2 table for sn and sp#
############################

#for adult and children
a1<- snspresult (cdc, cdc$agegroup,1) #
tmp1 <- cbind( a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$agegroup,2) #
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<- snspresult (cdc, cdc$agegroup,3) #
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$agegroup,4) #
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
a5<-snspresult (cdc, cdc$agegroup,5) #
tmp5 <- cbind(a5[1,1:3],a5[2,1:4])
a6<-snspresult (cdc, cdc$agegroup,6) #
tmp6 <- cbind(a6[1,1:3],a6[2,1:4])

out <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6) #combining
out <-out[c(7,1:6)] #rearranging order
names(out) <- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("age 0-5","6-10","11-15","16-30","31-50","50+")  #renaming
table1p1 <- out

# for sex  
a1<- snspresult (cdc, cdc$male,1) #male
tmp1 <- cbind(a1[1,1:3],a1[2,1:4])
a2 <- snspresult (cdc, cdc$male,0) #female
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
out <- rbind(tmp1,tmp2)
out <-out[c(7,1:6)]
#renaming
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("male","female")
table1p2 <- out

#for symptom onset 12,24,36,48 hrs 
a1<- snspresult (cdc, cdc$onsettime,1) #symptom onet <12 hours
tmp1 <- cbind(a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$onsettime,2) #symptom onet <24 hours
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<-snspresult (cdc, cdc$onsettime,3) #symptom onet <36 hours
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$onsettime,4) #symptom onet <48 hours
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
a5<- snspresult (cdc, cdc$onsettime,5) 
tmp5 <- cbind(a5[1,1:3],a5[2,1:4])
out <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5)

out <-out[c(7,1:6)]
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("onset <12hrs","12-24hrs","24-36hrs","36-48hrs",">48hrs")
table1p3 <- out

#071204 for experienced in qv done
a1<- snspresult (cdc, cdc$qvdonegp,1)
tmp1 <- cbind(a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$qvdonegp,2)
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<- snspresult (cdc, cdc$qvdonegp,3)
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$qvdonegp,4)
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
out <- rbind(tmp1,tmp2,tmp3,tmp4)
out <-out[c(7,1:6)]
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("QV <30","30-60","60-90",">90")
table1p4 <- out

#for different symptoms
a1<- snspresult (cdc, cdc$rnose,1) #
tmp1 <- cbind( a1[1,1:3],a1[2,1:4])
a2<- snspresult (cdc, cdc$cough,1) #
tmp2 <- cbind(a2[1,1:3],a2[2,1:4])
a3<- snspresult (cdc, cdc$sthroat,1) #
tmp3 <- cbind(a3[1,1:3],a3[2,1:4])
a4<- snspresult (cdc, cdc$headache,1) #
tmp4 <- cbind(a4[1,1:3],a4[2,1:4])
a5<-snspresult (cdc, cdc$pmuscle,1) #
tmp5 <- cbind(a5[1,1:3],a5[2,1:4])
a6<-snspresult (cdc, cdc$measure_fever,1) #
tmp6 <- cbind(a6[1,1:3],a6[2,1:4])
a7<-snspresult (cdc, cdc$fcs,1) #
tmp7 <- cbind(a7[1,1:3],a7[2,1:4])
out <- rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7) #combining
out <-out[c(7,1:6)] #rearranging order
names(out)<- c("Total","Sn","lowerCI","upperCI","Sp","lowerCI","upperCI")
row.names(out)<- c("runny nose","cough","sore throat","headache","muscle pains","fever","fever+c/s")
table1p5 <- out

rbind(table1p1,table1p2,table1p3,table1p4,table1p5)



#add p value (fisher's test)

# sensitivity

# age
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$agegroup[cdc$cultpos==1]))$p.value,2)
# sex
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$male[cdc$cultpos==1]))$p.value,2)
# onsettime
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$onsettime[cdc$cultpos==1]))$p.value,2)
# qvdonegp
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$qvdonegp[cdc$cultpos==1]))$p.value,2)
# symptoms
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$rnose[cdc$cultpos==1]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$cough[cdc$cultpos==1]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$sthroat[cdc$cultpos==1]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$headache[cdc$cultpos==1]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$pmuscle[cdc$cultpos==1]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$measure_fever[cdc$cultpos==1]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==1], cdc$fcs[cdc$cultpos==1]))$p.value,2)


# specificity

round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$agegroup[cdc$cultpos==0]))$p.value,2)
# sex
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$male[cdc$cultpos==0]))$p.value,2)
# onsettime
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$onsettime[cdc$cultpos==0]))$p.value,2)
# qvdonegp
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$qvdonegp[cdc$cultpos==0]))$p.value,2)
# symptoms
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$rnose[cdc$cultpos==0]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$cough[cdc$cultpos==0]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$sthroat[cdc$cultpos==0]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$headache[cdc$cultpos==0]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$pmuscle[cdc$cultpos==0]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$measure_fever[cdc$cultpos==0]))$p.value,2)
round(fisher.test(table(cdc$QVpos[cdc$cultpos==0], cdc$fcs[cdc$cultpos==0]))$p.value,2)


# END OF SCRIPTS

