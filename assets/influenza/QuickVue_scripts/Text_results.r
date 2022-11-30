#
# R syntax to reproduce information of results in main text from:
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

# Main text - Result - Paragraph 2

# overall by QV vs culture

oritab <-table(cdc[c("QVpos","cultpos")])
maketab2(oritab)

# stratified by flu type

fluA <- cdc[cdc$cultposB!=1,]
maketab2(table(fluA[c("QVposA","cultposA")]))

fluB <- cdc[cdc$cultposA!=1,]
maketab2(table(fluB[c("QVposB","cultposB")]))


# Main text - Result - Paragraph 3

# t test for viral load between QV +ve/-ve, separated by influenza A and influenza B

dat <-cdc[!is.na(cdc$qPCR),] 

dat$logqPCR <- log(dat$qPCR,10)
dat$logqPCR[dat$logqPCR == -Inf] <- 0

#seperate into 3 groups
group1 <- dat[dat$QVres == 3 & dat$culture==0,] # both -ve
group2 <- dat[dat$QVres == 3 & dat$culture!=0,] # QV-ve , culture +ve
group3 <- dat[dat$QVres != 3 & dat$culture!=0,] # QV +ve , culture +ve

group2a <- group2[group2$culture =="A",] 
group2b <- group2[group2$culture =="B",]
group3a <- group3[group3$culture =="A",]
group3b <- group3[group3$culture =="B",]

groupa <-rbind (group2a,group3a)
groupb <-rbind (group2b,group3b)

# Influenza A
groupa$QVres[groupa$QVres==3] <-0
round(t.test(groupa$logqPCR~groupa$QVres)$p.value,3)
dim(groupa)[1]

# Influenza B
groupb$QVres[groupb$QVres==3] <-0
groupb$QVres[groupb$QVres==2] <-1
round(t.test(groupb$logqPCR~groupb$QVres)$p.value,3)
dim(groupb)[1]


# Main text - Result - Paragraph 3 - Last sentence

cdc2 <- cdc
cdc2$goldpos[cdc2$cultpos ==1 | (!is.na(cdc2$qPCR)&cdc2$qPCR !=0) ] <-1 
cdc2$goldpos[is.na(cdc2$goldpos)] <-0

goldtab <-table(cdc2[c("QVpos","goldpos")])
maketab2(goldtab)  #compare QV with new gold standard 


# END OF SCRIPTS

