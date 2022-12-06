#
# R syntax to reproduce information for Table 5 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 21, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe_cross.r")

# tabulate ili by infected strains in winter
w.s1<-table(index.pre$iliw[index.pre$sh1w==1], index.pre$intervention[index.pre$sh1w==1])
w.s3<-table(index.pre$iliw[index.pre$sh3w==1], index.pre$intervention[index.pre$sh3w==1])
w.b<-table(index.pre$iliw[index.pre$flo.w==1], index.pre$intervention[index.pre$flo.w==1])
w.n<-table(index.pre$iliw[!(index.pre$sh1w==1 |index.pre$sh3w==1| index.pre$flo.w==1)], index.pre$intervention[!(index.pre$sh1w==1 |index.pre$sh3w==1| index.pre$flo.w==1)])

# tabulate ili by infected strains in summer
s.s1<-table(index.pre$ilis[index.pre$sh1s==1], index.pre$intervention[index.pre$sh1s==1])
s.s3<-table(index.pre$ilis[index.pre$sh3s==1], index.pre$intervention[index.pre$sh3s==1])
s.b<-table(index.pre$ilis[index.pre$flo.s==1], index.pre$intervention[index.pre$flo.s==1])
s.p<-table(index.pre$ilis[index.pre$ph1s==1], index.pre$intervention[index.pre$ph1s==1])
s.n<-table(index.pre$ilis[!(index.pre$sh1s==1 |index.pre$sh3s==1| index.pre$flo.s==1 |index.pre$ph1s==1)], index.pre$intervention[!(index.pre$sh1s==1 |index.pre$sh3s==1| index.pre$flo.s==1 |index.pre$ph1s==1)])

#Create blank table
tab5 <- matrix(rep(NA, 7*9), nrow=9,
               dimnames=list(c("winter.sH1", "winter.sH3", "winter.B", "winter.none",
                                "summer.sH1", "summer.sH3", "summer.B", "summer.pH1", "summer.none"),
                              c("TIV.n",  "TIV.N",  "TIV.%", "plac.n", "plac.N", "plac.%", "p.value")))

#fill values into blank table
temp <- c(w.s1, w.s3, w.b[1], 0,w.b[2],0, w.n, s.s1, s.s3,s.b[1],0, s.b[2],0, s.p, s.n)

for(i in 1:9){
tab5[i,1] <- temp[4+(i-1)*4]
tab5[i,2] <- temp[4+(i-1)*4]+temp[3+(i-1)*4]
tab5[i,4] <- temp[2+(i-1)*4]
tab5[i,5] <- temp[2+(i-1)*4]+temp[1+(i-1)*4]
tab5[i,7] <- round(fisher.test(matrix(c(temp[1+(i-1)*4],temp[2+(i-1)*4],temp[3+(i-1)*4],temp[4+(i-1)*4]), nrow=2, byrow=F))$p.value,2)
}

tab5[,3] <- round(tab5[,1]/tab5[,2]*100,1) #calculate % of n/N
tab5[,6] <- round(tab5[,4]/tab5[,5]*100,1) #calculate % of n/N

tab5

#
# End of script.
#
