#
# R syntax including functions needed in the analyses described in:
#
# Cheng CKY, Cowling BJ, Chan KH, Fang VJ, Seto WH, et al.
# Factors affecting QuickVue Influenza A+B rapid test performance in
# the community setting.
# Diagnostic Microbiology and Infectious Disease, 2009; 65(1): 35-41.
#
# Last updated by Calvin Cheng and Vicky Fang
# Sep 08, 2009

# generate random number
rand <-function(x){
set.seed(100)
sort(sample(1:200, x))
}
#end

#
#function for display results when put in different parameters 
#

snspresult <- function(cdc, entry, parameter) {
  tab <- cdc[c(entry==parameter & !is.na(entry)),]
  tab <- table(tab[c("QVpos","cultpos")])
  maketab2(tab)
}
#end

#function for rearrange the 2 by 2 table to ++ +- -+ -- and display the result
maketab2 <-function(oritab){
  tab <-matrix(NA, ncol=2, nrow=2)
  tab[1,1]<-oritab[2,2]
  tab[1,2]<-oritab[2,1]
  tab[2,1]<-oritab[1,2]
  tab[2,2]<-oritab[1,1]
  print(tab)
  print(sum(tab))
  
  snsp(tab)
}

#
# function for calculating sens, spec and exact 95% CIs "snsp"
#

snsp <- function(tab){
cdc <-tab
total <- sum(tab)
# 95% CI for sensitivity
sn <- binom.test(x=cdc[1,1], n=cdc[1,1] + cdc[2,1])

# 95% CI for specificity
sp <- binom.test(x=cdc[2,2], n=cdc[1,2] + cdc[2,2])

#95% CI for PPV
ppv <- binom.test(x=cdc[1,1], n=cdc[1,1] + cdc[1,2])
#95% CI for NPV
npv <- binom.test(x=cdc[2,2], n=cdc[2,1] + cdc[2,2])


# Reformat and present results
output <- data.frame(
  Value = c(sn[[1]] / sn[[2]], sp[[1]] / sp[[2]] ,  
  ppv[[1]] / ppv[[2]], npv[[1]] / npv[[2]]),
 lowerCI = c(sn[[4]][1], sp[[4]][1], ppv[[4]][1], npv[[4]][1]),
  upperCI = c(sn[[4]][2], sp[[4]][2], ppv[[4]][2], npv[[4]][2]),
  Total = c(total, total, total,total),
  row.names=c("Sensitivity", "Specificity","PPV","NPV"))

round(output, 2)
}

# END OF SCRIPTS


