#
# R syntax to reproduce Table S2 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

require(chron)
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_scripts/dataframe.r")

#infect <- merge(demog[c("hhID","member","age")],random,all.x=T)
symp$date <- as.Date(symp$date,format="%d/%m/%Y")
symp$symp <- 1*(rowSums(symp[,c(5:10,12)])>=1)

# pH1N1
ph1 <- swab[swab$Swine.H1=="P",1:3]
ph1 <- ph1[order(ph1$hhID,ph1$member,ph1$date),]
ph1u <- unique(ph1[1:2])
ph1u$ILI <- ph1u$ARI <- ph1u$asymp <- ph1u$musclepain <- ph1u$runnynose <- ph1u$cough <- ph1u$sorethroat <-
ph1u$headache <- ph1u$chills <- ph1u$fever <- ph1u$fluend <- ph1u$flustart <- NA
for(i in 1:nrow(ph1u)){
   tmp <- as.numeric(as.character(ph1$date[ph1$hhID==ph1u$hhID[i]&ph1$member==ph1u$member[i]]))
   ph1u$flustart[i] <- min(tmp)
   ph1u$fluend[i] <- max(tmp)
   symp.tmp <- symp[symp$hhID==ph1u$hhID[i]&symp$member==ph1u$member[i]&symp$date>=as.Date(as.character(ph1u$flustart[i]),format="%Y%m%d")-5&
                           symp$date<=as.Date(as.character(ph1u$fluend[i]),format="%Y%m%d")+5,c(1,2,12,5:10,15,13,14)]
   ph1u[i,5:14] <- 1*(colSums(symp.tmp[,3:12])>=1)
   ph1u[i,12] <- 1-ph1u[i,12]
}
nrow(ph1u)

# sH3N2
sh3 <- swab[swab$H3=="P",1:3]
sh3 <- sh3[order(sh3$hhID,sh3$member,sh3$date),]
sh3u <- unique(sh3[1:2])
sh3u$ILI <- sh3u$ARI <- sh3u$asymp <- sh3u$musclepain <- sh3u$runnynose <- sh3u$cough <- sh3u$sorethroat <-
sh3u$headache <- sh3u$chills <- sh3u$fever <- sh3u$fluend <- sh3u$flustart <- NA
for(i in 1:nrow(sh3u)){
   tmp <- as.numeric(as.character(sh3$date[sh3$hhID==sh3u$hhID[i]&sh3$member==sh3u$member[i]]))
   sh3u$flustart[i] <- min(tmp)
   sh3u$fluend[i] <- max(tmp)
   symp.tmp <- symp[symp$hhID==sh3u$hhID[i]&symp$member==sh3u$member[i]&symp$date>=as.Date(as.character(sh3u$flustart[i]),format="%Y%m%d")-5&
                           symp$date<=as.Date(as.character(sh3u$fluend[i]),format="%Y%m%d")+5,c(1,2,12,5:10,15,13,14)]
   sh3u[i,5:14] <- 1*(colSums(symp.tmp[,3:12])>=1)
   sh3u[i,12] <- 1-sh3u[i,12]
}
nrow(sh3u)

# B
flub <- swab[swab$FluB=="P",1:3]
flub <- flub[order(flub$hhID,flub$member,flub$date),]
flubu <- unique(flub[1:2])
flubu$ILI <- flubu$ARI <- flubu$asymp <- flubu$musclepain <- flubu$runnynose <- flubu$cough <- flubu$sorethroat <-
flubu$headache <- flubu$chills <- flubu$fever <- flubu$fluend <- flubu$flustart <- NA
for(i in 1:nrow(flubu)){
   tmp <- as.numeric(as.character(flub$date[flub$hhID==flubu$hhID[i]&flub$member==flubu$member[i]]))
   flubu$flustart[i] <- min(tmp)
   flubu$fluend[i] <- max(tmp)
   symp.tmp <- symp[symp$hhID==flubu$hhID[i]&symp$member==flubu$member[i]&symp$date>=as.Date(as.character(flubu$flustart[i]),format="%Y%m%d")-5&
                           symp$date<=as.Date(as.character(flubu$fluend[i]),format="%Y%m%d")+5,c(1,2,12,5:10,15,13,14)]
   flubu[i,5:14] <- 1*(colSums(symp.tmp[,3:12])>=1)
   flubu[i,12] <- 1-flubu[i,12]
}
nrow(flubu)

# calculate p-values
tab <- as.data.frame(cbind(colSums(ph1u[5:14]),colSums(sh3u[5:14]),colSums(flubu[5:14]))); names(tab) <- c("pH1","sH3","B")
tab$pH1.prop <- round(tab$pH1/nrow(ph1u),2);tab$sH3.prop <- round(tab$sH3/nrow(sh3u),2);tab$B.prop <- round(tab$B/nrow(flubu),2)
tab$p.value <- NA
for(i in 1:nrow(tab)){
  tab$p.value[i] <- round(fisher.test(rbind(tab[i,1:3],c(nrow(ph1u),nrow(sh3u),nrow(flubu))-tab[i,1:3]))$p.value,2)
}
tab <- tab[c(order(tab$pH1[1:7],decreasing=T),8:10),c(1,4,2,5,3,6,7)]
tab

# asymp changed to 2,0,0 according to telephone records

#
# End of script.
#

