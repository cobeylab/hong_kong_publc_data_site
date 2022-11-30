#
# R syntax to reproduce information for eTable 2 from:
#
# Klick B, Nishiura H, Ng S, et al.
# Transmissibility of seasonal and pandemic influenza in a cohort
# of households in Hong Kong in 2009
# Epidemiology. 2011 Nov;22(6):793-6.
#
# Last updated by Klick B, Fang VJ and Cowling BJ.
# January 3, 2012
#

dir <- "http://sph.hku.hk/data/KiddivaxPilotV1/"

sero <- read.csv(paste(dir, "serology_m.csv", sep=""))
demog <- read.csv(paste(dir, "demographic_m.csv", sep=""))
code <- read.csv(paste(dir, "randomcode_h.csv", sep=""))

tab <- matrix(rep(NA,9*12),ncol=12,
              dimnames=list(c("All:pH1","sH1","sH3","<=1:160:pH1","sH1","sH3","<=1:40:pH1","sH1","sH3"),
                            c("CPI-child","CI_low","CI_up","CPI-adult","CI_low","CI_up","SAP-child","CI_low","CI_up","SAP-adult","CI_low","CI_up")))

sero$ph1.post.titer[(sero$hhID==9186|sero$hhID==9194|sero$hhID==9213|sero$hhID==9218)&sero$member==0] <- 20

#convert ph1 titer to numeric & NA for Positive PositiveOne & 5 for Negative
sero$ph1.mid.titer[sero$ph1.mid.titer=="Negative"] <-5
sero$ph1.post.titer[sero$ph1.post.titer=="Negative" | sero$ph1.post.titer=="PositiveOne"] <-5

sero$ph1s <- 1*(as.numeric(as.character(sero$ph1.post.titer))/as.numeric(as.character(sero$ph1.mid.titer))>=4&as.numeric(as.character(sero$ph1.post.titer))>=40)
sero$ph1s[is.na(sero$ph1s)&sero$ph1.post.titer=="Positive"] <- 1
sero$ph1s[is.na(sero$ph1s)&sero$ph1.post.titer==5] <- 0
sero$sh1s <- 1*(sero$sh1.posts/sero$sh1.mids>=4&sero$sh1.posts>=40)
sero$sh3s <- 1*(sero$sh3.posts/sero$sh3.mids>=4&sero$sh3.posts>=40)

c_data<-merge(x=sero,y=demog,by=c("hhID","member"),all.x=T,all.y=T)

c_data$sh1s[c_data$hhID==9101 & c_data$member==0]<-0
c_data$sh1s[c_data$hhID==9110 & c_data$member==2]<-0
c_data$sh3s[c_data$hhID==9122 & c_data$member==4]<-0
c_data$sh1s[c_data$hhID==9128 & c_data$member==3]<-0
c_data$ph1s[c_data$hhID==9128 & c_data$member==3]<-0
c_data$sh1s[c_data$hhID==9161 & c_data$member==1]<-0
c_data$sh1s[c_data$hhID==9167 & c_data$member==3]<-0
c_data$sh3s[c_data$hhID==9167 & c_data$member==3]<-0
c_data$sh1s[c_data$hhID==9193 & c_data$member==0]<-0
c_data$sh3s[c_data$hhID==9193 & c_data$member==0]<-0
c_data$sh1s[c_data$hhID==9195 & c_data$member==3]<-0
c_data$sh1s[c_data$hhID==9203 & c_data$member==3]<-0
c_data$sh3s[c_data$hhID==9206 & c_data$member==3]<-0
c_data$ph1s[c_data$hhID==9208 & c_data$member==1]<-0
c_data$sh1s[c_data$hhID==9211 & c_data$member==0]<-0
c_data$sh3s[c_data$hhID==9214 & c_data$member==3]<-0

c_data<-subset(c_data, hhID!=9108)
c_data<-subset(c_data, hhID!=9172)

a <- c_data[!is.na(c_data$sh1.mids)&!is.na(c_data$sh3.mids),]
a$adult <- 1*(a$age>14)

tab <- matrix(rep(NA,10*8),ncol=8,
              dimnames=list(c("All:child","adult","sH1<=1:160:child","adult","sH3<=1:160:child","adult","sH1<=1:40:child","adult","sH3<=1:40:child","adult"),
                            c("N","%","pH1","%","sH1","%","sH3","%")))
                            
tab[1:2,1] <- table(a$adult); tab[1:2,3] <- table(a$adult[a$ph1s==1&!is.na(a$ph1s)])
tab[1:2,5] <- table(a$adult[a$sh1s==1&!is.na(a$sh1s)]); tab[1:2,7] <- table(a$adult[a$sh3s==1&!is.na(a$sh3s)])

tab[3:4,1] <- table(a$adult[a$sh1.mids<=160&!is.na(a$sh1.mids)])
tab[3:4,3] <- table(a$adult[a$sh1.mids<=160&!is.na(a$sh1.mids)&a$ph1s==1&!is.na(a$ph1s)])
tab[3:4,5] <- table(a$adult[a$sh1.mids<=160&!is.na(a$sh1.mids)&a$sh1s==1&!is.na(a$sh1s)])
tab[3:4,7] <- table(a$adult[a$sh1.mids<=160&!is.na(a$sh1.mids)&a$sh3s==1&!is.na(a$sh3s)])

tab[5:6,1] <- table(a$adult[a$sh3.mids<=160&!is.na(a$sh3.mids)])
tab[5:6,3] <- table(a$adult[a$sh3.mids<=160&!is.na(a$sh3.mids)&a$ph1s==1&!is.na(a$ph1s)])
tab[5:6,5] <- table(a$adult[a$sh3.mids<=160&!is.na(a$sh3.mids)&a$sh1s==1&!is.na(a$sh1s)])
tab[5:6,7] <- table(a$adult[a$sh3.mids<=160&!is.na(a$sh3.mids)&a$sh3s==1&!is.na(a$sh3s)])

tab[7:8,1] <- table(a$adult[a$sh1.mids<=40&!is.na(a$sh1.mids)])
tab[7:8,3] <- table(a$adult[a$sh1.mids<=40&!is.na(a$sh1.mids)&a$ph1s==1&!is.na(a$ph1s)])
tab[7:8,5] <- table(a$adult[a$sh1.mids<=40&!is.na(a$sh1.mids)&a$sh1s==1&!is.na(a$sh1s)])
tab[7:8,7] <- table(a$adult[a$sh1.mids<=40&!is.na(a$sh1.mids)&a$sh3s==1&!is.na(a$sh3s)])

tab[9:10,1] <- table(a$adult[a$sh3.mids<=40&!is.na(a$sh3.mids)])
tab[9:10,3] <- table(a$adult[a$sh3.mids<=40&!is.na(a$sh3.mids)&a$ph1s==1&!is.na(a$ph1s)])
tab[9:10,5] <- table(a$adult[a$sh3.mids<=40&!is.na(a$sh3.mids)&a$sh1s==1&!is.na(a$sh1s)])
tab[9:10,7] <- table(a$adult[a$sh3.mids<=40&!is.na(a$sh3.mids)&a$sh3s==1&!is.na(a$sh3s)])

tab[,4] <- round(tab[,3]/tab[,1],2)
tab[,6] <- round(tab[,5]/tab[,1],2)
tab[,8] <- round(tab[,7]/tab[,1],2)

tab[c(3,5,7,9),2] <- round(tab[c(3,5,7,9),1]/tab[1,1],2)
tab[c(4,6,8,10),2] <- round(tab[c(4,6,8,10),1]/tab[2,1],2)

tab

# End of script.


