#
# R syntax to reproduce Table 2 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective Efficacy Against Pandemic Influenza of Seasonal Influenza Vaccination 
# in Children in Hong Kong: A Randomized Controlled Trial
# CID, 2012.
#
# Last updated by Fang VJ, Sophia Ng, and Cowling BJ.
# August 20, 2012

require(Hmisc)
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_scripts/dataframe.r")

tab <-matrix(rep(NA,12*3), ncol=3,
              dimnames=list(c("sH1: pre-vax titer>=1:40", "sH1: post-vax titer>=1:40", "sH1: Mean titer rise",
                              "sH3: pre-vax titer>=1:40", "sH3: post-vax titer>=1:40", "sH3: Mean titer rise",
                              "B: pre-vax titer>=1:40", "B: post-vax titer>=1:40", "B: Mean titer rise",
                              "pH1: pre-vax titer>=1:40", "pH1: post-vax titer>=1:40", "pH1: Mean titer rise"),
                            c("TIV","Placebo","p.value")))

index.pre$agegp <- cut(index.pre$age,c(0,9,17))
immu <- index.pre[,c("hhID","member","intervention","male","vac0809","agegp","prevax.pH1","postvax.pH1","prevax.sH1","postvax.sH1",
                     "prevax.sH3","postvax.sH3","prevax.B.Brisbane","postvax.B.Brisbane")]

# ARR
arr <- data.frame(lapply(arr,function(x,...){x[is.na(x)] <- 0 ; x}))
arr <- arr[arr$returned==1,-2]
arr$pain <- arr$bruise <- arr$red <- arr$swell <-
            arr$mpain <- arr$cough <- arr$headache <- arr$tired <- arr$chill <-arr$fever <-  NA
for (i in 1:nrow(arr)){
  for (j in 46:55){
   arr[i,j] <- max(arr[i,(j-46)*4+6:9],na.rm=TRUE)
  }
}
immu <- merge(immu,arr,by="hhID",all.x=T)

immu$ph1.pre.c <-1*(immu$prevax.pH1>=40);immu$ph1.postv.c <-1*(immu$postvax.pH1>=40)
immu$ph1.rise <-immu$postvax.pH1/immu$prevax.pH1
immu$sh1.pre.c <- 1*(immu$prevax.sH1>=40);immu$sh1.postv.c <-1*(immu$postvax.sH1>=40)
immu$sh1.rise <-immu$postvax.sH1/immu$prevax.sH1
immu$sh3.pre.c <-1*(immu$prevax.sH3>=40);immu$sh3.postv.c <-1*(immu$postvax.sH3>=40)
immu$sh3.rise <-immu$postvax.sH3/immu$prevax.sH3
immu$b.pre.c <-1*(immu$prevax.B.Brisbane>=40);immu$b.postv.c <-1*(immu$postvax.B.Brisbane>=40)
immu$b.rise <-immu$postvax.B.Brisbane/immu$prevax.B.Brisbane

set.seed(12347)
immu.i <- aregImpute( ~ factor(intervention) +
                         sh1.pre.c + sh3.pre.c + sh1.postv.c + sh3.postv.c + ph1.pre.c + ph1.postv.c + I(sh1.rise) + I(sh3.rise) + I(ph1.rise) +
                         b.pre.c + b.postv.c + I(b.rise) +
                         I(fever)+I(chill)+I(tired)+I(headache)+I(cough)+I(mpain)+I(swell)+I(red)+I(bruise)+I(pain)+
                         male + I(agegp) + vac0809, data=immu, n.impute=10)

immu.nomiss <- list(immu, immu, immu, immu, immu, immu, immu, immu, immu, immu)

for(i in 1:10){
    immu.nomiss[[i]]$b.pre.c[is.na(immu.nomiss[[i]]$b.pre.c)] <- immu.i$imputed$b.pre.c[,i]
    immu.nomiss[[i]]$b.postv.c[is.na(immu.nomiss[[i]]$b.postv.c)] <- immu.i$imputed$b.postv.c[,i]
    immu.nomiss[[i]]$b.rise[is.na(immu.nomiss[[i]]$b.rise)] <- immu.i$imputed$b.rise[,i]
    immu.nomiss[[i]]$sh1.pre.c[is.na(immu.nomiss[[i]]$sh1.pre.c)] <- immu.i$imputed$sh1.pre.c[,i]
    immu.nomiss[[i]]$sh1.postv.c[is.na(immu.nomiss[[i]]$sh1.postv.c)] <- immu.i$imputed$sh1.postv.c[,i]
    immu.nomiss[[i]]$sh1.rise[is.na(immu.nomiss[[i]]$sh1.rise)] <- immu.i$imputed$sh1.rise[,i]
    immu.nomiss[[i]]$sh3.pre.c[is.na(immu.nomiss[[i]]$sh3.pre.c)] <- immu.i$imputed$sh3.pre.c[,i]
    immu.nomiss[[i]]$sh3.postv.c[is.na(immu.nomiss[[i]]$sh3.postv.c)] <- immu.i$imputed$sh3.postv.c[,i]
    immu.nomiss[[i]]$sh3.rise[is.na(immu.nomiss[[i]]$sh3.rise)] <- immu.i$imputed$sh3.rise[,i]
    immu.nomiss[[i]]$ph1.pre.c[is.na(immu.nomiss[[i]]$ph1.pre.c)] <- immu.i$imputed$ph1.pre.c[,i]
    immu.nomiss[[i]]$ph1.postv.c[is.na(immu.nomiss[[i]]$ph1.postv.c)] <- immu.i$imputed$ph1.postv.c[,i]
    immu.nomiss[[i]]$ph1.rise[is.na(immu.nomiss[[i]]$ph1.rise)] <- immu.i$imputed$ph1.rise[,i]
    immu.nomiss[[i]]$fever[is.na(immu.nomiss[[i]]$fever)] <- immu.i$imputed$fever[,i]
    immu.nomiss[[i]]$chill[is.na(immu.nomiss[[i]]$chill)] <- immu.i$imputed$chill[,i]
    immu.nomiss[[i]]$tired[is.na(immu.nomiss[[i]]$tired)] <- immu.i$imputed$tired[,i]
    immu.nomiss[[i]]$headache[is.na(immu.nomiss[[i]]$headache)] <- immu.i$imputed$headache[,i]
    immu.nomiss[[i]]$cough[is.na(immu.nomiss[[i]]$cough)] <- immu.i$imputed$cough[,i]
    immu.nomiss[[i]]$mpain[is.na(immu.nomiss[[i]]$mpain)] <- immu.i$imputed$mpain[,i]
    immu.nomiss[[i]]$swell[is.na(immu.nomiss[[i]]$swell)] <- immu.i$imputed$swell[,i]
    immu.nomiss[[i]]$red[is.na(immu.nomiss[[i]]$red)] <- immu.i$imputed$red[,i]
    immu.nomiss[[i]]$bruise[is.na(immu.nomiss[[i]]$bruise)] <- immu.i$imputed$bruise[,i]
    immu.nomiss[[i]]$pain[is.na(immu.nomiss[[i]]$pain)] <- immu.i$imputed$pain[,i]
    immu.nomiss[[i]]$vac0809[is.na(immu.nomiss[[i]]$vac0809)] <- immu.i$imputed$vac0809[,i]
}

# function to combined the imputed results
combine.mi <- function(data,m){       # data: hhID, group, event
  nTIV <- sum(data[[1]]$intervention=="TIV")
  npla <- sum(data[[1]]$intervention=="placebo")

  w <- rep(NA,m)
  ptab <- array(NA,c(length(unique(data[[1]]$event)),length(unique(data[[1]]$intervention)),m))
  for (i in 1:m){
    ctab <- table(data[[i]]$event,data[[i]]$intervention=="placebo")
    ptab[,,i] <- prop.table(ctab,2)
    w[i] <- chisq.test(ctab)$statistic
  }

  k <- (nrow(ctab)-1)*(ncol(ctab)-1)
  r2 <- (m+1)*sum((sqrt(w)-sqrt(mean(w)))^2)/(m*(m-1))
  W2 <- (mean(w)/k-(m+1)*r2/(m-1))/(1+r2)
  v2 <- k^(-3/m)*(m-1)*(1+1/r2)^2     # degree of freedom (combined chi-square test)
  output.per <- round((ptab[,,1]+ptab[,,2]+ptab[,,3]+ptab[,,4]+ptab[,,5]+ptab[,,6]+ptab[,,7]+ptab[,,8]+ptab[,,9]+ptab[,,10])/m,2)
	rownames(output.per) <- rownames(ctab)

  output <- as.data.frame(output.per)
  output <- output[rownames(output)!=0,]
  colnames(output) <- c("TIV%","placebo%")
	output$p.value = round(pf(W2,df1=k,df2=v2,lower.tail=FALSE),2)
	output
}

combine2.mi <- function(e1m,e1sd,e0m,e0sd,m){
  nTIV <- sum(immu.nomiss[[1]]$intervention=="TIV")
  npla <- sum(immu.nomiss[[1]]$intervention=="placebo")
	mean1 <- mean(e1m)
	mean0 <- mean(e0m)
	meanQ <- mean(e1m-e0m)
	Ubar <- mean(e1sd^2/nTIV+e0sd^2/npla)
	B <- sum((e1m-e0m-meanQ)^2)/(m-1)
	T <- (1 + 1/m) * B + Ubar
	degf <- (m-1)*(1+Ubar/((1+1/m)*B))^2
	data.frame( TIV = round(exp(mean1),1), placebo = round(exp(mean0),1),
              p.value = round(2*(1 - pt(abs(meanQ)/sqrt(T), df=degf)),2) )
}

# output

# sH1N1
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$sh1.pre.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[1,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$sh1.postv.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[2,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
event.TIV.mean <- event.TIV.sd <- event.placebo.mean <- event.placebo.sd <- rep(NA,10)
for (i in 1:10){
    temp[[i]]$sh1.rise <- log(temp[[i]]$sh1.rise)
    event.TIV.mean[i] <- mean(temp[[i]]$sh1.rise[temp[[i]]$intervention=="TIV"])
    event.TIV.sd[i] <- sd(temp[[i]]$sh1.rise[temp[[i]]$intervention=="TIV"])
    event.placebo.mean[i] <- mean(temp[[i]]$sh1.rise[temp[[i]]$intervention=="placebo"])
    event.placebo.sd[i] <- sd(temp[[i]]$sh1.rise[temp[[i]]$intervention=="placebo"])
}
tab[3,] <- as.numeric(combine2.mi(event.TIV.mean,event.TIV.sd,event.placebo.mean,event.placebo.sd,10))

# sH3N2
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$sh3.pre.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[4,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$sh3.postv.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[5,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
event.TIV.mean <- event.TIV.sd <- event.placebo.mean <- event.placebo.sd <- rep(NA,10)
for (i in 1:10){
    temp[[i]]$sh3.rise <- log(temp[[i]]$sh3.rise)
    event.TIV.mean[i] <- mean(temp[[i]]$sh3.rise[temp[[i]]$intervention=="TIV"])
    event.TIV.sd[i] <- sd(temp[[i]]$sh3.rise[temp[[i]]$intervention=="TIV"])
    event.placebo.mean[i] <- mean(temp[[i]]$sh3.rise[temp[[i]]$intervention=="placebo"])
    event.placebo.sd[i] <- sd(temp[[i]]$sh3.rise[temp[[i]]$intervention=="placebo"])
}
tab[6,] <- as.numeric(combine2.mi(event.TIV.mean,event.TIV.sd,event.placebo.mean,event.placebo.sd,10))

# Flu B Brisbane
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$b.pre.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[7,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$b.postv.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[8,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
event.TIV.mean <- event.TIV.sd <- event.placebo.mean <- event.placebo.sd <- rep(NA,10)
for (i in 1:10){
    temp[[i]]$b.rise <- log(temp[[i]]$b.rise)
    event.TIV.mean[i] <- mean(temp[[i]]$b.rise[temp[[i]]$intervention=="TIV"])
    event.TIV.sd[i] <- sd(temp[[i]]$b.rise[temp[[i]]$intervention=="TIV"])
    event.placebo.mean[i] <- mean(temp[[i]]$b.rise[temp[[i]]$intervention=="placebo"])
    event.placebo.sd[i] <- sd(temp[[i]]$b.rise[temp[[i]]$intervention=="placebo"])
}
tab[9,] <- as.numeric(combine2.mi(event.TIV.mean,event.TIV.sd,event.placebo.mean,event.placebo.sd,10))

# ph1N1
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$ph1.pre.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[10,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
for (i in 1:10){
    temp[[i]]$event <- 1*(temp[[i]]$ph1.postv.c)
    temp[[i]] <- temp[[i]][c("hhID","intervention","event")]
}
tab[11,] <- as.numeric(combine.mi(temp,10))
#
temp <- immu.nomiss
event.TIV.mean <- event.TIV.sd <- event.placebo.mean <- event.placebo.sd <- rep(NA,10)
for (i in 1:10){
    temp[[i]]$ph1.rise <- log(temp[[i]]$ph1.rise)
    event.TIV.mean[i] <- mean(temp[[i]]$ph1.rise[temp[[i]]$intervention=="TIV"])
    event.TIV.sd[i] <- sd(temp[[i]]$ph1.rise[temp[[i]]$intervention=="TIV"])
    event.placebo.mean[i] <- mean(temp[[i]]$ph1.rise[temp[[i]]$intervention=="placebo"])
    event.placebo.sd[i] <- sd(temp[[i]]$ph1.rise[temp[[i]]$intervention=="placebo"])
}
tab[12,] <- as.numeric(combine2.mi(event.TIV.mean,event.TIV.sd,event.placebo.mean,event.placebo.sd,10))

tab

#
# End of script.
#

