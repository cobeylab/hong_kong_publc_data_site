#
# R syntax to reproduce information for Table 2 from:
#
# Zhou Y, Ng DMW, Seto WH, et al.
# Seroprevalence of antibody to pandemic influenza A (H1N1) 2009 among
# healthcare works after the first wave in Hong Kong.
# Journal of Hospital Infection, 2011.
#
# Last updated by Zhou Y and Fang VJ.
# September 20, 2011

dir <- "http://sph.hku.hk/data/"

require(Hmisc)
require(Design)
data <- read.csv(paste(dir,"hcws2010.csv",sep=""))
tab2 <- matrix(NA,ncol=6,nrow=11)

unvax <- data[data$vacpH1N1==0,]
unvax$agegp <- as.character(factor(cut(unvax$age,c(0,24,34,44,54,100)),labels=c(2,1,3,4,5)))
unvax$hospital <- substr(unvax$Record_ID,1,1)

# MULTIPLE IMPUTATION for missing data

set.seed(123456)
hhw.i <- aregImpute( ~ VNantibody+male+I(agegp)+occupation+Numpatients+department+I(Numchildren)+hospital+vac0910, data=unvax, n.impute=10)
hh.nomiss <- list(unvax, unvax, unvax, unvax, unvax, unvax, unvax, unvax, unvax, unvax)
for(i in 1:10){
    hh.nomiss[[i]]$VNantibody[is.na(hh.nomiss[[i]]$VNantibody)] <- hhw.i$imputed$VNantibody[,i]
    hh.nomiss[[i]]$q3[is.na(hh.nomiss[[i]]$male)] <- hhw.i$imputed$male[,i]
    hh.nomiss[[i]]$agegp[is.na(hh.nomiss[[i]]$agegp)] <- hhw.i$imputed$agegp[,i]
    hh.nomiss[[i]]$q2[is.na(hh.nomiss[[i]]$occupation)] <- hhw.i$imputed$occupation[,i]
    hh.nomiss[[i]]$Numpatients[is.na(hh.nomiss[[i]]$Numpatients)] <- hhw.i$imputed$Numpatients[,i]
    hh.nomiss[[i]]$department[is.na(hh.nomiss[[i]]$department)] <- hhw.i$imputed$department[,i]
    hh.nomiss[[i]]$Numchildren[is.na(hh.nomiss[[i]]$Numchildren)] <- hhw.i$imputed$Numchildren[,i]
    hh.nomiss[[i]]$hospital[is.na(hh.nomiss[[i]]$hospital)] <- hhw.i$imputed$hospital[,i]
    hh.nomiss[[i]]$q10a[is.na(hh.nomiss[[i]]$vac0910)] <- hhw.i$imputed$vac0910[,i]
}

combine.mi <- function(model, n.impute){
	betas <- matrix(c(model[[1]]$coef, model[[2]]$coef, model[[3]]$coef, model[[4]]$coef, model[[5]]$coef,
	                	model[[6]]$coef, model[[7]]$coef, model[[8]]$coef, model[[9]]$coef, model[[10]]$coef), byrow=FALSE, ncol=10)
	vars <- matrix(c(diag(model[[1]][[6]]), diag(model[[2]][[6]]), diag(model[[3]][[6]]), diag(model[[4]][[6]]),
                   diag(model[[5]][[6]]), diag(model[[6]][[6]]), diag(model[[7]][[6]]), diag(model[[8]][[6]]),
		               diag(model[[9]][[6]]), diag(model[[10]][[6]])), byrow=FALSE, ncol=10)
	coef.names <- names(model[[1]]$coef)
	mean.coefs <- rowMeans(betas)
	Ubar <- rowMeans(vars)
	B <- rowSums((betas - mean.coefs)*(betas-mean.coefs) /
		(n.impute - 1))
	T <- (1 + 1/n.impute) * B + Ubar
	degf <- (n.impute - 1)*(1 + Ubar / ((1 + 1/n.impute)*B))*
		(1 + Ubar / ((1 + 1/n.impute)*B))
	output <- data.frame(OR = exp(mean.coefs),
		                   lowerCI = exp(mean.coefs - qt(0.975, df=degf)*sqrt(T)),
	                     upperCI = exp(mean.coefs + qt(0.975, df=degf)*sqrt(T)),
	                     p.value = 2*pt(-1*abs(mean.coefs/sqrt(T)),df=degf),
	                     row.names=coef.names)
  round(output,2)
}

# univariate model
risk.fit <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
risk.fit[[i]] <- lrm(VNantibody~factor(agegp),data=hh.nomiss[[i]])
}
tab2[c(1,3:5),1:3] <- as.matrix(combine.mi(risk.fit,10))[-1,1:3]

risk.fit <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
risk.fit[[i]] <- lrm(VNantibody~factor(department),data=hh.nomiss[[i]])
}
tab2[7:11,1:3] <- as.matrix(combine.mi(risk.fit,10))[-1,1:3]

# multivariate model
risk.fit <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
risk.fit[[i]] <- lrm(VNantibody~factor(agegp)+factor(department),data=hh.nomiss[[i]])
}
tab2[c(1,3:5,7:11),4:6] <- as.matrix(combine.mi(risk.fit,10))[-1,1:3]

tab2[c(2,6),c(1,4)] <- 1
colnames(tab2) <- c("Crude OR","CI_low","CI_up","Adjusted OR","CI_low","CI_up")
rownames(tab2) <- c("Age:19-24","25-34","35-44","45-54","55-64","Department:Medicine","Surgery","Emergency","Paediatrics","Other clinical","Non-clinical")
tab2

# End of script.
