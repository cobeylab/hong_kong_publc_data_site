#
# R syntax to reproduce information for Table 7 from:
#
# Cowling BJ, Ng S, Ma ESK, et al.
# Protective efficacy of seasonal influenza vaccination against seasonal and
# pandemic influenza virus infection during 2009 in Hong Kong
# CID, 2010 (in press).
#
# Last updated by Fang VJ, Ng S and Cowling BJ.
# October 21, 2010

source("http://www.hku.hk/bcowling/influenza/kiddivaxPilot_scripts/dataframe_cross.r")

require(Hmisc)

tab7.index <- matrix(rep(NA, 4*10), nrow=10,
               dimnames=list(c("Age6-8", "Age9-15", "Female", "Male", "No flu A", "Flu A",
                                "no vax", "seasonal vax", "Before Oct 1", "After Oct 1"),
                              c("n",  "OR",  "CI_low", "CI_up")))

tab7.member <- matrix(rep(NA, 4*11), nrow=11,
               dimnames=list(c("Age<16", "Age16-45", "Age>45", "Female", "Male", "No flu A", "Flu A",
                                "not contact of  svax", "contact of  svax", "Before Oct 1", "After Oct 1"),
                              c("n",  "OR",  "CI_low", "CI_up")))

##################################################
#Multiple Imputation for missing index info      #
##################################################
index.pre$agegp<-cut(index.pre$age,c(0,8,15))

set.seed(2)
index.i <- aregImpute( ~ factor(intervention) +
                         sh1w + sh1s + sh3w + sh3s + ph1s + ariw + iliw + aris + ilis + ilip + I(day.p) +
                          flo + flo.w + flo.s + male+ I(agegp)+ vac09, data=index.pre,n.impute=10)


index.nomiss <- list(index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre, index.pre)

for(i in 1:10){
    index.nomiss[[i]]$flo[is.na(index.nomiss[[i]]$flo)] <- index.i$imputed$flo[,i]
    index.nomiss[[i]]$flo.w[is.na(index.nomiss[[i]]$flo.w)] <- index.i$imputed$flo.w[,i]
    index.nomiss[[i]]$flo.s[is.na(index.nomiss[[i]]$flo.s)] <- index.i$imputed$flo.s[,i]
    index.nomiss[[i]]$sh1w[is.na(index.nomiss[[i]]$sh1w)] <- index.i$imputed$sh1w[,i]
    index.nomiss[[i]]$sh3w[is.na(index.nomiss[[i]]$sh3w)] <- index.i$imputed$sh3w[,i]
    index.nomiss[[i]]$sh3s[is.na(index.nomiss[[i]]$sh3s)] <- index.i$imputed$sh3s[,i]
    index.nomiss[[i]]$ph1s[is.na(index.nomiss[[i]]$ph1s)] <- index.i$imputed$ph1s[,i]
    index.nomiss[[i]]$sh1s[is.na(index.nomiss[[i]]$sh1s)] <- index.i$imputed$sh1s[,i]
    index.nomiss[[i]]$iliw[is.na(index.nomiss[[i]]$iliw)] <- index.i$imputed$iliw[,i]
    index.nomiss[[i]]$ariw[is.na(index.nomiss[[i]]$ariw)] <- index.i$imputed$ariw[,i]
    index.nomiss[[i]]$ilis[is.na(index.nomiss[[i]]$ilis)] <- index.i$imputed$ilis[,i]
    index.nomiss[[i]]$aris[is.na(index.nomiss[[i]]$aris)] <- index.i$imputed$aris[,i]
    index.nomiss[[i]]$ilip[is.na(index.nomiss[[i]]$ilip)] <- index.i$imputed$ilip[,i]
    index.nomiss[[i]]$day.p[is.na(index.nomiss[[i]]$day.p)] <- index.i$imputed$day.p[,i]
    }

###################################################
#Multiple Imputation for missing member info      #
###################################################
# Re-do agegp for regression
member.pre$agegp<-cut(member.pre$age,c(0,15,45,100))

set.seed(1)
member.i <- aregImpute( ~ factor(intervention)+
                         sh1w + sh1s + sh3w + sh3s + ph1s + ariw + iliw + aris + ilis + ilip + I(day.p) +
                         flo + flo.w + flo.s + male+ I(agegp)+ vac09, data=member.pre,n.impute=10)

member.nomiss <- list(member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre, member.pre)

for(i in 1:10){
    member.nomiss[[i]]$flo.s[is.na(member.nomiss[[i]]$flo.s)] <- member.i$imputed$flo.s[,i]
    member.nomiss[[i]]$flo.w[is.na(member.nomiss[[i]]$flo.w)] <- member.i$imputed$flo.w[,i]
    member.nomiss[[i]]$flo[is.na(member.nomiss[[i]]$flo)] <- member.i$imputed$flo[,i]
    member.nomiss[[i]]$sh1w[is.na(member.nomiss[[i]]$sh1w)] <- member.i$imputed$sh1w[,i]
    member.nomiss[[i]]$sh3w[is.na(member.nomiss[[i]]$sh3w)] <- member.i$imputed$sh3w[,i]
    member.nomiss[[i]]$sh3s[is.na(member.nomiss[[i]]$sh3s)] <- member.i$imputed$sh3s[,i]
    member.nomiss[[i]]$ph1s[is.na(member.nomiss[[i]]$ph1s)] <- member.i$imputed$ph1s[,i]
    member.nomiss[[i]]$sh1s[is.na(member.nomiss[[i]]$sh1s)] <- member.i$imputed$sh1s[,i]
    member.nomiss[[i]]$iliw[is.na(member.nomiss[[i]]$iliw)] <- member.i$imputed$iliw[,i]
    member.nomiss[[i]]$ariw[is.na(member.nomiss[[i]]$ariw)] <- member.i$imputed$ariw[,i]
    member.nomiss[[i]]$ilis[is.na(member.nomiss[[i]]$ilis)] <- member.i$imputed$ilis[,i]
    member.nomiss[[i]]$aris[is.na(member.nomiss[[i]]$aris)] <- member.i$imputed$aris[,i]
    member.nomiss[[i]]$ilip[is.na(member.nomiss[[i]]$ilip)] <- member.i$imputed$ilip[,i]
    member.nomiss[[i]]$vac09[is.na(member.nomiss[[i]]$vac09)] <- member.i$imputed$vac09[,i]
    member.nomiss[[i]]$day.p[is.na(member.nomiss[[i]]$day.p)] <- member.i$imputed$day.p[,i]
    }

# Logistic regression

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
	                     row.names=coef.names)
  round(output,2)
}


require(Design)
risk.fit <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
risk.fit[[i]] <- lrm(ph1s~factor(agegp)+male+(sh1w==1|sh1s==1|sh3w==1|sh3s==1)+vac09+cut(day.p,c(0,30,60)),data=index.nomiss[[i]])
}
out.tab7.index <- combine.mi(risk.fit,10)

risk.fit <- list(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
for (i in 1:10){
risk.fit[[i]] <- lrm(ph1s~factor(agegp)+male+(sh1w==1|sh1s==1|sh3w==1|sh3s==1)+intervention+cut(day.p,c(0,30,60)),data=member.nomiss[[i]])
}
out.tab7.member <- combine.mi(risk.fit,10)

tab7.index[c(1,3,5,7,9),2] <- 1.00
tab7.index[,1] <- c(table(index.pre$agegp),table(index.pre$male),table(index.pre$sh1w==1|index.pre$sh1s==1|index.pre$sh3w==1|index.pre$sh3s==1),
                      table(index.pre$vac09),table(cut(index.pre$day.p,c(0,30,60))))
tab7.index[c(2,4,6,8,10),2:4] <- as.matrix(out.tab7.index[-1,])
tab7.index

tab7.member[c(1,4,6,8,10),2] <- 1.00
tab7.member[,1] <- c(table(member.pre$agegp),table(member.pre$male),table(member.pre$sh1w==1|member.pre$sh1s==1|member.pre$sh3w==1|member.pre$sh3s==1),
                      table(member.pre$intervention),table(cut(member.pre$day.p,c(0,30,60))))
tab7.member[c(2,3,5,7,9,11),2:4] <- as.matrix(out.tab7.member[-1,])
tab7.member

#
# End of script.
#
