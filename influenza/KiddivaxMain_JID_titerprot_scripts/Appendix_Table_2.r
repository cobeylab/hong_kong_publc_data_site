#
# R syntax to reproduce Appendix Table 2 from:
#
# Ng S, Fang VJ, Ip DKM, et al.
# Estimation of the association between antibody titers and
# protection against confirmed influenza virus infection in children
# JID, 2013.
#
# Last updated by Fang VJ, Ng Sophia, and Cowling BJ.
# December 2014

source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_2.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_5.r")
source("http://www.hku.hk/bcowling/influenza/KiddivaxMain_JID_titerprot_scripts/source_6.r")

#function to convert z-value to p-value by putting to the model summary
getp <- function (model) {
                      x <-data.frame(summary(model)[7])[4,5]
                      if (x<0)
                      {p <-1- pnorm(abs(x), mean=0, sd=1)}
                      else
                      {p <- pnorm(x, mean=0, sd=1)}
                      
                      round(p,2)}

# Create tables for waning rates
# Main results
b1.y <-cbind(model.b1.age$coefficients[2]+model.b1.age$coefficients[4] ,
          model.b1.age$coefficients[2]+model.b1.age$coefficients[4] + 1.96*sqrt(sum(model.b1.age$robust.variance[c(2,4),c(2,4)])) ,
          model.b1.age$coefficients[2]+model.b1.age$coefficients[4] - 1.96*sqrt(sum(model.b1.age$robust.variance[c(2,4),c(2,4)])))

b1.o <-cbind(model.b1.age$coefficients[2] ,
          model.b1.age$coefficients[2] + 1.96*sqrt(sum(model.b1.age$robust.variance[2,2])) ,
          model.b1.age$coefficients[2] - 1.96*sqrt(sum(model.b1.age$robust.variance[2,2])))

b0 <-cbind(model.b0$coefficients[2] ,
          model.b0$coefficients[2] + 1.96*sqrt(model.b0$robust.variance[2,2]) ,
           model.b0$coefficients[2] - 1.96*sqrt(model.b0$robust.variance[2,2])
           )

p1 <-cbind(model.p1$coefficients[2] ,
          model.p1$coefficients[2] + 1.96*sqrt(model.p1$robust.variance[2,2]) ,
          model.p1$coefficients[2] - 1.96*sqrt(model.p1$robust.variance[2,2]))

p0 <-cbind(model.p0$coefficients[2] ,
          model.p0$coefficients[2] + 1.96*sqrt(model.p0$robust.variance[2,2]) ,
          model.p0$coefficients[2] - 1.96*sqrt(model.p0$robust.variance[2,2]),
          getp(model.pi.m))

main <-round(rbind( cbind(n.m[1],p1,NA,n.m[2],p0),
                    cbind(NA,NA,NA,NA,NA,n.m[5],b0,NA),
                    cbind(n.m[3],b1.y,NA,NA,NA,NA,NA,NA),
                    cbind(n.m[4],b1.o,NA,NA,NA,NA,NA,NA)),
                    4)

1/abs(rbind( cbind(p1,p0),
         cbind(NA,NA,NA,b0,NA),
         cbind(b1.y,NA,NA,NA,NA),
         cbind(b1.o,NA,NA,NA,NA)))    # half-life of HI antibody

# Sen 50% results
b1.30.y <-cbind(model.b1.30.age$coefficients[2]+model.b1.30.age$coefficients[4] ,
          model.b1.30.age$coefficients[2]+model.b1.30.age$coefficients[4] + 1.96*sqrt(sum(model.b1.30.age$robust.variance[c(2,4),c(2,4)])) ,
          model.b1.30.age$coefficients[2]+model.b1.30.age$coefficients[4] - 1.96*sqrt(sum(model.b1.30.age$robust.variance[c(2,4),c(2,4)])))

b1.30.o <-cbind(model.b1.30.age$coefficients[2] ,
          model.b1.30.age$coefficients[2] + 1.96*sqrt(sum(model.b1.30.age$robust.variance[2,2])) ,
          model.b1.30.age$coefficients[2] - 1.96*sqrt(sum(model.b1.30.age$robust.variance[2,2])))

b0.30 <-cbind(model.b0.30$coefficients[2] ,
          model.b0.30$coefficients[2] + 1.96*sqrt(model.b0.30$robust.variance[2,2]) ,
           model.b0.30$coefficients[2] - 1.96*sqrt(model.b0.30$robust.variance[2,2])
           )

p1.30 <-cbind(model.p1.30$coefficients[2] ,
          model.p1.30$coefficients[2] + 1.96*sqrt(model.p1.30$robust.variance[2,2]) ,
          model.p1.30$coefficients[2] - 1.96*sqrt(model.p1.30$robust.variance[2,2]))

p0.30 <-cbind(model.p0.30$coefficients[2] ,
          model.p0.30$coefficients[2] + 1.96*sqrt(model.p0.30$robust.variance[2,2]) ,
          model.p0.30$coefficients[2] - 1.96*sqrt(model.p0.30$robust.variance[2,2]),
          getp(model.pi.m))

sen30 <-round(rbind( cbind(n30[1],p1.30,NA,n30[2],p0.30),
                    cbind(NA,NA,NA,NA,NA,n30[5],b0.30,NA),
                    cbind(n30[3],b1.30.y,NA,NA,NA,NA,NA,NA),
                    cbind(n30[4],b1.30.o,NA,NA,NA,NA,NA,NA)),
                    4)

1/abs(rbind( cbind(p1.30,p0.30),
         cbind(NA,NA,NA,b0.30,NA),
         cbind(b1.30.y,NA,NA,NA,NA),
         cbind(b1.30.o,NA,NA,NA,NA)))

# Sen 5xpcr results
b1.xpcr.y <-cbind(model.b1.xpcr.age$coefficients[2]+model.b1.xpcr.age$coefficients[4] ,
          model.b1.xpcr.age$coefficients[2]+model.b1.xpcr.age$coefficients[4] + 1.96*sqrt(sum(model.b1.xpcr.age$robust.variance[c(2,4),c(2,4)])) ,
          model.b1.xpcr.age$coefficients[2]+model.b1.xpcr.age$coefficients[4] - 1.96*sqrt(sum(model.b1.xpcr.age$robust.variance[c(2,4),c(2,4)])))

b1.xpcr.o <-cbind(model.b1.xpcr.age$coefficients[2] ,
          model.b1.xpcr.age$coefficients[2] + 1.96*sqrt(sum(model.b1.xpcr.age$robust.variance[2,2])) ,
          model.b1.xpcr.age$coefficients[2] - 1.96*sqrt(sum(model.b1.xpcr.age$robust.variance[2,2])))

b0.xpcr <-cbind(model.b0.xpcr$coefficients[2] ,
          model.b0.xpcr$coefficients[2] + 1.96*sqrt(model.b0.xpcr$robust.variance[2,2]) ,
           model.b0.xpcr$coefficients[2] - 1.96*sqrt(model.b0.xpcr$robust.variance[2,2])
           )

p1.xpcr <-cbind(model.p1.xpcr$coefficients[2] ,
          model.p1.xpcr$coefficients[2] + 1.96*sqrt(model.p1.xpcr$robust.variance[2,2]) ,
          model.p1.xpcr$coefficients[2] - 1.96*sqrt(model.p1.xpcr$robust.variance[2,2]))

p0.xpcr <-cbind(model.p0.xpcr$coefficients[2] ,
          model.p0.xpcr$coefficients[2] + 1.96*sqrt(model.p0.xpcr$robust.variance[2,2]) ,
          model.p0.xpcr$coefficients[2] - 1.96*sqrt(model.p0.xpcr$robust.variance[2,2]),
          getp(model.pi.m))

sen.xpcr <-round(rbind( cbind(n.xpcr[1],p1.xpcr,NA,n.xpcr[2],p0.xpcr),
                    cbind(NA,NA,NA,NA,NA,n.xpcr[5],b0.xpcr,NA),
                    cbind(n.xpcr[3],b1.xpcr.y,NA,NA,NA,NA,NA,NA),
                    cbind(n.xpcr[4],b1.xpcr.o,NA,NA,NA,NA,NA,NA)),
                    4)

1/abs(rbind( cbind(p1.xpcr,p0.xpcr),
         cbind(NA,NA,NA,b0.xpcr,NA),
         cbind(b1.xpcr.y,NA,NA,NA,NA),
         cbind(b1.xpcr.o,NA,NA,NA,NA)))
         
AT2 <-rbind(main, NA,NA, sen30, NA,NA,sen.xpcr)

# Convert daily change to monthly change
AT3 <- AT2

AT3[,5] <- ((AT3[,3]-AT3[,2])/1.96)^2
AT3[,2] <- AT3[,2]*30
AT3[,3] <-AT3[,2] + 1.96 * sqrt(AT3[,5] * 30^2)
AT3[,4] <-AT3[,2] - 1.96 * sqrt(AT3[,5] * 30^2)

AT3[,5] <- ((AT3[,8]-AT3[,7])/1.96)^2
AT3[,7] <- AT3[,7]*30
AT3[,8] <-AT3[,7] + 1.96 * sqrt(AT3[,5] * 30^2)
AT3[,9] <-AT3[,7] - 1.96 * sqrt(AT3[,5] * 30^2)

AT3[,5] <- NA
AT3

#
# End of script.
#



